# -----------------------------------------------------------------------------
# Description:  Various functions relating calculation of CPL-inequalities,
#               generation of the secondary cone based off of 2-face
#               triangulation data, and generating FRSTs from said data.
# -----------------------------------------------------------------------------

# standard imports
import multiprocessing as mp
import numpy as np
import os
import pickle
import qpsolvers
import random
from scipy import sparse
try:
    import sympy as smp
except:
    import subprocess
    import sys
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'sympy'])
    import sympy as smp

import collections, itertools, math, warnings, time

import atexit
from tqdm import tqdm

# CYTools imports
from cytools import config
from cytools.cone import Cone, feasibility
from cytools.polytope import Polytope
from cytools.polytopeface import PolytopeFace
from cytools.triangulation import Triangulation

# repo imports
try:
    from .matrix import *
except:
    from matrix import *

# face triangulation methods
# --------------------------
def compute_face_triangs(poly: Polytope,
                 dim: int=2,
                 which: list=None,
                 only_regular: bool=True,
                 max_npts: int=None,
                 N_face_triangs: int=1000,
                 triang_method: str="grow2d",
                 seed: int=None,
                 verbosity=0):
    """
    **Description:**
    For each dim-face of a given polytope, generate the list of all fine
    dim-face triangulations of it.

    **Arguments:**
    - `poly`: Polytope of interset.
    - `dim`: The dimension of the faces to generate triangulations for.
    - `only_regular`: Whether to only generate the regular 2-face
        triangulations.
    - `max_npts`: The maximum number of points in a 2-face s.t. we try to get
        all triangulations.
    - `N_face_triangs`: If we aren't trying to get all triangulations, how
        many to grab per face (upper bound).
    - `triang_method`: The method to generate random triangulations. Allowed
        are "fast", "fair", and "grow2d".
    - `seed`: Random seed if grabbing only some triangulations.
    - `verbosity` *(int, optional)*: Verbosity level.

    **Returns:**
    *(list of lists of Triangulation)* List of faces. Each face is a list of
        Triangulation. May have different order than in the input (i.e.,
        faces.)
    *(list of dictionaries)* List of faces. Each face is dictionary from point
        to polytopal index
    """
    # output variable
    triangs = []

    # the faces
    if isinstance(which,int):
        which=[which]

    faces = poly.faces(dim)
    if which is not None:
        faces = [faces[i] for i in which]

    # iterate over faces
    ind = 0
    for face in faces:
        p = face.as_poly()      # convert to Polytope to get all triangulations

        if (max_npts is not None) and (len(p.points())>max_npts):
            if verbosity>=1:
                print(f"compute_face_triangs: 2face #{ind} has {len(p.points())} points! "
                           f"This is >={max_npts} (user-set limit). "
                           f"Will only request <={N_face_triangs} random "
                            "samples.")
            
            allowed_methods = ["fast","fair","grow2d"]
            if triang_method not in allowed_methods:
                raise ValueError(f"triang_method={triang_method} was not an "\
                                 f"allowed method... Allowed are {allowed_methods}.")
            elif triang_method=="fast":
                triangs.append(p.random_triangulations_fast(N=N_face_triangs,\
                                                as_list=True,make_star=False,\
                                                include_points_interior_to_facets=True,\
                                                seed=seed))
            elif triang_method=="fair":
                if verbosity>=1:
                    print("compute_face_triangs: warning... fair never worked well...")
                triangs.append(p.random_triangulations_fair(N=N_face_triangs,\
                                                as_list=True,make_star=False,\
                                                include_points_interior_to_facets=True,\
                                                seed=seed))
            elif triang_method=="grow2d":
                triangs.append(list(p.grow_frt(N=N_face_triangs,seed=seed)))

            if verbosity>=2:
                print(f"compute_face_triangs: found {len(triangs[-1])} triangs...")
                if len(triangs[-1]) < N_face_triangs:
                    print(f"compute_face_triangs: This is < {N_face_triangs}... "
                           "maybe you got them all? (not guaranteed)")
        else:
            if verbosity>=1:
                    print("compute_face_triangs: computing all face triangulations...")
            triangs.append(p.all_triangulations(only_fine=True, only_star=False,\
            only_regular=only_regular, include_points_interior_to_facets=True,
            as_list=True))
        
        # increment ind (only used for warning message)
        ind += 1

    return triangs
Polytope.compute_face_triangs = compute_face_triangs

def n_2face_triangs(p, only_regular=True):
    """
    **Description:**
    Return the count of all 2-face triangulations of the input polytope

    **Arguments:**
    - `p` *(Polytope)*: The polytope of interest

    **Returns:**
    *(integer)* The count of distinct sets of 2-face triangulations
    """
    triangs = p.compute_face_triangs(dim=2, only_regular=only_regular)
    return math.prod([len(f) for f in triangs])
Polytope.n_2face_triangs = n_2face_triangs
Polytope.num_2face_triangs = n_2face_triangs

# (low-level) 2-face inequality functions
# ---------------------------------------
# prefix with '_' to indicate that these shouldn't directly be called by user

# a large slowdown in _2d_frt_cone_ineqs is calculating nullspaces...
# cache them here...
file = 'cached_twoface_ineqs.p'
if os.path.isfile(file):
    with gzip.open(file, 'rb') as f:
        _ineq_cached = pickle.load(f)
else:
    _ineq_cached = dict()

def _save_ineqs():
    with gzip.open(file, 'wb') as f:
        pickle.dump(_ineq_cached, f)

atexit.register(lambda:_save_ineqs)

def _2d_frt_cone_ineqs(self, secondary_dim):
    """
    **Description:**
    Compute the secondary cone for a 2-face with input triangulation. This is
    the cone whose interior gives the height-vectors which would lead to said
    triangulation.

    **Overview:**
    The hyperplane inequalities/normals are calculated by looking at each pair
    of simplices that share an edge. For each pair, there are thus 4 relevant
    points, p0, p1, p2, and p3. Order the points such that p0 and p1 define
    the shared edge.

    The associated inequalities/normals are then calculated as the (basis
    vectors of the) null-space of the matrix
        M = [[p0_x, p1_x, p2_x, p3_x],
             [p0_y, p1_y, p2_y, p3_y],
             [   1,    1,    1,    1]],
    (where we fix the sign of these basis vectors, v, such that v[2]<0... This
    ensures that the hyperplane inequality 'points in the right way'). This matrix
    is the 'homogenization' of our points.

    **Explanation:**
    The hyperplane inequalities/normals can be thought of simply as ensuring
    that the height of the line segment (p0,p1) is lower than the height of the
    line segment (p2,p3) at the point of intersection.
    
    Why? Because the existence of the line segment (p0,p1) is equivalent to
    the fact that it is lower than the plane (p0,p2,p3) and the plane
    (p1,p2,p3)... that is, adding one of the endpoints LOWERS the convex hull,
    making new faces. To check this, we can check if the (interior of the)
    line segment is lower than the plane. The intersection point *of the two
    lines* is just a convenient point to do this comparison. Further, due to
    linearity, this logic works even if the intersection point is outside the
    polygon.

    This leads to a simple calculation: find (A,B) such that
        A*p0 + (1-A)*p1 == B*p2 + (1-B)*p3.
    This is the point of intersection. Then, the vector
        n = (A, 1-A, -B, -(1-B))
    defines our inequality simply:
        h.n > 0 <=> A*h0 + (1-A)*h1 + (-B)*h2 + (-1+B)*h3 > 0
                <=> A*h0 + (1-A)*h1 > B*h2 + (1-B)*h3
                <=> -A*h0 + (A-1)*h1 < -B*h2 + (B-1)*h3
    (We have to introduce the sign convention so that the inequality goes the
    right way.)

    It's convenient to reformulate this calculation in terms of matrices. To do
    this, we rewrite A=c0/(c0+c1) and B=c2/(c0+c1), giving
        n = (c0/(c0+c1), c1/(c0+c1), c2/(c0+c1), -(c0+c1+c2)/(c0+c1)).
    Note r*n defines the same triangulation (up to sign) for any scalar r, so
    we can scale
        n = (c0, c1, c2, -(c0+c1+c2))
    Thus, to find the intersection of the four points, all we need to do is
    find an n of this form, such that [p0, p1, p2, p3]@n = 0. This is exactly
    the nullspace calculation (the row of 1s ensures that c3=-c0-c1-c2...)

    **Arguments:**
    - `triang` *(Triangulation)*: The 2-face triangulation
    - `secondary_dim` *(int)*: The dimension of the secondary-cone space
        (i.e., the number of points in the polytope).

    **Returns:**
    *(LIL)* Each row is an inwards-facing hyperplane normal...
        represents a CPL inequality
    """
    # the output variable (doesn't need to be LIL object, but that is nice...)
    ineqs = LIL(dtype=np.int8, width=secondary_dim)

    # relevant inputs
    simps = self.simplices(as_labels=True)

    # for each element, find indices of all rows that include it
    ele_to_row = collections.defaultdict(list)
    for i, row in enumerate(simps):
        for elem in row:
            ele_to_row[elem].append(i)

    # for pair of rows/simplices, calculate the shared elements
    # (this contains quads, along with simplices sharing only 1 ele)
    row_to_shared = collections.defaultdict(set)
    for ele, indices in ele_to_row.items():
        for pair in itertools.combinations(indices, 2):
            row_to_shared[pair].add(ele)

    # Find pairs of rows that share at least two common elements
    for rows, s in row_to_shared.items():
        if len(s)<=1: continue
        else: s = list(s)

        simp1, simp2 = simps[rows[0]], simps[rows[1]]

        # calculate the not-shared points
        n_s = [lab for lab in list(simp1)+list(simp2) if (lab not in s)]

        # find the null-space of the following matrix
        M = self.points(which=n_s+s, optimal=True).T
        M_tup = tuple(tuple(row[1:]-row[0]) for row in M)

        # Grab/calculate the nullspace
        # (want exact answers, hence sympy...)
        ineq = _ineq_cached.get(M_tup,None)
        ineq = None
        if ineq is None:
            # calculate the nullspace
            null = smp.Matrix(M.tolist() + [[1,1,1,1]]).nullspace()[0]
            denom = [term.q for term in null]
            null *= math.lcm(*denom)

            # ensure the not-shared points have positive coordinates
            if null[0]<0: ineq = [-int(x) for x in null]
            else:         ineq = [int(x) for x in null]
            
            # cache this answer
            _ineq_cached[M_tup] = ineq

        # define the associated hyperplane normal
        ineqs.new_row()
        if ineq[0]!=0: ineqs[-1,n_s[0]] = ineq[0]
        if ineq[1]!=0: ineqs[-1,n_s[1]] = ineq[1]
        if ineq[2]!=0: ineqs[-1,  s[0]] = ineq[2]
        if ineq[3]!=0: ineqs[-1,  s[1]] = ineq[3]

    return ineqs
Triangulation._2d_frt_cone_ineqs = _2d_frt_cone_ineqs

def _2d_s_cone_ineqs(self, poly, secondary_dim):
    """
    Enforces star-ness...
    """

    # the output variable (doesn't need to be LIL object, but that is nice...)
    ineqs = LIL(dtype=np.int8, width=secondary_dim)

    # relevant inputs
    pts = self.points()

    # find each facet containing each 2d simplex
    containing_facets = collections.defaultdict(list)
    for s in self.simplices(2, as_labels=True):
        for f in poly.faces(3):
            if set(s).issubset(set(f.labels)):
                containing_facets[tuple(s)].append(f)

    # for each 2d simplex, enforce that it (with origin) appears for each
    # 4d circuit
    o = poly.label_origin
    for s in self.simplices(2, as_labels=True):
        s = s.tolist()
        for f1,f2 in itertools.combinations(containing_facets[tuple(s)],2):
            f1_only = set(f1.labels_bdry)-set(f2.labels_bdry)-set(s)
            f2_only = set(f2.labels_bdry)-set(f1.labels_bdry)-set(s)
            for p1,p2 in itertools.product(f1_only,f2_only):
                # calculate the not-shared points
                n_s = [p1,p2]

                # find the null-space of the following matrix
                M = poly.points(which=n_s+s+[o], optimal=True).T

                # Grab/calculate the nullspace
                # (want exact answers, hence sympy...)
                # calculate the nullspace
                null = smp.Matrix(M.tolist() + [[1,1,1,1,1,1]]).nullspace()[0]
                denom = [term.q for term in null]
                null *= math.lcm(*denom)

                # ensure the not-shared points have positive coordinates
                if null[0]<0: ineq = [-int(x) for x in null]
                else:         ineq = [int(x) for x in null]

                # define the associated hyperplane normal
                ineqs.new_row()
                if ineq[0]!=0: ineqs[-1,p1]   = ineq[0]
                if ineq[1]!=0: ineqs[-1,p2]   = ineq[1]
                if ineq[2]!=0: ineqs[-1,s[0]] = ineq[2]
                if ineq[3]!=0: ineqs[-1,s[1]] = ineq[3]
                if ineq[4]!=0: ineqs[-1,s[2]] = ineq[4]
                if ineq[5]!=0: ineqs[-1,o]    = ineq[5]

    return ineqs
Triangulation._2d_s_cone_ineqs = _2d_s_cone_ineqs

# generate secondary cone
# -----------------------
def permissible_heights(triangs: [Triangulation],
                        npts: int,
                        nworkers: int = 1,
                        no_duplicates: bool = True,
                        as_cone: bool = True,
                        require_star: bool = False) -> "LIL | Cone":
    """
    **Description:**
    See/cite https://arxiv.org/abs/2309.10855

    For an input set of ***2-face*** triangulations, generate the hyperplanes
    defining the 'expanded-secondary' cone. That is, the cone whose interior
    gives height vectors leading to the corresponding FRTs of its 2-faces.

    We call this 'expanded' since we don't control the origin.

    **Arguments:**
    - `triangs` The triangulations for each face.
    - `npts`: The number of points in the 4D polytope.
    - `nworkers`: Number of processors to use. If None or <=1, then runs
        sequentially
    - `no_duplicates`: Whether to only add unique ineqs.
    - `as_cone`: Whether to return a formal Cone object.

    **Returns:**
    The expanded secondary cone, either as hyperplanes or as a formal Cone
    object.
    """
    # the output variable (doesn't need to be LIL object, but that is nice...)
    ineqs = LIL(dtype=np.int8, width=npts, nworkers=nworkers)

    # iterate over face triangulations
    for f_triang in triangs:
        # CPL inequalities associated with ith triangulation
        # (normally, this is the triangulation of the ith face, but it doesn't
        # need to be... you can decide to pass a subset of faces)
        f_ineqs = _2d_frt_cone_ineqs(f_triang, npts)
        if require_star:
            f_ineqs.append(_2d_s_cone_ineqs(f_triang, p, npts))

        ineqs.append(f_ineqs, tocopy=False)    

    if no_duplicates:   ineqs.unique_rows(allow_shuffle=True)
    if as_cone:
        return hypers_to_cone(ineqs, dim=npts)
    else:
        return ineqs

def expanded_secondary_cone(self,
                            no_duplicates: bool = True,
                            as_cone: bool = True) -> "LIL | Cone":
    """
    **Description:**
    See/cite https://arxiv.org/abs/2309.10855

    For a 4D triangulation, generate the 'expanded-secondary' cone. That is,
    the cone whose interior gives height vectors leading to the corresponding
    FRTs of its 2-faces.

    We call this 'expanded' since we don't control the origin.

    **Arguments:**
    - `no_duplicates`: Whether to only add unique ineqs.
    - `as_cone`: Whether to return a formal Cone object.

    **Returns:**
    The expanded secondary cone, either as hyperplanes or as a formal Cone
    object.
    """
    return permissible_heights(self.restrict(dim=2,as_poly=True),
                               npts=len(self.points()),
                               no_duplicates=no_duplicates,
                               as_cone=as_cone)
Triangulation.expanded_secondary_cone = expanded_secondary_cone

# conversion functions
# --------------------
def hypers_to_frt(hypers,
                   p,
                   to_skip=None,
                   reduced_pts=None,
                   make_star=False,
                   backend=None,
                   verbosity=0):
    """
    **Description:**
    Convert (expanded) secondary cone hyperplanes to the corresponding FR(S)T.

    **Arguments:**
    - `hypers` *(array-like)*: Hyperplanes defining a expanded-secondary cone.
    - `p` *(Polytope)*: The corresponding polytope.
    - `to_skip` *((list of) ints, optional)*: The indices to skip in making the
        FRST. These are non-origin points not in any 2-face.
    - `reduced_pts` *(list of lists, optional)*: The coords of vertices in p,
        after removing the vertices listed in to_skip.
    - `backend` *(string, optional)*: The backend to use for cone calculations.
    - `verbosity` *(int, optional)*: Verbosity level.

    **Returns:**
    *(Triangulaton)* The FRST.
    """
    ambient_dim = len(p.points())

    # set the backend
    if backend is None:
        if config.mosek_is_activated() and (ambient_dim>=25):
            backend = "mosek"
        else:
            backend = "glop"

    # input checking
    if verbosity>=1:
        print("hypers_to_frt: Calculating 'to_skip' and 'reduced_pts'...")

    if to_skip is None:
        to_skip = p.labels_facet

    if reduced_pts is None:
        pts = p.points()
        reduced_pts = np.delete(pts, to_skip, 0)

    # calculate an interior point
    if backend=="mosek":
        # The problem is defined as:
        # Minimize (1/2) x.P.x + q.x
        # Subject to G.x <= h
        P = 2*sparse.identity(hypers.shape[1], dtype=float, format="csc")
        q = np.zeros(hypers.shape[1], dtype=float)
        h = np.full(hypers.shape[0], -1, dtype=float)
        G = -1*sparse.csc_matrix(hypers, dtype=float)
        heights = qpsolvers.solve_qp(P,q,G,h, solver="mosek", max_iter=10**6, verbose=verbosity>=2)
    else:
        heights = feasibility(hyperplanes=hypers,
                              c=1,
                              ambient_dim=ambient_dim,
                              backend=backend,
                              verbose=verbosity>=1)

    # calculate the frst
    if (heights is not None):
        reduced_heights = np.delete(heights, to_skip, 0)
        try:
            return Triangulation(p,
                                 p.labels_not_facet,
                                 heights=reduced_heights,
                                 make_star=make_star,
                                 check_heights=False)
        except:
            warnings.warn("This will be outdated with new CYTools pushes")
            return Triangulation(reduced_pts,
                                 poly=p,
                                 heights=reduced_heights,
                                 make_star=make_star,
                                 check_heights=False)
    else:
        if verbosity >= 1:
            print(f"Couldn't find heights for {cone.ambient_dim()}D cone with "
                    f"hyperplanes = {cone._hyperplanes}. Likely not-solid...")
        return None

def hypers_to_cone(hypers, dim=None, parse_inputs=False):
    """
    **Description:**
    Wrapper for Cone(hyperplanes=...) to be more resilient to failures.

    **Arguments:**
    - `hypers_in` *(array-like)*: The defining hyperplanes
    - `dim` *(integer, optional)*: The dimension of the hyperplane normals.
        Only used/needed if hypers is empty...
    - `parse_inputs` *(bool, optional)*: Whether to tell CYTools to 'parse' the
        input hyperplanes... only needed if using LIL or LIL_stack.

    **Returns:**
    *(Cone or None)* The corresponding cone (if possible). Otherwise, None.
    """
    # input checking
    if (isinstance(hypers, LIL_stack) and hypers.is_empty) or\
            (not isinstance(hypers, LIL_stack) and len(hypers)==0):
        if dim is None:
            raise ValueError("If no hyperplanes given, must set dimension")

        # no hyperplanes... define the trivial cone!
        rays = []

        for i in range(dim):
            unit = np.array([1 if j==i else 0 for j in range(dim)])
            rays.append(unit)
            rays.append(-unit)

        return Cone(rays=rays)

    # return the cone
    try:    return Cone(hyperplanes=hypers, parse_inputs=parse_inputs)
    except: raise Exception(f"Error in building cone with hypers={hypers}!")

# prefix with '_' to indicate that these shouldn't directly be called by user
def _cone_to_frt(cone=None,
                 p=None,
                 to_skip=None,
                 reduced_pts=None,
                 make_star=False,
                 backend=None,
                 verbosity=0):
    """
    **Description:**
    Convert a expanded-secondary cone to the corresponding FRST.

    **Arguments:**
    - `cone` *(Cone)*: The expanded-cpl cone.
    - `p` *(Polytope)*: The corresponding polytope.
    - `to_skip` *((list of) ints, optional)*: The indices to skip in making the
        FRST. These are non-origin points not in any 2-face.
    - `reduced_pts` *(list of lists, optional)*: The coords of vertices in p,
        after removing the vertices listed in to_skip.
    - `backend` *(string, optional)*: The backend to use for cone calculations.
    - `verbosity` *(int, optional)*: Verbosity level.

    **Returns:**
    *(Triangulaton)* The FRST.
    """
    # input checking
    if p is None:
        # if no poly is provided, we need
        #   1) reduced_pts
        #   2) either cone==None OR to_skip!=None
        assert reduced_pts is not None
        
        if to_skip is None:
            assert (cone is None)
        elif not isinstance(to_skip,list):
            to_skip = [to_skip]
    else:
        # p is provided... get to_skip and reduced_pts
        if verbosity>=1:
            print("_cone_to_frt: Calculating 'to_skip' and 'reduced_pts'...")

        if to_skip is None:
            to_skip = p.labels_facet
    
    # ensure to_skip is of the right format
    if (to_skip is not None) and (not isinstance(to_skip,list)):
        to_skip = [to_skip]

    if reduced_pts is None:
        pts = p.points()
        reduced_pts = np.delete(pts, to_skip, 0)

    # grab height vector
    if verbosity>=1:
        print("_cone_to_frt: Grabbing height vector...")

    if cone is None:
        if verbosity>=1:
            print("_cone_to_frt: No input cone -> assume all space!")
        if pts is None:
            pts = p.points()
        h = list(range(len(pts)))   # pick arbitrary height vector
    else:
        h = cone.find_interior_point(check=False,
                                     backend=backend,
                                     show_hints=False)

    # build the FRST
    if verbosity>=1:
        print("_cone_to_frt: Building the FRST...")

    if (h is not None):
        reduced_h = np.delete(h, to_skip, 0)
        try:
            return Triangulation(p, p.labels_not_facet, heights=reduced_h,\
                                        make_star=make_star, check_heights=False)
        except:
            warnings.warn("This will be outdated with new CYTools pushes")
            return Triangulation(reduced_pts, poly=p, heights=reduced_h,\
                                        make_star=make_star, check_heights=False)
    else:
        if verbosity >= 1:
            print(f"Couldn't find heights for {cone.ambient_dim()}D cone with "
                    f"hyperplanes = {cone._hyperplanes}. Likely not-solid...")
        return None

# extend face-triangulations to FRST
# ----------------------------------
def triangfaces_to_frt(triangs,
                       p,
                       make_star=False,
                       dicts=None,
                       npts=None,
                       backend=None,
                       verbosity=0):
    """
    **Description:**
    See/cite https://arxiv.org/abs/2309.10855

    Given a list of 2-face triangulations, construct an FRST that reduces to
    said triangulations.

    You can decide to not specify some of the 2-face triangulations. For this,
    just leave the associated element in triangs as None.

    **Arguments:**
    - `triangs` *(list of Triangulation)*: The 2-face triangulations. Elements
        can be None, in which case said 2-face is free.
    - `p` *(Polytope)* The full polytope of interst.
    - `dicts` *(list of dictionaries, optional)*: The maps from 2-face indices
        to entire polytope indices. Can be calculated automatically.
    - `npts` *(int, optional)*: The number of points in the entire polytope.
        Can be calculated automatically.
    - `backend` *(string, optional)*: The backend to use for cone calculations.
    - `verbosity` *(int, optional)*: Verbosity level.

    **Returns:**
    *(Triangulation)* The FRST obeying the specified 2-face triangulations.
    """
    # grab info if not provided
    if npts is None:
        npts = len(p.points())

    # remove the omitted triangs
    ts = [t for t in triangs if t is not None]

    # build the frst
    ineqs = permissible_heights(ts,npts=npts,no_duplicates=False,as_cone=False)
    frst = hypers_to_frt(ineqs, p, make_star=make_star,
                         backend=backend, verbosity=verbosity-1)
    return frst

def triangfaces_to_frst(triangs, p, dicts=None, npts=None, backend=None, verbosity=0):
    return triangfaces_to_frt(triangs=triangs,
                              p=p,
                              make_star=True,
                              dicts=dicts,
                              npts=npts,
                              backend=backend,
                              verbosity=verbosity)

# generate all 2-face inequivalent hyperplanes/cones/FRSTs
# --------------------------------------------------------
def triangface_ineqs(p, max_npts: int=17,
                     face_triangs: list=None,
                     N_face_triangs: int=1000,
                     triang_method = "grow2d",
                     return_triangs_dicts = False,
                     require_star: bool = False,
                     verbosity=0):
    """
    **Description:**
    Calculate all the 2-face triangulations.

    For each 2-face triangulation, generate the 'expanded-secondary' cone
    according to the restrictions imposed by 2-faces.

    **Arguments:**
    - `p` *(Polytope)*: The polytope of interest.

    **Returns:**
    *(list of lists of LIL)* List of faces. For each face, a list of all
        expanded-cpl cones.
    """
    npts = len(p.points())

    # find all 2-face triangulations
    if face_triangs is None:
        if verbosity>1:
            print(f"triangface_ineqs: Calculating face triangulations")
        face_triangs = p.compute_face_triangs(dim=2,only_regular=True,
                                    max_npts=max_npts,
                                    N_face_triangs=N_face_triangs,
                                    triang_method=triang_method,
                                    verbosity=verbosity-1)

    # iterate over faces
    if verbosity>1:
        print(f"triangface_ineqs: Calculating hyperplane inequalities")
    ineqs = []
    iter_wrapper = tqdm if verbosity>=1 else lambda x:x # (for progress bars)
    for f_triangs in iter_wrapper(face_triangs):
        ineqs.append([])

        # iterate over triangulations of this face
        for f_triang in f_triangs:
            tmp_ineqs = _2d_frt_cone_ineqs(f_triang, npts)
            if require_star:
                tmp_ineqs.append(_2d_s_cone_ineqs(f_triang, p, npts))
            ineqs[-1].append(tmp_ineqs)

    if return_triangs_dicts==False:
        return ineqs
    else:
        return ineqs, face_triangs
Polytope.triangface_ineqs = triangface_ineqs

def ntfe_hypers(p, N=None, seed=None,
                max_npts: int=17,
                face_ineqs: list=None,
                face_triangs: list=None,
                N_face_triangs: int=1000,
                triang_method = "grow2d",
                require_star = False,
                as_generator = False,
                separate_boring = True,
                verbosity=0):
    """
    **Description:**
    See/cite https://arxiv.org/abs/2309.10855

    Generate the hyperplane normals corresponding to all expanded-secondary
    cones for each set of 2-face triangulations of our polytope.

    **Arguments:**
    - `p` *(Polytope)*: The polytope of interest
    - `n` *(int, optional)*: Number of 2-face triangulation sets to consider,
        if not all are wanted
    - `seed` *(int, optional)*: The random seed used for randomly choosing
        hypers (if n<)
    - `verbosity` *(int, optional)*: Verbosity level.

    **Returns:**
    *(list of LIL)* The list of hyperplanes defining each
        'expanded-secondary cone'
    """
    if seed is None:
        seed = time.time_ns()%(2**32)

    # grab the cpl-cone inequalities
    if face_ineqs is None:
        if verbosity>=1:
            print("ntfe_hypers: Constructing hyperplanes for faces...")
        ineqs_array = p.triangface_ineqs(max_npts=max_npts,
                                         face_triangs=face_triangs,
                                         N_face_triangs=N_face_triangs,
                                         triang_method=triang_method,
                                         require_star = require_star,
                                         verbosity=verbosity-1)
    else:
        ineqs_array = face_ineqs.copy()

    if separate_boring:
        ineqs_boring = []

        i=0
        while i<len(ineqs_array):
            if len(ineqs_array[i])==1:
                ineqs_boring.append(ineqs_array.pop(i)[0])
            else:
                i+=1

        if len(ineqs_boring):
            ineqs_boring = sum(ineqs_boring[1:],ineqs_boring[0])
            ineqs_array.append([ineqs_boring])

    # get number of triangulations per 2-face
    if verbosity>=1:
        print("ntfe_hypers: Calculating total number of ineqs...")
    choices_counts = list(map(len,ineqs_array))
    choices = list(map(range, choices_counts))

    # for each set of 2-face triangulations, group the inequalities
    if verbosity>=1:
        print("ntfe_hypers: Intersecting face H-cones...", end=' ')
        print(f"(there are {math.prod(choices_counts)} total)")
    if (N is None) or (N >= math.prod(choices_counts)):
        if verbosity>=2:
            print(f"ntfe_hypers: Calculating all N={math.prod(choices_counts)} "\
                   "intersections...")
        chosen = range(math.prod(choices_counts))
    else:
        if True:
            # choose uniformly on chromosones
            chosen = set()

            # set the seed
            np.random.seed(seed)

            # choose the hypers
            while len(chosen) < N:
                choice = tuple(np.random.choice(x) for x in choices)
                chosen.add(choice)
        else:
            chosen = np.random.choice(math.prod(choices_counts), N, replace=False)
    
    # grab/return hyperplaners
    if as_generator:
        def gen():
            for choice in chosen:
                yield LIL_stack(ineqs_array, choice, choices_counts)
        return gen()
    else:
        if verbosity>=1:
            hypers = [LIL_stack(ineqs_array,choice,choices_counts) for choice in tqdm(chosen)]
        else:
            hypers = [LIL_stack(ineqs_array,choice,choices_counts) for choice in chosen]

        return hypers
Polytope.ntfe_hypers = ntfe_hypers

def ntfe_cones(p, N:int=None, hypers=None, seed: int=None,
               max_npts: int=17,
               face_ineqs: list=None,
               face_triangs: list=None,
               N_face_triangs: int=1000,
               triang_method = "grow2d",
               require_star = False,
               as_generator = False,
               verbosity=0):
    """
    **Description:**
    See/cite https://arxiv.org/abs/2309.10855

    Generate the expanded-secondary cones corresponding to the input
    hyperplanes.

    **Arguments:**
    - `hypers` *(list of LIL, optional)*: The hyperplanes defining the
        cones. If no hyperplanes are input, these are automatically calculated.
    - `p` *(Polytope, optional)*: The polytope of interest. Only used if
        hypers is None.
    - `verbosity` *(boolean, optional)*: Verbosity level.

    **Returns:**
    *(list of Cones)* The list of the expanded-cpl cones.
    """
    if seed is None:
        seed = time.time_ns()%(2**32)

    # input checking
    if hypers is None:
        if p is None:
            raise ValueError("No hyperplanes or polytope was input... "+\
                                            "One of the two must be input...")
        else:
            hypers = p.ntfe_hypers(N, max_npts=max_npts,
                                   face_ineqs=face_ineqs,
                                   face_triangs=face_triangs,
                                   N_face_triangs=N_face_triangs,
                                   seed=seed,
                                   triang_method=triang_method,
                                   require_star = require_star,
                                   as_generator = as_generator,
                                   verbosity=verbosity-1)
            dim = len(p.points())
    else:
        dim = None
        if isinstance(hypers[0],LIL_stack):
            if not hypers[0].is_empty:
                dim = len(hypers[0][0])
        elif len(hypers[0]):
            dim = len(hypers[0][0])
        
        if dim is None:
            dim = len(p.points())

    # if returning a generator, just do so here
    if as_generator:
        def gen():
            for hyper in hypers:
                yield hypers_to_cone(hyper, dim=dim, parse_inputs=False)
        return gen()

    # for each hyperplane, calculate the corresponding cone
    n_hypers = len(hypers)

    # randomly sample hypers:
    if (N is not None) and (N<n_hypers):
        hyper_inds = list(range(n_hypers))

        # shuffle the indices and select the first N
        random.seed(seed)
        random.shuffle(hyper_inds)
        hyper_inds = hyper_inds[:N]

        iterator = [hypers[i] for i in hyper_inds]
    else:
        iterator = hypers

    # convert hyperplanes to cones
    cones = []
    iter_wrapper = tqdm if verbosity>=1 else lambda x:x # (for progress bars)
    for hyper in iter_wrapper(iterator):
        cones.append(hypers_to_cone(hyper, dim=dim, parse_inputs=False))

    return cones
Polytope.ntfe_cones = ntfe_cones

# (hack to get multiprocessing to work)
_func = None
def worker_init(func):
    global _func
    _func = func
def worker(x):
    return _func(x)

def ntfe_frts(p: "Polytope",
              N: int=None,
              make_star: bool=False,
              cones: list=None,
              hypers: list=None,
              face_ineqs: list=None,
              face_triangs: list=None,
              max_npts: int=17,
              seed: int=None,
              N_face_triangs: int=1000,
              triang_method: str="fast",
              as_generator: bool=False,
              nproc: int=mp.cpu_count(),
              backend: str=None,
              verbosity: int=0):
    """
    **Description:**
    See/cite https://arxiv.org/abs/2309.10855

    Generate NTFE FRTs of a polytope.

    **Arguments:**
    - `p`: The polytope of interest.
    - `N`: The number of NTFEs to generate. If None, then all are generated.
    - `make_star`: Whether to make the triangulations star.
    - `cones`: The list of the expanded-cpl cones. If no cones are input,
        these are automatically calculated.
    - `hypers`: The hyperplanes defining the cones. Only needed if cones is
        None. If no hyperplanes are input, these are automatically calculated.
    - `face_ineqs`: The inequalities associated to each 2-face triangulation.
        Basically the same as `hypers`...
    - `face_triangs`: The triangulations of each 2-face.
    - `max_npts`: The maximum number of points in a 2-face s.t. we try to
        calculate *all* FRTs of it. Otherwise, just sample.
    - `seed`: The random seed to use for generating 2-face triangulations
    - `triang_method`: The 2-face triangulation method
    - `as_generator`: Whether to return the FRTs via a generator instead of a
        list. Useful for memory concerns.
    - `nproc`: The number of processors to use when constructing FRTs.
    - `backend`: The backend to use for cone calculations.
    - `verbosity`: Verbosity level.

    **Returns:**
    The FRTs.
    """
    warnings.warn("Shift to not having to construct formal Cone object")
    if seed is None:
        seed = time.time_ns()%(2**32)

    # grab cones, if not provided
    if verbosity>=1:
        print("ntfe_frts: Calculating expanded-cpl cones...")

    if hypers is not None:
        data = hypers
        data_to_frst = hypers_to_frt
    elif cones is not None:
        data = cones
        data_to_frst = _cone_to_frt
    else:
        data = p.ntfe_hypers(N,
                             max_npts=max_npts,
                             seed=seed,
                             face_ineqs=face_ineqs,
                            face_triangs=face_triangs,
                             N_face_triangs=N_face_triangs,
                             triang_method=triang_method,
                             as_generator=as_generator,
                             verbosity=verbosity-1)
        data_to_frst = hypers_to_frt


    # randomly select N cones/hyperplanes
    # (might get fewer than N FRSTs, in case the cones aren't all solid)
    if N is not None:
        random.seed(seed)
        random.shuffle(data)
        data = data[:N]

    # find indices not in any 2-face
    if verbosity>=1:
        print("ntfe_frts: Calculating to_skip, reduced_pts"+\
                "(same for all FRSTs, so it's good to do it here...)")

    to_skip = p.labels_facet
    reduced_pts = np.delete(p.points(), to_skip, 0)

    # if returning a generator, just do so here
    if as_generator:
        def gen():
            for datum in data:
                frst = data_to_frst(datum, p, to_skip=to_skip,\
                            reduced_pts=reduced_pts,\
                            make_star=make_star,\
                            backend=backend, verbosity=verbosity-1)
                if (frst is None) or (frst==False):
                    continue
                else:
                    yield frst
        return gen()

    # for each expanded-cpl cone, calculate the corresponding FRST
    frsts = []
    def func(datum):
        return data_to_frst(datum, p, to_skip=to_skip,\
                        reduced_pts=reduced_pts,\
                        make_star=make_star,\
                        backend=backend, verbosity=verbosity-1)

    with mp.Pool(nproc, initializer=worker_init, initargs=(func,)) as p:
        for frst in p.imap(worker, data, chunksize=max(1, (len(data)+1)//nproc)):
            if (frst is None) or (frst==False):
                continue

            frsts.append(frst)

    return frsts
Polytope.ntfe_frts = ntfe_frts

def ntfe_frsts(p: "Polytope",
               N: int=None,
               cones: list=None,
               hypers: list=None,
               face_ineqs: list=None,
               face_triangs: list=None,
               max_npts: int=17,
               seed: int=None,
               N_face_triangs: int=1000,
               triang_method: str="fast",
               as_generator: bool=False,
               nproc: int=mp.cpu_count(),
               backend: str=None,
               verbosity: int=0):
    """
    **Description:**
    See/cite https://arxiv.org/abs/2309.10855

    Generate NTFE FRTs of a polytope.

    **Arguments:**
    - `p`: The polytope of interest.
    - `N`: The number of NTFEs to generate. If None, then all are generated.
    - `cones`: The list of the expanded-cpl cones. If no cones are input,
        these are automatically calculated.
    - `hypers`: The hyperplanes defining the cones. Only needed if cones is
        None. If no hyperplanes are input, these are automatically calculated.
    - `face_ineqs`: The inequalities associated to each 2-face triangulation.
        Basically the same as `hypers`...
    - `face_triangs`: The triangulations of each 2-face.
    - `max_npts`: The maximum number of points in a 2-face s.t. we try to
        calculate *all* FRTs of it. Otherwise, just sample.
    - `seed`: The random seed to use for generating 2-face triangulations
    - `triang_method`: The 2-face triangulation method
    - `as_generator`: Whether to return the FRTs via a generator instead of a
        list. Useful for memory concerns.
    - `nproc`: The number of processors to use when constructing FRTs.
    - `backend`: The backend to use for cone calculations.
    - `verbosity`: Verbosity level.

    **Returns:**
    The FRTs.
    """
    return ntfe_frts(p,
                     N=N,
                     make_star=True,
                     cones=cones,
                     hypers=hypers,
                     face_ineqs=face_ineqs,
                     face_triangs=face_triangs,
                     max_npts=max_npts,
                     seed=seed,
                     N_face_triangs=N_face_triangs,
                     triang_method=triang_method,
                     as_generator=as_generator,
                     nproc=nproc,
                     backend=backend,
                     verbosity=verbosity)
Polytope.ntfe_frsts = ntfe_frsts
