# -----------------------------------------------------------------------------
# Description:  Wrappers for cygv. No new functionality, just ease of use.
# -----------------------------------------------------------------------------

from cytools.calabiyau import CalabiYau

try:
    import cygv
except:
    import subprocess
    import sys
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'cygv'])
    import cygv

# generic function
def _compute_gvs_gws(self,
                     gv_or_gw: str,
                     grading_vec: "ArrayLike"=None,
                     max_deg: bool=None,
                     min_points: bool=None):
    """
    **Description:**
    Wrapper for cygv GV and GW computations. A method of cytools.CalabiYau
    
    NOT INTENDED TO BE CALLED DIRECTLY!

    **Arguments:**
    - `gv_or_gw`: String specifying whether 'gv' or 'gw' computations are
        performed.
    - `grading_vec`: The grading vector to use in the computations. A default
        is chosen if none is provided.
    - `max_deg`: The maximum degree to compute GVs/GWs to. Must be specified
        iff min_points=None.
    - `min_points`: The minimum number of GVs/GWs to compute. Must be
        specified iff max_deg=None.

    **Returns:**
    The GV/GW invariants.
    """
    # check that the user either passed max_deg or min_points
    if not (max_deg is None) ^ (min_points is None):
        raise ValueError("Either max_deg or min_points must be set!")
    
    # get basics
    kappa = self.intersection_numbers(in_basis=True, format='coo')
    glsm = self.curve_basis(include_origin=False, as_matrix=True)
    mori = self.toric_mori_cone(in_basis=True)
    generators = mori.rays()
    
    # compute a grading vector if none is provided
    if grading_vec is None:
        grading_vec = mori.find_grading_vector()
        
    # compute the GVs
    if gv_or_gw=='gv':
        fct = cygv.compute_gv
    else:
        fct = cygv.compute_gw

    return dict(fct(generators = mori.rays(),
                    grading_vector = grading_vec,
                    q = self.curve_basis(include_origin=False, as_matrix=True),
                    intnums = self.intersection_numbers(in_basis=True, format='coo'),
                    max_deg = max_deg,
                    min_points = min_points))
CalabiYau._compute_gvs_gws = _compute_gvs_gws

# GVs
def compute_gvs(self,
                grading_vec: "ArrayLike"=None,
                max_deg: bool=None,
                min_points: bool=None):
    """
    **Description:**
    Wrapper for cygv GV computations. A method of cytools.CalabiYau

    **Arguments:**
    - `grading_vec`: The grading vector to use in the computations. A default
        is chosen if none is provided.
    - `max_deg`: The maximum degree to compute GVs to. Must be specified iff
        min_points=None.
    - `min_points`: The minimum number of GVs/GWs to compute. Must be
        specified iff max_deg=None.

    **Returns:**
    The GV invariants.
    """
    return self._compute_gvs_gws('gv', grading_vec, max_deg, min_points)
CalabiYau.compute_gvs = compute_gvs

# GWs
def compute_gws(self, grading_vec=None, max_deg=None, min_points=None):
    """
    **Description:**
    Wrapper for cygv GW computations. A method of cytools.CalabiYau

    **Arguments:**
    - `grading_vec`: The grading vector to use in the computations. A default
        is chosen if none is provided.
    - `max_deg`: The maximum degree to compute GWs to. Must be specified iff
        min_points=None.
    - `min_points`: The minimum number of GWs to compute. Must be specified
        iff max_deg=None.

    **Returns:**
    The GW invariants.
    """
    return self._compute_gvs_gws('gw', grading_vec, max_deg, min_points)
CalabiYau.compute_gws = compute_gws
