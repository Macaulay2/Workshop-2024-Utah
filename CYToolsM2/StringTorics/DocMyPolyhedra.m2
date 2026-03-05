-- Documentation for polytope utility functions (MyPolyhedra.m2)

doc ///
  Key
    vertexList
    (vertexList, Polyhedron)
  Headline
    sorted list of vertices of a polytope
  Usage
    L = vertexList P
  Inputs
    P:Polyhedron
  Outputs
    L:List
      of lists of integers, the vertices of $P$ sorted lexicographically
  Description
    Text
      Returns the vertices of a polytope as a sorted list of lists of integers.
      This provides a canonical ordering of the vertices, unlike the @TO (vertices, Polyhedron)@
      function from the Polyhedra package which returns a matrix whose column order may vary.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      vertexList P
      vertexMatrix P
  SeeAlso
    vertexMatrix
    (vertices, Polyhedron)
///

doc ///
  Key
    vertexMatrix
    (vertexMatrix, Polyhedron)
  Headline
    matrix whose columns are the sorted vertices of a polytope
  Usage
    M = vertexMatrix P
  Inputs
    P:Polyhedron
  Outputs
    M:Matrix
      whose columns are the vertices of $P$, in the order given by @TO vertexList@
  Description
    Text
      Returns the matrix whose columns are the vertices of $P$, sorted in the
      canonical order given by @TO vertexList@.  This is the transpose of
      @TT "matrix vertexList P"@.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      vertexList P
      vertexMatrix P
  SeeAlso
    vertexList
///

doc ///
  Key
    faceDimensionHash
    (faceDimensionHash, Polyhedron)
  Headline
    hash table of faces and their dimensions
  Usage
    H = faceDimensionHash P
  Inputs
    P:Polyhedron
  Outputs
    H:HashTable
      whose keys are faces (as sorted lists of vertex indices) and values are their dimensions
  Description
    Text
      Returns a hash table mapping each face of $P$ (given as a sorted list of
      indices into @TO vertexList@) to its dimension.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      H = faceDimensionHash P
  SeeAlso
    faceList
    vertexList
///

doc ///
  Key
    faceList
    (faceList, Polyhedron)
    (faceList, ZZ, Polyhedron)
  Headline
    list of faces of a polytope
  Usage
    L = faceList P
    L = faceList(d, P)
  Inputs
    d:ZZ
      the dimension of faces to return (optional)
    P:Polyhedron
  Outputs
    L:List
      of faces, each given as a sorted list of vertex indices
  Description
    Text
      With one argument, returns all faces of $P$ sorted by dimension.
      With two arguments, returns only the faces of dimension $d$.
      Each face is represented as a sorted list of indices into @TO vertexList@.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      faceList P
      faceList(0, P)
      faceList(1, P)
  SeeAlso
    faceDimensionHash
    annotatedFaces
///

doc ///
  Key
    dualFace
    (dualFace, Polyhedron, List)
  Headline
    the dual face in the polar polytope
  Usage
    g = dualFace(P, f)
  Inputs
    P:Polyhedron
      a reflexive polytope
    f:List
      a sorted list of vertex indices giving a face of $P$
  Outputs
    g:List
      a sorted list of vertex indices giving the dual face of @TT "polar P"@
  Description
    Text
      Given a face $f$ of $P$ (as a sorted list of vertex indices),
      returns the dual face of the polar polytope.
      For a reflexive polytope, if $f$ has dimension $d$, the dual face has
      dimension $n - 1 - d$, where $n = $ @TT "dim P"@.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      vertexList P
      vertexList polar P
      f = {0,1,2}
      dualFace(P, f)
  SeeAlso
    minimalFace
    annotatedFaces
///

doc ///
  Key
    minimalFace
    (minimalFace, Polyhedron, List)
  Headline
    the smallest face of a polytope containing given points
  Usage
    f = minimalFace(P, pts)
  Inputs
    P:Polyhedron
    pts:List
      either a single lattice point (list of integers) or a list of lattice points
  Outputs
    f:List
      a sorted list of vertex indices giving the minimal face containing @TT "pts"@
  Description
    Text
      Returns the minimal face of $P$ that contains the given lattice point(s).
      The face is returned as a sorted list of indices into @TO vertexList@.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      vertexList P
      minimalFace(P, {-1,-1,-1}) 
      minimalFace(P, {0,0,-1}) -- WRONG?!
      minimalFace(P, {0,0,0})
  SeeAlso
    dualFace
    latticePointList
///

doc ///
  Key
    latticePointList
    (latticePointList, Polyhedron)
    (latticePointList, Polyhedron, List)
  Headline
    lattice points of a polytope, ordered by face dimension
  Usage
    L = latticePointList P
    L = latticePointList(P, f)
  Inputs
    P:Polyhedron
    f:List
      a sorted list of vertex indices giving a face of $P$ (optional)
  Outputs
    L:List
      of lattice points (if no face given) or of indices into the lattice point list
  Description
    Text
      With one argument, returns all lattice points of $P$ as a list of lists of integers,
      sorted first by the dimension of the smallest face containing the point.
      Vertices come first, then points on edges, then points on 2-faces, etc.

      With a face $f$ given, returns the list of indices (into @TT "latticePointList P"@)
      of lattice points lying on that face.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      latticePointList P
      f = {0,1,3,4}  -- a 2-face,
      latticePointList(P, f)
      f = {0,1,4} -- contained in a 2-face
      latticePointList(P, f)
  SeeAlso
    latticePointHash
    interiorLatticePointList
    (latticePoints, Polyhedron)
///

doc ///
  Key
    latticePointHash
    (latticePointHash, Polyhedron)
  Headline
    hash table from lattice points to their indices
  Usage
    H = latticePointHash P
  Inputs
    P:Polyhedron
  Outputs
    H:HashTable
      mapping each lattice point (as a list of integers) to its index in @TO latticePointList@
  Description
    Text
      Returns a hash table that maps each lattice point of $P$ to its index in
      the list returned by @TO latticePointList@.  This is useful for quickly
      looking up the index of a given lattice point.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      latticePointList P
      latticePointHash P
  SeeAlso
    latticePointList
///

doc ///
  Key
    interiorLatticePointList
    (interiorLatticePointList, Polyhedron, List)
  Headline
    lattice points in the relative interior of a face
  Usage
    L = interiorLatticePointList(P, f)
  Inputs
    P:Polyhedron
    f:List
      a sorted list of vertex indices giving a face of $P$
  Outputs
    L:List
      of indices into @TO latticePointList@ of points in the relative interior of $f$
  Description
    Text
      Returns the list of indices (into @TT "latticePointList P"@)
      of lattice points lying in the relative interior of the given face $f$.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      faceList(3, P)
      interiorLatticePointList(P, {0,1,3,4})
  SeeAlso
    latticePointList
    annotatedFaces
///

doc ///
  Key
    latticePointsAndDimensions
    (latticePointsAndDimensions, Polyhedron)
  Headline
    lattice points paired with their minimal face dimensions
  Usage
    (LP, dims) = latticePointsAndDimensions P
  Inputs
    P:Polyhedron
  Outputs
    LP:List
      the lattice points of $P$ (same as @TO latticePointList@)
    dims:List
      the dimension of the minimal face containing each lattice point
  Description
    Text
      Returns a pair: the lattice points of $P$ (as from @TO latticePointList@)
      and a parallel list of integers giving the dimension of the smallest face
      containing each lattice point.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      (LP, dims) = latticePointsAndDimensions P
      netList(prepend({"point", "face dimension"}, transpose {LP, dims}))
  SeeAlso
    latticePointList
    minimalFace
///

doc ///
  Key
    (annotatedFaces, Polyhedron)
  Headline
    annotated face data for a polytope
  Usage
    L = annotatedFaces P
  Inputs
    P:Polyhedron
      a reflexive polytope
  Outputs
    L:List
      each entry is a list $\{d, f, lps, nint, ndualint\}$
  Description
    Text
      Returns a sorted list of all faces of $P$, each annotated with:
      the dimension $d$ of the face, the vertex indices $f$, the lattice point
      indices $lps$, the number of interior lattice points $nint$, and the
      number of interior lattice points of the dual face in @TT "polar P"@.

      This is the @TO Polyhedron@ version.  See also @TO (annotatedFaces, ReflexivePolytope)@.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      netList annotatedFaces P
  SeeAlso
    (annotatedFaces, ReflexivePolytope)
    (annotatedFaces, ZZ, Polyhedron)
    faceList
    interiorLatticePointList
///

doc ///
  Key
    (annotatedFaces, ZZ, Polyhedron)
  Headline
    annotated face data for faces of a given dimension
  Usage
    L = annotatedFaces(d, P)
  Inputs
    d:ZZ
      the dimension of faces to annotate
    P:Polyhedron
      a reflexive polytope
  Outputs
    L:List
      each entry is a list $\{f, lps, nint, ndualint\}$
  Description
    Text
      Returns a sorted list of faces of dimension $d$ in $P$, each annotated with:
      the vertex indices $f$, the lattice point indices $lps$, the number of
      interior lattice points $nint$, and the number of interior lattice points
      of the dual face in @TT "polar P"@.

      Note: the dimension entry is omitted since it is given as input.
    Example
      P = convexHull transpose matrix {{-1, -1, -1}, {-1, 3, -1}, {0, -1, 2}, {1, 1, -1}, {2, -1, -1}}
      netList annotatedFaces(1, P) -- edges
      netList annotatedFaces(2, P) -- facets
  SeeAlso
    (annotatedFaces, Polyhedron)
    faceList
    interiorLatticePointList
///

doc ///
  Key
    (automorphisms, Polyhedron)
  Headline
    lattice automorphisms of a polytope
  Usage
    L = automorphisms P
  Inputs
    P:Polyhedron
  Outputs
    L:List
      of matrices, each an integer matrix giving a lattice automorphism of $P$
  Description
    Text
      Returns the list of all integer matrices $M$ such that $M$ maps the
      vertex set of $P$ to itself (i.e., lattice isomorphisms from $P$ to $P$).
    Example
      P = convexHull matrix {{1, -1, -1, -1, -1, -1, -1},
          {-1, 3, 3, 3, -1, -1, -1},
          {0, -1, 0, -1, 0, 1, 1},
          {0, -2, -1, -1, 1, 1, 2}}
      (dim P, isReflexive P) -- a reflexive polytope of dimension 4
      automorphisms P
      # automorphisms P -- 48 automorphisms.
  SeeAlso
    (automorphisms, ReflexivePolytope)
    automorphismsAsPermutations
///
