-- TODO? change ReflexivePolytope, CYPolytope to CanonicalPolytope?

doc ///
  Key
    ReflexivePolytope
  Headline
    data for a reflexive polytope
  Description
    Text
      This type encapsulates the data in a reflexive polytope.  In a toric variety setting,
      you should think of this polytope as being in the "N" lattice.  An object of this
      type can contain cached information, which is useful when constructing Calabi-Yau varieties
      from this polytope (usually as hypersurfaces or complete intersections in associated
      toric varieties, see @TO "Batyrev construction"@.
    Text
      A {\tt ReflexivePolytope} can be constructed from vertices of a reflexive polytope.
    Text
      For example, let's start with a dimension 3 reflexive polytope: the cube.
      We use @TO cyPolytope@ to create the corresponding Macaulay2 object.
    Example
      verts = {
          {-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
          {-1, -1, 1}, {1, -1, 1}, {-1, 1, 1}, {1, 1, 1}}
      Q = reflexivePolytope verts
      dim Q
      --isWellDefined Q
    Text
      Given a {\tt ReflexivePolytope}, the method @TO (rays,ReflexivePolytope)@
      returns a list of the boundary lattice points
      not in facets.  This might seem like a peculiar definition, but we often consider a triangulation of
      this point set (e.g. of the reflexive polytope), and then consider fans whose cones are the cones over
      the faces in this triangulation.  The rays in this fan are precisely the output of {\tt rays}.
      Each lattice point is contained in some minimal face of the polytope.
      {\tt faceDimensions} gives the list of the dimensions of this minimal face, over all "rays" of $Q$.
    Example
      rays Q
      faceDimensions Q
      tally oo
    Text
      Note that there are 8 vertices, and 12 lattice points in edges.
    Text
      We can obtain the polyhedron (as an object from the @TO Polyhedra@ package).
      TODO: downplay this?
    Example
      P2 = polytope Q
      isReflexive P2
      vertices P2
      matrix {latticePoints P2}
      entries transpose oo
      latticePoints Q
    Text
      Notice that latticePoints of the {\tt Polyhedra} object, {\tt P2} are in a different order
      than the lattice points from the {\tt ReflexivePolytope} $Q$  All indices for a ReflexivePolytope
      object (faces, maximal cones, etc) {\it all} refer to the order of rays/lattice points from $Q$.
  SeeAlso
      CalabiYauInToric
      reflexivePolytope
///
 
doc ///
  Key
    reflexivePolytope
    (reflexivePolytope, Matrix)
    (reflexivePolytope, List)
    (reflexivePolytope, Polyhedron)
    (reflexivePolytope, KSEntry)
  Headline
    construction of a ReflexivePolytope (reflexive polytope)
  Usage
    Q = reflexivePolytope A
  Inputs
    A:{Matrix,List,Polyhedron,KSEntry}
  Outputs
    Q:ReflexivePolytope
  Description
    Text
      A ReflexivePolytope contains the data associated to a reflexive
      polytope.  One is most easily constructed from vertices for this
      polytope (given either an integer matrix whose columns are these
      vertices, or a list of these vertices, or a @TO "Polyhedron"@
      object from the @TO "Polyhedra"@ package.

      This function does {\it not} check that the result is a reflexive polytope.
      You must do this yourself if you are want to check.
      Use "TO (isWellDefined, ReflexivePolytope)" or "TO (isReflexive, ReflexivePolytope)" to check this.
      
      One direct way to construct this reflexive polytope is to give it a
      list of (integer) vertices or various points.  The function takes the
      convex hull of these points, and returns the reflexive polytope.
    Example
      verts = {{-1, -1, -1, 1}, {-1, -1, -1, 0}, {-1, -1, 0, 0}, {-1, 0, -1, 0},
          {0, -1, -1, 0}, {0, 0, -1, 1}, {1, 1, 2, -1}}
      Q = reflexivePolytope verts
      vertices Q
      P2 = polytope Q
      isReflexive P2
    Text
      Instead of a list of points, you may also give a matrix whose columns are
      integral points.  The convex hull of these is constructed.
    Example
      pts = transpose matrix latticePoints Q
      Q1 = reflexivePolytope pts
      vertices Q1 === vertices Q
    Text
      If you call reflexivePolytope with a polytope as its argument (which happens to be
      reflexive), you get a ReflexivePolytope referring to the same polytope.
    Example
      P2 = polytope Q
      Q2 = reflexivePolytope P2
      vertices Q2 == vertices Q
    Text
      Finally, perhaps the easiest way to generate reflexive polytopes (of dimension 4)
      is to use the Kreuzer-Skarke database.  This one assumes that the entry is for the
      polytope on the "M" lattice side, so it computes the polar dual on the "N" lattice side.

      This particular example happens to match the input above.
      For dimension 4 polytopes, the Kreuzeer-Skarke database is generally
      arranged by a number $h^{1,1}$, which is information
      corresponding to a Calabi-Yau hypersurface in a toric variety
      constructed via this polytope.
    Example
      topes = kreuzerSkarke 3;
      tope = topes_100
      Q3 = reflexivePolytope tope
      vertices Q3 == vertices Q
      --(hh^(1,1) Q3, hh^(1,2) Q3)
  SeeAlso
    (vertices, ReflexivePolytope)
    (latticePoints, ReflexivePolytope)
    (rays, ReflexivePolytope)
    (annotatedFaces, ReflexivePolytope)
///

doc ///
  Key
    (latticePoints, ReflexivePolytope)
  Headline
    lattice points on a reflexive polytope
  Usage
    latticePoints Q
  Inputs
    Q:ReflexivePolytope
  Outputs
    :List
      of lists of integers: the lattice points, sorted by (increasing) dimension of the minimal face they sit on
  Description
    Text
      This function returns the list of lattice points of a @TO ReflexivePolytope@.
      The lattice points are sorted in some manner, such that the minimalface dimension of a lattice
      point is monotone increasing in this list.  As the origin is the only interior lattice
      point, its face dimension is the dimension of the polytope.

      Here is an example, taken from the Kreuzer-Skarke database.
    Example
      verts = {
          {0, 0, -1, -1},
          {-1, -1, -1, 0},
          {-1, -1, -1, 1},
          {-1, -1, 0, 0},
          {-1, 0, -1, -1},
          {0, -1, -1, 0},
          {0, 0, -1, 2},
          {0, 0, 2, -1},
          {0, 1, -1, -1},
          {1, 1, -1, 2}
          }
      Q = reflexivePolytope verts
      transpose matrix latticePoints Q
      faceDimensions Q -- this should match latticePoints...
      for i from 0 to #latticePoints Q - 1 list faceDimension(Q, i)
      #rays Q === #latticePoints Q - 3
      sort vertices Q == sort verts
      #rays Q
  SeeAlso
///

///
  Key
    (polytope, ReflexivePolytope, String)
    (polytope, ReflexivePolytope)
  Headline
  Usage
  Inputs
  Outputs
  Consequences
  Description
    Text
    Example
  Caveat
  SeeAlso
///

doc ///
  Key
    (annotatedFaces, ReflexivePolytope)
  Headline
    a list of faces of a reflexive polytope together with lattice point information
  Usage
    annotatedFaces Q
  Inputs
    Q:ReflexivePolytope
  Outputs
    :List
      each entry is a list containing: the dimension of the face, the indices of the
      vertices, the indices of all (boundary) lattice points in the face, the number
      of interior points in the face, and the number of interior points in the dual face
  Description
    Text
    Example
      verts = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, -1}, {0, -1, 0}, {-1, 0, 0}, {-1, 1, 0}}
      Q = reflexivePolytope verts
      transpose matrix latticePoints Q
      faceDimensions Q
      netList annotatedFaces Q
      polar Q
      faceDimensions polar Q
      netList annotatedFaces polar Q
    Text
      This is how we obtained the vertices to create Q.
    Example
      rays smoothFanoToricVariety(3, 12)
  SeeAlso
    (latticePoints, ReflexivePolytope)
    (faceDimensions, ReflexivePolytope)
    (polar, ReflexivePolytope)
///

doc ///
  Key
    (dim, ReflexivePolytope)
  Headline
    dimension of a reflexive polytope
  Usage
    dim Q
  Inputs
    Q:ReflexivePolytope
  Outputs
    :ZZ
      the dimension of this polytope
  Description
    Text
    Example
      verts = rays smoothFanoToricVariety(4, 30)
      --verts1 = for v in verts list append(v, sum v) -- 
      Q = reflexivePolytope verts
      rays Q
      isReflexive polytope Q
      hh^(1,1) Q
      hh^(1,2) Q
  SeeAlso
///

///
  --Some tests of ReflexivePolytopes stuff
restart
needsPackage "StringTorics"
  -- need to be in StringTorics package currently for this:
  topes = value get "topes-h11-5";
  #topes
  Q = reflexivePolytope topes#2000
  vertices Q
  rays Q
  netList annotatedFaces Q
  findAllFRSTs Q
  findAllFRVTs Q
  partitionFRSTsByDFaceEquivalence(2, Q)
///
