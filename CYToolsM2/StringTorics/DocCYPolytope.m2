doc ///
    Key 
        "working with reflexive polytopes"
    Headline
        information about the basic constructors
    Description
        Text
          A reflexive polytope is a lattice polytope (all vertices are
          integer points), whose polar dual is also a lattice
          polytope.  A reflexive polytope always has exactly one
          lattice point in the interior, which we take to be the
          origin.  The polar dual of a reflexive polytope is also a
          reflexive polytope.
        Text
          An important way to make reflexive polytopes is to
          
    	Text
    	    @SUBSECTION "Constructors for reflexive polytopes"@
	Text
    	    @UL {
        	TO reflexivePolytope
    	    }@
    	Text
    	    @SUBSECTION "Useful functions"@
	Text
    	    @UL {
--                TO (normalForm, ReflexivePolytope),
        	TO (polar, ReflexivePolytope)
    	    }@
    	Text
    	    @SUBSECTION "Face structure and combinatorics"@
	Text
    	    @UL {
                TO (vertices, ReflexivePolytope),
                TO (dim, ReflexivePolytope),
        	TO (annotatedFaces, ReflexivePolytope),
                TO (latticePoints, ReflexivePolytope),
                TO (faceDimensions, ReflexivePolytope),
                TO (faceDimension, ReflexivePolytope, ZZ),
                TO (automorphisms, ReflexivePolytope),
                TO (automorphismsAsPermutations, ReflexivePolytope)
--                TO (genus, ReflexivePolytope, List)
--                should have Erhart polynomial too!
            }@
    	Text
    	    @SUBSECTION "Triangulations of a reflexive polytope"@
        Text
            A (Batyrev) Calabi-Yau n-fold is determined by a (n+1)-dimensional reflexive
            polytope, and a fine, regular, star triangulation of the set of all
            lattice points not in the interior of facets.
	Text
    	    @UL {
                TO (findOneFRST, ReflexivePolytope),
        	TO (findAllFRSTs, ReflexivePolytope),
                TO isTriangulationOfPolytope,
    	    }@
        Text
    	    @SUBSECTION "Features common to all Calabi-Yau hypersurfaces arising from a polytope"@
	Text
    	    @UL {
        	TO (degrees, ReflexivePolytope),
                TO (rays, ReflexivePolytope),
                TO (basisIndices, ReflexivePolytope),
                TO (isFavorable, ReflexivePolytope),
                TO (hh, Sequence, ReflexivePolytope),
    	    }@
        
    SeeAlso
        "creating and using Batyrev Calabi-Yau hypersurfaces"
        "working with intersection numbers and intersection rings"
        "working with complete intersections in toric varieties"
        "working with cohomology of line bundles and sheaves in torics"
        "topological equivalence of Calabi-Yau 3-folds"
        "Gopakumar-Vafa invariants of (general in complex moduli) Calabi-Yau 3-folds"
        "birational geometry of Calabi-Yaus"
        "creating and using CYDatabase's"
        "worked examples and example workflows using StringTorics"
///

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
      We use @TO reflexivePolytope@ to create the corresponding Macaulay2 object.
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
      For dimension 4 polytopes, the Kreuzer-Skarke database is generally
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
    (dump, ReflexivePolytope)
  Headline
    read and write a ReflexivePolytope to string
  Usage
    str = dump Q
  Inputs
    Q: ReflexivePolytope
  Outputs
    :String
      a string suitable for placing in a Macaulay2 database,
      designed so that reading it back intakes very little time and
      recomputation
  Description
    Text
      This dumps a ReflexivePolytope to a string (generally destined to be written to a file or
      database).  Only computed items are dumped.

      Here we see what happens if we dump a ReflexivePolytope that has just
      been constructed, and the same one after we have computed time-consuming
      features.
    Example
      verts = {{-1, -1, -1, 1}, {-1, -1, -1, 0}, {-1, -1, 0, 0}, {-1, 0, -1, 0},
          {0, -1, -1, 0}, {0, 0, -1, 1}, {1, 1, 2, -1}}
      Q = reflexivePolytope verts
      elapsedTime computeBasics Q
      dump Q
      Q1 = reflexivePolytope dump Q
      Q1 === Q
      elapsedTime computeBasics Q1
  SeeAlso
    (dump, CalabiYauInToric)
    (reflexivePolytope, String)
///

doc ///
  Key
    (reflexivePolytope, String)
  Headline
    read and write a ReflexivePolytope to string
  Usage
    Q = reflexivePolytope str
  Inputs
    str: String
      created using @TO (dump, ReflexivePolytope)@
  Outputs
    :ReflexivePolytope
      recreates the object without expensive recomputation
  Description
    Text
      This is the inverse to @TO (dump, ReflexivePolytope)@.
    Example
      verts = {{-1, -1, -1, 1}, {-1, -1, -1, 0}, {-1, -1, 0, 0}, {-1, 0, -1, 0},
          {0, -1, -1, 0}, {0, 0, -1, 1}, {1, 1, 2, -1}}
      Q = reflexivePolytope verts
      str1 = dump Q
      elapsedTime computeBasics Q
      str2 = dump Q
      Q1 = reflexivePolytope str1
      Q2 = reflexivePolytope str2
      dump Q1
      dump Q2
  Caveat
    Note that the name of the reconstructed object is the same as $Q$.  However, it is {\it not}
    the same object: it has different elements cached! (Strictly speaking it is the same, but
    realistically you might be surprised when it recomputes something!)
  SeeAlso
    (dump, CalabiYauInToric)
    (reflexivePolytope, String)
    (calabiYau, String, Function)
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
      Q = reflexivePolytope verts
      dim Q
  SeeAlso
    ReflexivePolytope
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
      This function returns the list of lattice points of a @TO
      ReflexivePolytope@.  The lattice points are sorted in some
      manner, such that the minimal face dimension of a lattice point
      is monotone increasing in this list.  As the origin is an
      interior lattice point (the only one too), its face dimension is the dimension of
      the polytope.

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
    (faceDimensions, ReflexivePolytope)
    (faceDimension, ReflexivePolytope, ZZ)
    (annotatedFaces, ReflexivePolytope)
///

doc ///
  Key
    faceDimensions
    (faceDimensions, ReflexivePolytope)
  Headline
    minimal face dimensions of the lattice points
  Usage
    faceDimensions Q
  Inputs
    Q:ReflexivePolytope
  Outputs
    :List
  Description
    Text
      Each lattice point of $Q$ lies on a unique smallest face of $Q$.
      The dimension of this face is its {\bf face dimension}.

      This function returns the list of the face dimension of each lattice point.
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
    Text
      Thus, there are 10 vertices, 4 lattice points which are interior to
      an edge, 2 which are interior to facets, and the origin.
  SeeAlso
    (faceDimension, ReflexivePolytope, ZZ)
    (latticePoints, ReflexivePolytope)
    (annotatedFaces, ReflexivePolytope)
///

doc ///
  Key
    faceDimension
    (faceDimension, ReflexivePolytope, ZZ)
  Headline
    minimal face dimension of a lattice point
  Usage
    faceDimension(Q, i)
  Inputs
    Q:ReflexivePolytope
    i:ZZ
  Outputs
    :ZZ
  Description
    Text
      Each lattice point of $Q$ lies on a unique smallest face of $Q$.
      The dimension of this face is its {\bf face dimension}.

      This function returns the list of the face dimension of the $i$-th lattice point.
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
  SeeAlso
    (faceDimensions, ReflexivePolytope)
    (latticePoints, ReflexivePolytope)
    (annotatedFaces, ReflexivePolytope)
///

doc ///
  Key
    annotatedFaces
    (annotatedFaces, ReflexivePolytope)
  Headline
    a list of faces of a reflexive polytope together with lattice point and face information
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
      netList annotatedFaces Q
    Text
      You can see how these are constructed.  First, we compute also the polar dual of Q
    Example
      P = polar Q
      # latticePoints P -- P has 29 lattice points
      transpose matrix latticePoints P
      # vertices Q -- Q has 10 vertices
      # vertices P -- Q has 11 maximal faces
      (matrix vertices P) * (transpose matrix vertices Q)


  SeeAlso
    (latticePoints, ReflexivePolytope)
    (faceDimensions, ReflexivePolytope)
    (polar, ReflexivePolytope)
///

doc ///
  Key
    findAllFRSTs
    (findAllFRSTs, ReflexivePolytope)
    (findAllFRSTs, Matrix)
    (findAllFRSTs, NormalToricVariety)
  Headline
    compute all fine, regular, star triangulations of a point set, using topcom
  Usage
    findAllFRSTs Q
  Inputs
    Q:{ReflexivePolytope, Polyhedron, Matrix, NormalToricVariety}
  Outputs
    :List
      of all triangulations
  Consequences
    Item 
      For ReflexivePolytope's, this information is cached in $Q$, under {\tt Q.cache.AllFRSTs}
  Description
    Text
      Given a matrix whose columns are points, this returns all fine (i.e. using all the columns), regular,
      star triangulations.  The returned triangulations are lists of lists of integer indices into the
      rays corresponding to the columns, vertices, etc. of the underlying object.

      Given a reflexive polytope (or a canonical polytope, that is, a lattice polytope whose only
      interior lattice point is the origin), this function uses all of the "rays" (lattice points not interior
      to facets or equal to the origin).
    Example
      tope = KSEntry "4 6  M:33 6 N:11 6 H:6,30 [-48] id:8
        1   0   0   3  -3  -3
        0   1   0   0  -2   2
        0   0   1   4  -6   2
        0   0   0   6  -6   0
        "
      Q = reflexivePolytope tope
      #rays Q == 10
      time Ts = findAllFRSTs Q -- takes very small amount of time
      #Ts
      V = normalToricVariety(rays Q, Ts#0)
      isWellDefined V
      isSimplicial V
      isProjective V
    Example
      A = transpose matrix rays Q
      findAllFRSTs A
  Caveat
    This function has no way to only return some of the triangulations.  For that, see
    @TO generateTriangulations@.
  SeeAlso
    findOneFRST
    reflexivePolytope
    normalToricVariety
///

doc ///
  Key
    "normal form of a reflexive polytope"
  Headline
    algorithm for computing normal forms of polytopes
  Description
    Text
      The algorithm is described in kreuzer and skarke, arxiv 9805190v1.

      The basic idea: Find a normal form for a lattice polytope $\Delta$ in $\RR^d$, having
      $n$ vertices, and let's assume that the polytope $\Delta$ has dimension $d$ as well.
      Let $A$ be the $d \times n$ integer matrix whose columns are these vertices.

      The goal: find a normal form for $A$, such that two such lattice polytopes
      are the same up to an invertible integer change of basis
      if and only if they have the same normal form.

      The group $GL_d(\ZZ)$ acts on matrices $A$ on the left, and the
      symmetric group $S_n$ on the right.  We want to find a normal form under this
      group.

      Step 1: Compute the matrix $B$ ($m \times d$) whose rows are the vertices
      of the polar dual of $\Delta$.  This is determined up to a permutation in
      $S_m$, given $A$.

      Step 2: Compute the vertex pairing matrix $N = B \cdot A$. This matrix is unqiue
      up to a permutation of rows and a permutation of columns. (Does this require $\Delta$
      is reflexive?  Probably not...)

      Step 3: Choose the lexicographically largest matrix $N$ under swapping rows
      and columns.  What does this actually mean?  In any case, we get
      a normal form for the vertex pairing matrix.

      Step 4: Now we need to start changing $A$ itself (we still have the group
      $GL_d(\ZZ)$ to work with, and perhaps a subgroup of the symmetric group (s).
      After swapping columns of $A$ to match the vertex pairing matrix,
      we use Hermite reduction to make positive diagonal entries on $A$, and each
      upper triangular entry above is minimal.

      Step 4: Actually, we do step 4 for each possible choice of permutation of
      the vertices of $\Delta$, and we choose the lexicographically smallest one
      of all, and that is our normal form.  This could be a problem if the subgroup is
      large!

      Note: we have to define more carefully what largest in lexicographic order
      means.  We need to order two matrices $A$, $A'$ of the same size.  Do we use
      the usual lex order first via columns or via rows?  It might matter...
    Text
      Let's try this on an example reflexive polytope (first).  We will first see what PALP
      returns, then we will do it "by hand" ourselves.
    Example
      needsPackage "PALPInterface"
      M = matrix KSEntry "4 15  M:27 15 N:16 12 H:11,22 [-22] id:100
        1   0   0   0   0   1   0   1   1   0   0   0   0  -2  -1
        0   1   0   0   0   1   1   0  -1   0  -1  -2  -2   0  -1
        0   0   1   0  -1  -1   0  -1   0  -1  -1   0   2   1  -1
        0   0   0   1   0   0  -1   0   0  -1  -1  -1   1   2  -1
        "
      Q = reflexivePolytope M
      P = polar Q
      A = transpose matrix vertices Q
      B = matrix vertices P
      B*A
      tally flatten entries (B*A)

      S = matrix {{0, 0, -3, 2}, {-5, 3, 0, 0}, {-3, -6, 1, 0}, {-3, -2, -1, 1}}  
      Q2 = reflexivePolytope (S*M)
      P2 = polar Q2
      A2 = transpose matrix vertices Q2
      B2 = matrix vertices P2
      B2*A2
      tally flatten entries (B2*A2)
      normalForm M
  SeeAlso
///

-*
VPM NF (v=15 f=12):
  2  0  2  1  0  0  0  1  3  0  1  3  5  0  0 =>  5  4  4  3  2  1  1  1  0  0  0  0  0  0  0
  1  2  2  0  0  1  3  0  0  1  0  0  0  0  0 =>  4  3  3  2  2  1  1  0  2  1  1  0  0  0  0
  2  1  0  2  2  3  0  3  2  1  1  0  0  0  0 =>  3  3  3  2  2  2  2  0  0  1  1  1  1  0  0
  1  1  2  1  0  0  1  0  1  0  0  1  3  2  0 =>  3  2  2  2  1  0  0  2  1  0  0  0  0  1  1
  1  1  1  2  1  1  0  1  1  0  0  0  2  3  0 =>  2  1  0  1  2  1  0  0  4  2  0  1  0  1  0
  0  0  0  0  2  0  1  1  1  3  4  4  0  0  5 =>  1  1  2  1  0  0  1  2  0  0  2  0  1  1  2
  1  1  2  0  0  0  2  0  1  1  1  2  2  0  1 =>  1  0  0  0  1  0  0  0  5  2  2  0  0  1  1
  1  0  0  1  2  1  0  2  2  2  3  3  1  0  3 =>  0  1  3  0  0  1  3  0  0  1  5  0  2  0  2
  0  1  1  0  1  0  2  0  0  2  2  2  0  1  3 =>  0  1  0  1  2  3  2  0  0  2  0  3  2  1  0
  0  1  1  2  1  0  0  0  0  0  0  0  2  5  1 =>  0  0  1  0  0  0  1  1  2  1  3  0  1  1  2
  0  0  0  1  2  0  0  1  1  2  3  3  1  2  4 =>  0  0  0  1  0  0  0  3  0  0  0  1  1  2  2
  0  1  0  2  2  1  0  1  0  1  1  0  0  4  2 =>  0  0  0  0  1  1  1  0  3  2  2  1  1  1  1

Poly NF try[0]:   C=0123456789abcde
   1   0   0   0   0   1   0   1   1   0   0   0   0  -2  -1 =>   1   0   0   0   0   1   0   1   1   0   0   0   0  -2  -1
   0   1   0   0   0   1   1   0  -1   0  -1  -2  -2   0  -1 =>   0   1   0   0   0   1   1   0  -1   0  -1  -2  -2   0  -1
   0   0   1   0  -1  -1   0  -1   0  -1  -1   0   2   1  -1 =>   0   0   1   0  -1  -1   0  -1   0  -1  -1   0   2   1  -1
   0   0   0   1   0   0  -1   0   0  -1  -1  -1   1   2  -1 =>   0   0   0   1   0   0  -1   0   0  -1  -1  -1   1   2  -1

Poly NF:  NormalForm=try[0]  #Sym(VPM)=1  #Sym(Poly)=1

V_perm made by Poly_Sym (order refers to VertNumList):
0123456789abcde
4 15  Normal form of vertices of P    perm=0123456789abcde
   1   0   0   0   0   1   0   1   1   0   0   0   0  -2  -1
   0   1   0   0   0   1   1   0  -1   0  -1  -2  -2   0  -1
   0   0   1   0  -1  -1   0  -1   0  -1  -1   0   2   1  -1
   0   0   0   1   0   0  -1   0   0  -1  -1  -1   1   2  -1
*-

///
-- how many different "tally" do all the h11=3 ones have.
-- make the list {#vertices, #facets, tally of values}
R = ZZ[a,b,c]
DB3 = StringTorics#"auxiliary files" | "cy3-h11-3.dbm"
(Qs, Xs) = readCYDatabase(DB3, Ring => R);
for k in keys Qs list (
    Q := Qs#k;
    P := polar Q;
    A := transpose matrix vertices Q;
    B := matrix vertices P;
    {# vertices Q, # vertices P, tally flatten entries (B*A)}
    )

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




-- This one is still a mess! (April 2025)
-- Problems: isProjective V -- takes how long?!
--           generateTriangulations is strange and hard to use?
///
  Key
    "example: generating some triangulations"
  Headline
    example of finding triangulations
  Description
    Text
      Let's start with one that we can compute with topcom all of the
      triangulations, say at h11=5.
    Example
      tope = KSEntry "4 10  M:38 10 N:10 7 H:5,33 [-56] id:6
        1   1   1   0   1  -2  -1   0   1   1
        0   2   0   0  -1   2  -3   4  -1   2
        0   0   2   0  -1  -2   1   0   3   2
        0   0   0   1  -1   1  -1   1  -1   0
        "
      Q = reflexivePolytope tope
      Ts = findAllFRSTs Q; 
      assert(#Ts == 9) -- 9 of these triangulations
      Vs = for t in Ts list normalToricVariety(rays Q, t);
      Vs/isWellDefined -- all is good
      Vs/isComplete
      Vs/isProjective
      findAllCYs Q -- all are equivalent?? looks like it.
      for t in Ts list restrictTriangulation(2, Q, t)
      assert(# unique oo == 1)-- yes, only one of them!
    Text
      What about generateTriangulations?? 
    Example
      A = transpose matrix latticePoints Q
      T's = generateTriangulations A
      #T's -- 112
      stars = select(T's, t -> isStar t)
      #stars == 9
      newset = set(stars/max/(t -> (t/(t0 -> drop(t0, -1)))))
      oldset = set Ts
      oldset == newset -- yes!
    Text
      We start with a reflexive polytope that has many triangulations.  We want to find,
      say, 100 of these.
    Example
      tope = KSEntry "4 9  M:13 9 N:23 10 H:20,11 [18] id:6
        1    0    0    0    0   -1   -1    0    2
        0    1    0    0   -1    2    1   -2    0
        0    0    1    0    1    0    1   -1   -1
        0    0    0    1    1   -1   -1    1   -1
        "
      Q = reflexivePolytope tope
      tri = findOneFRST Q
      V = normalToricVariety(rays Q, tri)
      isWellDefined V
      isComplete V
      elapsedTime isProjective V
      A = transpose matrix latticePoints Q
      faceDimensions Q -- this should match latticePoints...
      (pts, tri) = regularStarTriangulation(2, polytope Q)
      #pts
      pts == rays Q -- true!! why?
      V = normalToricVariety(rays Q, tri)
      time assert isWellDefined V -- not cheap at h11=20... (8 sec)
      isComplete V
      time assert isProjective V -- takes some time: 
      tri1 = triangulation(transpose matrix rays Q, tri)
      isFine tri1
      isRegularTriangulation tri1 -- code in Topcom crashes...
      
      first1 = regularFineTriangulation A
      tris = generateTriangulations(first1, Limit => 1000); -- hmm, gives 15...
      stars = select(tris, isStar)
      V = normalToricVariety(rays Q, max tris#0)
      isWellDefined V -- not a good sign...
      generateTriangulations(A, Limit => 100) -- hmm, gives 102...
      -- recall, these will be fine
      isFavorable Q -- not favorable
      generateTriangulations(A, Homogenize => false, Limit => 10) -- hmm, gives 15...
      for i from 0 to #latticePoints Q - 1 list faceDimension(Q, i)
      #rays Q === #latticePoints Q - 3
      sort vertices Q == sort verts
      #rays Q
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


-*
o120 = {0 => (dump, ReflexivePolytope)                                     }
       {1 => (rays, ReflexivePolytope)                                     }
       {4 => (label, ReflexivePolytope)                                    }
       {5 => (addToCYDatabase, String, ReflexivePolytope)                  }
       {11 => (vertices, ReflexivePolytope)                                }
       {12 => (latticePoints, ReflexivePolytope)                           }
       {13 => (faceDimensions, ReflexivePolytope)                          }
       {14 => (faceDimension, ReflexivePolytope, ZZ)                       }
       {15 => (dim, ReflexivePolytope)                                     }
       {16 => (polytope, ReflexivePolytope, String)                        }
       {17 => (polytope, ReflexivePolytope)                                }
       {18 => (polar, ReflexivePolytope)                                   }
       {19 => (annotatedFaces, ReflexivePolytope)                          }
       {20 => (findTwoFaceInteriorDivisors, ReflexivePolytope)             }
       {21 => (degrees, ReflexivePolytope)                                 }
       {22 => (basisIndices, ReflexivePolytope)                            }
       {23 => (isFavorable, ReflexivePolytope)                             }
       {24 => (hh, Sequence, ReflexivePolytope)                            }
       {25 => (automorphisms, ReflexivePolytope)                           }
       {26 => (automorphismsAsPermutations, ReflexivePolytope)             }
       {29 => (findAllFRSTs, ReflexivePolytope)                            }
       {30 => (findAllFRVTs, ReflexivePolytope)                            }
       {31 => (findOneFRST, ReflexivePolytope)                             }
       {32 => (restrictTriangulation, ZZ, ReflexivePolytope, List)         }
       {33 => (partitionFRSTsByDFaceEquivalence, ZZ, ReflexivePolytope)    }

       {2 => (calabiYau, ReflexivePolytope, List)                          }
       {3 => (makeCY, ReflexivePolytope)                                   }
       {34 => (makeCYs, ReflexivePolytope)                                 }
       {35 => (findAllCYs, ReflexivePolytope)                              }
       {27 => (isTriangulationOfPolytope, ReflexivePolytope, List)         }
       {28 => (isTriangulationOfPolytope, ReflexivePolytope, Triangulation)}
       {6 => (calabiYau, Database, ReflexivePolytope, Sequence)            }
       {7 => (calabiYau, String, ReflexivePolytope, Sequence)              }

       {8 => (expression, ReflexivePolytope)                               }
       {9 => (describe, ReflexivePolytope)                                 }
       {10 => (net, ReflexivePolytope)                                     }
*-       

-- TODO? change ReflexivePolytope, CYPolytope to CanonicalPolytope?

doc ///
  Key
    CYPolytope
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
      A {\tt CYPolytope} can be constructed from vertices of a reflexive polytope.
    Text
      For example, let's start with a dimension 3 reflexive polytope: the cube.
      We use @TO cyPolytope@ to create the corresponding Macaulay2 object.
    Example
      verts = {
          {-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
          {-1, -1, 1}, {1, -1, 1}, {-1, 1, 1}, {1, 1, 1}}
      Q = cyPolytope verts
      dim Q
    Text
      Given a {\tt CYPolytope}, the method {\tt rays} returns a list of the boundary lattice points
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
    Example
      P2 = polytope Q
      isReflexive P2
      vertices P2
      matrix {latticePoints P2}
      entries transpose oo
      latticePoints Q
    Text
      Notice that latticePoints of the {\tt Polyhedra} object, {\tt P2} are in a different order
      than the lattice points from the {\tt CYPolytope} $Q$  All indices for a CYPolytope
      object (faces, maximal cones, etc) {\it all} refer to the order of rays/lattice points from $Q$.
  SeeAlso
      CalabiYauInToric
      cyPolytope
///

doc ///
  Key
    "facilities available for working with triangulations"
  Headline
    facilities available for working with triangulations
  Description
    Text
      In this introduction, we consider triangulations of the following square
      in the plane.
    Example
      square = transpose matrix{{1,1},{-1,1},{-1,-1},{1,-1}}
      regularFineTriangulation square
      assert(# allTriangulations square == 2)
    Text
      Now consider all of the lattice points of the square.
    Example
      P = convexHull square
      LP = latticePointList P
      sq9 = transpose matrix LP
    Text
      We could have entered sq9 by hand as so:
    Example
      sq9 = matrix {{-1, -1, 1, 1, -1, 0, 0, 1, 0},
                    {-1, 1, -1, 1, 0, -1, 1, 0, 0}}
    Text
      We first show some functions from Topcom that are useful.
    Example
      t1 = regularFineTriangulation sq9
      regularTriangulationWeights t1
      fineStarTriangulation(sq9, max t1)
      delaunaySubdivision sq9 -- not a triangulation (4 squares).
      orientedCircuits sq9 -- many of these are not useful when considering only fine triangulations.
    Text
      Let's generate all of the triangulations of the square.
      Really, we want all triangulations which are fine (involve all the lattice points)
      and are star (involve the origin), and are regular.
    Example
      regularFineStarTriangulation sq9 -- leaves out 8, the index of the origin in sq9.
      Ts = allTriangulations sq9;
      #Ts
      Ts = Ts/max;
      # select(Ts, t -> isFine(sq9,t))
      # select(Ts, t -> isStar(sq9,t))
      # select(Ts, t -> isStar(sq9,t) and isFine(sq9,t))
      # select(Ts, t -> isFine(sq9,t) and isRegularTriangulation(sq9,t))
    Text
      Regular triangulations and subdivisions can be computed.
      In this example, we take the first 5 fine triangulations found above,
      find weights (they are all regular triangulations) giving these triangulations,
      then reconstruct the triangulation using @TO regularSubdivision@.
      These are the same as the triangulations we started with.
    Example
      fineT = take(select(Ts, isFine_sq9), 5)
      wts = for t in fineT list regularTriangulationWeights(sq9, t)
      fineT2 = for w in wts list regularSubdivision(sq9, matrix{w})
      fineT == fineT2
    Text
      We might want to check that these are indeed triangulations.
      I am not completely convinced that @TO topcomIsTriangulation@ always gives a
      correct answer, so we also implement a slower routine @TO naiveIsTriangulation@.
    Example
      starT = first select(Ts, t -> isStar(sq9,t))
      naiveIsTriangulation(sq9, starT)
      topcomIsTriangulation(sq9, starT)
      notSq9 = matrix {{-1, -1, 1, 1, -1, 0, 0, 2, 0},
                    {-1, 1, -1, 1, 0, -1, 1, 0, 0}}
      naiveIsTriangulation(notSq9, starT)
      topcomIsTriangulation(notSq9, starT)
      debug Triangulations
      isTriangulation(notSq9, starT)
    Text
      All regular triangulations fit into a polytope, whose vertices are the
      GKZ volume vectors (for each lattice point, consider the sum of the volumes
      of the simplices containing the point as a vertex).  This gives a vector in
      $\Z^d$, where $d$ is the number of lattice points, which is computed by the method
      @TO volumeVector@.
    Example
      volume convexHull sq9
      tri = fineT_0
      for f in tri list volume convexHull(sq9_f)
      sum oo == volume convexHull sq9
      volumeVector(sq9, tri)
      volumeVector(sq9, starT)
    Text
      Sometimes we want to generate only some of the triangulations, as there can be
      a huge number of them.  Unfortunately, the topcom functions do not allow this
      functionality.  Instead, use @TO generateTriangulations@.  Note that this function
      only generates fine triangulations.
    Example
      T4 = generateTriangulations(sq9, Limit => 100);
      T3 = select(Ts, t -> isFine(sq9,t));
      assert(set (T4/max) === set T3)
  SeeAlso
    generateTriangulations
    "Topcom::allTriangulations"
///
