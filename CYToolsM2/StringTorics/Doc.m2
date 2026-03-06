-- missing doc:
-- CYPolytope:
--  annotatedFaces
--  basisIndices
--  isFavorable
--  polar
--  degrees

-- ReflexivePolytope
--  polar
--  makeCY
--  findAllCYs
--  (makeCY, makeAllCYs)?

-- regularFineStarTriangulation

-- Database doc
-- seems a bunch are missing

-- Cohomology

-- OUT OF DATE, to be removed.
doc ///
    Key
        StringTorics
    Headline
        toric variety functions useful for investigations in string theory
    Description
        Text
            The purpose of this package is to provide facilities to
            explore the geometry of Calabi-Yau 3-folds arising as
            hypersurfaces in toric varieties, and the resulting
            applications in string theory.
        Text
            There are a large number of Calabi-Yau 3-folds which lie
            on a toric variety.  The main example is the construction
            of Batyrev, taking anti-canonical hypersurfaces in
            simplicial projective toric varieties given a
            triangulation of a reflexive polytope.  However, there are
            others too: complete intersections in products of
            projective spaces or other toric varieties; triangulations
            of the non-simplicial fan which do not induce a
            triangulation of a reflexive polytope: these are often
            called {\it vex} triangulations.
        Text
            This package allows one to work with these Calabi-Yau
            varieties.  Currently, the best supported are the Batyrev
            hypersurfaces.
        Text
            Although I wrote this code, much of the work to determine
            algorithms and use them in research was done in collaboration
            with Liam McAllister (Cornell physics), and our collaborators
            (Cody Long, Jakob Moritz, Elijah Sheridan, Naomi Gendler,
            Manki Kim, Andreas Braun, Andres Rios-Tascon).
    SeeAlso
        "working with reflexive polytopes"
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


-- REPLACEMENT for main doc node.
doc ///
   Key
     "StringTorics2"
   Headline
     toric variety functions useful for investigations in string theory
   Description
    Text
      The purpose of this package is to provide facilities to explore the
      geometry of Calabi-Yau 3-folds arising as hypersurfaces in toric varieties,
      and the resulting applications in string theory.
      

      Below, one sees some of the computations that are made available in this package.
      
      This package uses the software packages TOPCOM, CohomCalg, and PALP, together
      with facilities already present in Macaulay2, to provide the following 
      functionality.
    Text
      @SUBSECTION "Reflexive polytopes"@
    Text
      @SUBSECTION "Triangulations"@
    Text
      @SUBSECTION "Examples of use"@
    Text
      @UL {
          {TO "Example use", ", a first example showing basic usage of this package"}
          }@
    Text
      @SUBSECTION "Reflexive polytopes"@
    Text
      In this package, a key type is @TO ReflexivePolytope@ (this supercedes @TO CYPolytope@.
      We might still change the name to CanonicalPolytope).  Objects of this class
      represent a reflexive polytope.  However, it also computes and caches information
      about the Batyrev
      Calabi-Yau hypersurfaces that is independent of the FRST (i.e. fine regular star triangulation
      of this reflexive polytope).  This includes the GLSM charge matrix, and a description of
      the generators of the Picard group of the corresponding Calabi-Yau varieties.
    Text
      @UL {
          TO ReflexivePolytope,
          TO (latticePoints, ReflexivePolytope),
          TO (annotatedFaces, ReflexivePolytope),
          --TO (isFavorable, ReflexivePolytope),
          TO (polar, ReflexivePolytope)
          --TO (degrees, ReflexivePolytope)
          --TO basisIndices,
          }@
    Text
      @SUBSECTION "Batyrev Calabi-Yau hypersurfaces"@
    Text
      The data for a Calabi-Yau hypersurface is a reflexive polytope, as above, and a triangulation
      of the set of all lattice points not interior to facets, that is an "FRST" (fine: it uses all of these points,
      regular: it is a regular triangulation, and star: all simplices contain the origin).  Basic construction
      and use of these include the following.  For more involved algorithms, see the sections below.
    Text
      @UL {
          TO CalabiYauInToric,
          TO calabiYau,
          TO makeCY,
          TO makeCYs,
          TO (label, CalabiYauInToric),
          TO (max, CalabiYauInToric),
          TO (rays, CalabiYauInToric),
          TO (normalToricVariety, CalabiYauInToric),
          TO (picardRing, CalabiYauInToric),
          TO (dump, CalabiYauInToric)
          }@
    Text
      @SUBSECTION "Routines to access the Kreuzer-Skarke database"@
    Text
      @UL {
          {TO "kreuzerSkarke"},
          {TO "kreuzerSkarkeDim3"}
          }@
    Text
      @SUBSECTION "Triangulations"@
    Text
      @UL {
          TO "facilities available for working with triangulations",
          TO findAllFRSTs,
          TO findAllFRVTs,
          TO findOneFRST,
          "More general triangulation features",
          TO regularFineStarTriangulation,
          TO allTriangulations,
          TO generateTriangulations,
          TO "example: generating some triangulations"
          }@
    Text
      @SUBSECTION "Intersection numbers and intersection rings"@
    Text
      @UL {
          TO (intersectionNumbers, CalabiYauInToric)
          }@
    Text
      @SUBSECTION "Topological equivalence of Calabi-Yau 3-folds"@
    Text
      @UL {
          }@
    Text
      @SUBSECTION "Gopakumar-Vafa invariants of (general in complex moduli) Calabi-Yau 3-folds"@
    Text
      @UL {
          }@
    Text
      @SUBSECTION "Birational geometry of Calabi-Yaus"@
    Text
      @UL {
          }@
    Text
      @SUBSECTION "Calabi Yau hypersurfaces in toric varieties"@
    Text
      @UL {
          {TO "CalabiYauInToric"},
          {TO "makeCY"},
          {TO "findAllCYs"}
          }@

    Text
      @SUBSECTION "Creating and using CYDatabase's"@
      
      A CYDatabase is a file which contains precomputed data about a collection of
      (reflexive) 4D-polytopes and the resulting Calabi Yau hypersurfaces.
    Text
      @UL {
          TO "Accesing a CYDatabase file",
          TO "Creating a CYDatabase file",
          TO addToCYDatabase,
          TO combineCYDatabases,
          TO readCYDatabase,
          TO readCYPolytopes,
          TO readCYs
          }@
    Text
      @SUBSECTION "Cohomology"@
    Text
      @UL {
          TO "cohomology of line bundles on toric varieties",
          TO cohomologyVector,
          TO toricOrthants,
          TO cohomologyBasis,
          TO cohomologyMatrix,
          TO cohomologyMatrixRank
          }@
   Caveat
   SeeAlso
     "installing StringTorics"
///

-- CYPolytope doc node moved to DocCYPolytope.m2

-- Creation of a CYPolytope

-- the following function is not ready yet, maybe we don't even want it??
///
  Key
    cyPolytopeFromRays
  Headline
    the most basic constructor of a CYPolytope 
  Usage
    Q = cyPolytopeFromRays L
  Inputs
    L:List
  Outputs
    Q:CYPolytope
  Description
    Text
      This is a light-weight constructor once the rays are known.
    Text
      Let's start with a reflexive 4D polytope corresponding to a Calabi-Yau
      3-fold with $h^{1,1} = 3$.
    Example
      topes = kreuzerSkarke 3;
      tope = topes_30
      A = matrix tope
      P = convexHull A
      P2 = polar P
      allrays = matrix{select(latticePoints P2, lp -> lp != 0)}
      allrays = transpose entries allrays
      Q = cyPolytopeFromRays allrays
      assert(rays Q === allrays)
      vertices P
      (transpose vertices P) * transpose matrix allrays
      convexHull transpose matrix allrays
      vertices oo

      Q = cyPolytope tope
      A = transpose matrix rays Q
      vertices convexHull A
      latticePoints convexHull A
      transpose matrix degrees Q
      basisIndices Q
      automorphisms Q
      hh^(1,1) Q
      hh^(1,2) Q
      dim Q
      A
      netList annotatedFaces Q
    Text
      We want to make sure we can take rays of a fan, in order, that someone gives us,
      and *not* change the order of the generators...!
///

///
  Key
    cyPolytope
    (cyPolytope, Matrix)
    (cyPolytope, List)
    (cyPolytope, Polyhedron)
    (cyPolytope, KSEntry)
  Headline
    construction of a CYPolytope (reflexive polytope)
  Usage
    Q = cyPolytope A
  Inputs
    A:{Matrix,List,Polyhedron,KSEntry}
  Outputs
    Q:CYPolytope
  Description
    Text
      A CYPolytope contains the data associated to a reflexive
      polytope.  One is most easily constructed from vertices for this
      polytope (given either an integer matrix whose columns are these
      vertices, or a list of these vertices, or a @TO "Polyhedron"@
      object from the @TO "Polyhedra"@ package.

      This function does {\it not} check that the result is a reflexive polytope.
      You must do this yourself if you are want to check.
      Use "TO (isWellDefined, CYPolytope)" or "TO (isReflexive, CYPolytope)" to check this.
      
      One direct way to construct this reflexive polytope is to give it a
      list of (integer) vertices or various points.  The function takes the
      convex hull of these points, and returns the reflexive polytope.
    Example
      verts = {{-1, -1, -1, 1}, {-1, -1, -1, 0}, {-1, -1, 0, 0}, {-1, 0, -1, 0},
          {0, -1, -1, 0}, {0, 0, -1, 1}, {1, 1, 2, -1}}
      Q = cyPolytope verts
      vertices Q
      P2 = polytope Q
      isReflexive P2
    Text
      Instead of a list of points, you may also give a matrix whose columns are
      integral points.  The convex hull of these is constructed.
    Example
      pts = transpose matrix latticePoints Q
      Q1 = cyPolytope pts
      vertices Q1 === vertices Q
    Text
      If you call cyPolytope with a polytope as its argument (which happens to be
      reflexive), you get a CYPolytope referring to the same polytope.
    Example
      P2 = polytope Q
      Q2 = cyPolytope P2
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
      Q3 = cyPolytope tope
      vertices Q3 == vertices Q
      (hh^(1,1) Q, hh^(1,2) Q)
  SeeAlso
    (vertices, CYPolytope)
    (latticePoints, CYPolytope)
--    (annotatedFaces, CYPolytope)
///

///
  Key
    annotatedFaces
--    (annotatedFaces, CYPolytope)
    (annotatedFaces, Polyhedron)
  Headline
    a list of faces of a reflexive polytope together with lattice point information
  Usage
    annotatedFaces Q
  Inputs
    Q:{CYPolytope, Polyhedron}
  Outputs
    :List
      each entry is a list containing: the dimension of the face, the indices of the
      vertices, the indices of all (boundary) lattice points in the face, the number
      of interior points in the face, and the number of interior points in the dual face
  Description
    Text
      verts = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, -1}, {0, -1, 0}, {-1, 0, 0}, {-1, 1, 0}}
      QN = reflexivePolytope verts
      netList annotatedFaces QN
      vertices polytope(QN, "N")
      vertices polytope(QN, "M")
      latticePoints polytope(QN, "N")
      latticePoints polytope(QN, "M")
      polar QN
      rays oo
      netList annotatedFaces QN
      netList annotatedFaces polar QN
      PN = convexHull transpose matrix verts
      PM = polar PN
      vertices PM
      
      isReflexive PN
      vertices PN -- notice these are in a different order

      netList annotatedFaces QN
      rays smoothFanoToricVariety(3, 12) -- this is how we obtained these vertices.
    Example
  Caveat
  SeeAlso
    latticePointList
///



///
  Key
  Headline
  Usage
  Inputs
  Outputs
  Consequences
    Item
  Description
    Text
    Example
  Caveat
  SeeAlso
///

///
  Key
    cyPolytope
  Headline
    create a reflexive polytope pair
  Usage
    cyPolytope ks
    cyPolytope vertexlist
    cyPolytope m
    cyPolytope P
    cyPolytope str
  Inputs
    ks:KSEntry
    vertexList:List
    m:Matrix
      the vertices are the columns
    Q:CYPolytope
      or @ofClass Polyhedron@
    ID => ZZ
      :ZZ
        a label for this polytope (TODO: is this a string or integer?)
  Outputs
    :CYPolytope
      
  Description
    Text
      topes = kreuzerSkarke(3, Limit => 50);
      topes_40
      matrix topes_40
      Q = cyPolytope topes_40 -- this is a polytope
      rays Q
      Q1 = cyPolytope matrix topes_40
      Q2 = cyPolytope rays Q1
      verts = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, -1}, {0, -1, 0}, {-1, 0, 0}, {-1, 1, 0}}
      QN = cyPolytope verts
      netList annotatedFaces QN
      vertices polytope(QN, "N")
      vertices polytope(QN, "M")
      latticePoints polytope(QN, "N")
      latticePoints polytope(QN, "M")
      polar QN
      rays oo
      netList annotatedFaces QN
      netList annotatedFaces polar QN
      PN = convexHull transpose matrix verts
      PM = polar PN
      vertices PM
      
      isReflexive P
      vertices P -- notice these are in a different order

      netList annotatedFaces Q
      rays smoothFanoToricVariety(3, 12) -- this is how we obtained these vertices.
    Example
  Caveat
  SeeAlso
///

-*
    Text
      @SUBSECTION "Additional polyhedral functions"@
    Text
      @UL {
          {}
          }@

    Text
      @UL {
          {}
          }@
    Text
      @UL {
          {}
          }@
*-

doc ///
   Key
     "Example use"
   Headline
     An example of using the functionality of the package
   Description
    Text
      Let's analyze one particular toric variety, and the corresponding 
      Calabi-Yau hypersurface.

      First, take an example from the Kreuzer-Skarke database.  Note: you need to be online
      in order for this to work!
    Text
      We don't run the code below here, as it requires a network connection.
      But you should try it.
    Pre
      polytopes = kreuzerSkarke(5, 57, Limit=>200, Access => "wget");
      #polytopes == 197
    Text
      There are 197 reflexive polytopes in 4 dimensions, whose
      (resolution of an) anti-canonical divisor X (a Calabi-Yau),
      has $h^{1,1}(X) = 5$, and $h^{1,2}(X) = 57$.
      
      Let's consider the 11th one on this list.
    Example
      polytopes_10
      ks = KSEntry "4 10  M:56 10 N:11 7 H:5,57 [-104] id:10
        1   1   1   1   0   0   0  -1  -1  -1
        0   2   0   0  -1  -1   3   3  -2  -2
        0   0   2   0   3  -1  -1  -2  -2   3
        0   0   0   2  -1   3  -1  -2   3  -2
        "
      A = matrix ks
    Text
      The Calabi-Yau is (the resolution of) an anti-canonical hypersurface
      in the 4 dimensional projective toric variety whose polytope is $P_1$
      in the $M$ lattice.
      
      This polytope has 10 vertices, 7 facets, and 56 lattice points: 
      the origin and 55 on the boundary.
    Example
      P1 = convexHull A
      fVector P1
      vertices P1
      matrix{latticePoints P1}

      # latticePointList P1
      vertexList P1
      transpose matrix latticePointList P1
    Text
      Now let's take the dual polytope (in the $N$-lattice)
    Example
      P2 = polar P1
    Text
      The following comand returns the list of non-zero lattice points, and
      the vertices appear first on this list.
    Example
      latticePointList P2
      transpose matrix latticePointList P2
      vertexMatrix P2
      V0 = normalToricVariety normalFan P1
      isSimplicial V0
      isSmooth V0
      V = reflexiveToSimplicialToricVariety P1
      isSimplicial V
      isSmooth V
      rays V == (latticePointList P2)_{0..#rays V-1}
    Text
      This is the simplicial toric variety we want!  Actually, any other
      simplicial toric variety will work too (i.e. any other fine star regular
      triangulation of the point configuration {\tt latticePointList P2}.
    Text
      The Calabi-Yau $X$ is an anti-canonical divisor in $V$ (the same as the
          pull-back of an anti-canonical divisor in $V_0$, which is ample, to $V$
          (the result is no longer ample).
      Batryev (based on Danilov-Khovanskii) has provided formulae for the Hodge numbers
      of such a hypersurface.  The formula for $h^{1,1}(X)$ and $h^{2,1}(X)$, for $X$ 
      the desingularization of a hypersurface in a Fano toric variety of dimension 4, 
      depends only on the polytope.
    Example
      h11OfCY P1
      h21OfCY P1
    Text
      The reflexive polytope is called favorable if $H^{1,1}(X)$ is generated by the
      images of toric divisors on $V$.
    Example
      isFavorable P1
    Text
      Intersection ring of $X$
    Text
      First, let's work on $X$, which is a sufficiently generic hypersurface of $V$
      in the divisor class $|-K_V|$.
    Example
      KV = toricDivisor V
      X = completeIntersection(V, {-KV})
    Text
      The Hodge diamond of $X$, as long as $X$ is sufficiently generic,
      is determined via Danilov-Khovanskii (reference: ...).  However, 
      the following works whenever $X$ is a smooth hypersurface of $V$.
    Example
      hodgeDiamond X
    Text
      Thus $h^{1,1}(X) = 5$ and $h^{1,2}(X) = 57$.
      
      We can do intersection theory on $X$ and $V$.  Note that $V$ is not smooth, but
      since it is simplicial, as long as we work over the rationals, not the integers,
      the intersection ring is well-defined.  Since $X$ is smooth (or assumed to be),
      the intersection ring of $X$ is well-defined.
      
      The way we do all of this is we first create the abstract variety of a point.
      We include some parameters here (we will use them below).  Then, we define
      abstract varieties for $V$ and $X$.  These are designed to work well with the
      intersection theory implemented in Macaulay2 (the package @TO "Schubert2"@).
      An abstract variety contains information about its intersection ring (or, the numerical
      intersection ring), the chern class of its tangent bundle, and how to integrate 
      a cycle on the variety.    
    Example
      pt = base(symbol a, symbol b, symbol c, symbol d, symbol e)
      Va = abstractVariety(V, pt)
      Xa = abstractVariety(X, pt)
      IX = intersectionRing Xa
    Text
      For example, the (restrictions to $X$) of the 9 toric prime divisors:
    Example
      gens IX
      integral(t_2*t_3*t_4)
      integral(t_2^3)
    Text
      The advantage of putting parameters into the base, is that we can generate
      formulas.  For example, the euler characteristic of $L = OO_X(aD_4+bD_5+cD_6+dD_7+eD_8)$
      ($L$ is a line bundle with the given element as first Chern class. Note
      that $D_4, \ldots, D_8$ generate the intersection ring) can be found using
      the following code, which invokes Hirzebruch-Riemann-Roch to find the euler
      characteristic.
    Example
      L = a*t_4 + b*t_5 + c*t_6 + d*t_7 + e*t_8
      F1 = chi OO(L)
      6*(F1-1)
   SeeAlso
///

-- "Example: (3,3) hypersurface in P2 x P2" doc node moved to DocCI.m2



///
-- test: FRVTs -- fails some time...

findAllFRVTs CYPolytope := List => Q -> (
  cQ := chirotope(transpose matrix rays Q, Homogenize => false);
  "foo-topcomfoo.in" << toString cQ << endl << close;
  ctris := for L in lines get ("!chiro2allfinetriangs <foo-topcomfoo.in") list (
      matches := regex("\\{\\{[0-9,\\{\\}]*", L);
      if matches === null or #matches != 1 then error "my logic is wrong";
      value substring(matches_0, L)
      );
  ctris)

findAllFRVTs CYPolytope := List => Q -> (
    alltris := topcomAllTriangulations(transpose matrix rays Q, ConnectedToRegular => false, Fine => true, Homogenize => false, RegularOnly => false);
    -- now we need to choose the ones that are regular...
    select(alltris, t -> isProjective normalToricVariety(rays Q, t))
    )

-- Returns: all regular fine simplicial fans of the rays of Q
findAllFRVTs CYPolytope := List => Q -> (
    A := transpose matrix rays Q;
    findAllSimplicialFans A
    )

findAllFRVTs Qs#0
for k in sort keys Qs list elapsedTime findAllFRVTs Qs#k -- these work for h11=2 examples...

-- h11=3
      topes = kreuzerSkarke(3, Limit => 1000);
      assert(#topes == 244)
      elapsedTime addToCYDatabase("can-delete-me-ntfe-h11-3.dbm", topes)

      RZ = ZZ[a,b,c]
      (Qs, Xs) = readCYDatabase("can-delete-me-ntfe-h11-3.dbm", Ring => RZ);
      assert(keys Qs === toList(0..243))
      sort keys Xs
      #oo == 275

      for k in drop(sort keys Qs,10) list elapsedTime (tris = findAllFRVTs Qs#k; print (k => #tris); tris)

      findAllFRVTs Qs#101
      findAllFRVTs Qs#0
      findAllFRVTs Qs#1 -- fails
        V = normalToricVariety(rays Qs#1, first oo)
        isWellDefined V and isProjective V
      tris = findAllFRVTs Qs#2 -- needs to use ConnectedToFalse => false, crashes otherwise...!
      findAllFRVTs Qs#3
      findAllFRSTs Qs#3
      findAllFRVTs Qs#4
      findAllFRVTs Qs#5
      findAllFRVTs Qs#6
      findAllFRVTs Qs#7
      findAllFRVTs Qs#8

      elapsedTime for k in sort keys Qs list k => elapsedTime findAllFRVTs Qs#k; -- 25 sec.
      tally (oo/last)
      ooo/last/length//tally
      
///

///
  restart
  debug needsPackage "StringTorics"

  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  elapsedTime addReflexiveToCYDatabase("can-delete-me-reflexives-ntfe-h11-3.dbm", topes)

  F = openDatabase "can-delete-me-reflexives-ntfe-h11-3.dbm"

  topes = kreuzerSkarke(4, Limit => 10000);
  assert(#topes == 1197)
  elapsedTime addReflexiveToCYDatabase("can-delete-me-reflexives-ntfe-h11-4.dbm", topes)

  F = openDatabase "can-delete-me-reflexives-ntfe-h11-3.dbm"

  ///


///
-- scratch work trying to understand cohomology for bicubics.
      fcn = (b) -> (10 * binomial(b-3 + 2, 2), binomial(b+2,2))
      fcn 100
      fcn 1000
      
      T = ZZ/101[x,y,z]
      C = res ideal random(T^1, T^{10:-3})
      C.dd_2
      
      hh^*(OO_X(-3,3))
      hh^*(OO_X(-3,4)) -- 15 of these
      
      fcn = (a,b) -> (
          binomial(-a+2,2) * binomial(b-1,2),
          binomial(-a-1,2) * binomial(b+2,2))
      fcn(0,3)
      fcn(0,4)
      matrix first cohomologyMatrix(2, V, {-3,4}, first equations X)
      fcn(-3,4)
      fcn(-3,3)
      fcn(-3,2)
      for b from 3 to 20 list fcn(-3,b)
      for b from 3 to 20 list fcn(-4,b)
      for b from 3 to 20 list fcn(-5,b)
      for b from 3 to 20 list fcn(-6,b)
          binomial(-a+2,2) * binomial(b-1,2) -
          binomial(-a-1,2) * binomial(b+2,2)

      t = symbol t
      R = ZZ/101[t_0..t_9]
      M = matrix{{t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, 0, 0, 0, 0, 0},
          {0, t_0, 0, t_1, t_2, 0, t_3, t_4, t_5, 0, t_6, t_7, t_8, t_9, 0},
          {0,0, t_0, 0, t_1, t_2, 0, t_3, t_4, t_5, 0, t_6, t_7, t_8, t_9}
          }
      minimalBetti coker M
      C = res coker M
      
      S = ZZ/101[a,b,c]
      phi = map(S, R, random(S^1, S^{10:-3}))
      C = res coker phi M

///

///
   Key
     "Example: the tetraquadric Calabi-Yau 3-fold"
   Headline
     the tetraquadric
   Description
    Text
      Let's analyze one particular toric variety, and the corresponding 
      Calabi-Yau hypersurface.
      
      First, let's construct this toric variety, and the corresponding Calabi-Yau
      3-fold.  Let $V = \PP^1 \times \PP^1 \times \PP^1 \times \PP^1$, and let $X \subset V$ be defined
      by a random $(2,2,2,2)$ form in 8 variables.
    Example
      P1 = toricProjectiveSpace 1
      V0 = P1 ** P1 ** P1 ** P1
      isSmooth V0
      RZ = ZZ[a..d]
      Q = cyPolytope(rays V0)
      isFavorable Q
      hh^(1,1) Q
      hh^(1,2) Q
    Text
      $Q \subset N \otimes \RR$ is the reflexive polytope in the $N = \ZZ^4$ lattice, and 
      
      We now create the Calabi-Yau. The ring here should be in $h^(1,1)(X)$ variables (over the integers, or
      the rationals.
    Example
      X = makeCY(Q, PicardRing => RZ)
      V = normalToricVariety(X, CoefficientRing => ZZ/32003)
      dim X
      describe X
      cubicForm X
      c2Form X
    Text
      Eventually, put in the GV invariants and flop matrices one finds.
      
      
    Example
      M1 = matrix"-1,0,0,0;2,1,0,0;2,0,1,0;2,0,0,1"
      for a in {0,0,0}..{3,3,3} list a => hh^*(OO_X(-1,a#0,a#1,a#2))
      GV = gvInvariants(X, DegreeLimit => 10);
      hh^*(OO_X(-1,2,2,2))
      hh^*(OO_X(-1,2,2,3))
      hh^*(OO_X(-2,2,2,2))
      hh^*(OO_X(-2,4,4,3))
      F = first equations X
      cohomologyBasis(1, V, {-4,2,2,2})
      cohomologyBasis(1, V, {-2,4,4,4})
      cohomologyMatrix(1, V, {-2,4,4,4}, F)
      hh^*(OO_X(-4,5,0,0))
      (m, tar, src) = cohomologyMatrix(1, V, {-2,4,4,4}, F)
      matrix m;
      #tar
      #src
      rank m
///

///
-- scratch work for tetraquadric example.
    Text
      The following won't be in this tutorial.
    Example
      S = ZZ/32003[x_1,y_1,x_2,y_2,x_3,y_3, Degrees => {2:{1,0,0}, 2:{0,1,0}, 2:{0,0,1}}]
      G0 = random({2,2,2}, S)
      G1 = random({2,2,2}, S)
      G2 = random({2,2,2}, S)
      syz matrix{{G0,G1,G2}}
      degrees source oo

      A = symbol A
      C = symbol C
      S = ZZ/32003[x_1,y_1,A_1..C_3, Degrees => {2:{1,0},9:{0,1}}]
      G0 = x_1^2 * A_1 + x_1*y_1 * A_2 + y_1^2 * A_3
      G1 = x_1^2 * B_1 + x_1*y_1 * B_2 + y_1^2 * B_3
      G2 = x_1^2 * C_1 + x_1*y_1 * C_2 + y_1^2 * C_3
      syz matrix{{G0,G1,G2}}
      syz matrix{{G0,G1,G2,0},{0,G0,G1,G2}}
      M = coker matrix{{G0,G1,G2,0},{0,G0,G1,G2}}
      res M
o67_{0}
o67_{1}
      T = ZZ/32003[a,b,c]
      syz matrix{{a,b,c,0},{0,a,b,c}}
///

///
  Key
    
  Headline
  Usage
  Inputs
  Outputs
  Consequences
  Description
    Text
      This is the first way to use the package.  Grab a reflexive polytope from the Kreuzer-Skarke database.
    Example
      restart
      needsPackage "StringTorics"
      topes = kreuzerSkarke(3, Limit => 100);
      Q = cyPolytope(topes_50, ID => 50)
      RZ = ZZ[x,y,z]
      X = makeCY(Q, PicardRing => RZ, ID => 0)
      label X
      V = normalToricVariety(X, CoefficientRing => ZZ/101)
      V === ambient X
      assert(coefficientRing ring ambient X === ZZ/101)
      L = OO_X(1,2,3)
      L = OO_X(1,2,-1)
      hh^* L
      cubicForm X
      c2Form X
      intersectionNumbers X
      basisIndices X
      dim X
      aX = abstractVariety(X, base(x,y,z))
      use intersectionRing aX

      L = OO_X(1,2,-33)
      hh^* L
      chi OO(t_0  + 2*t_1 - 33*t_2)

      for d in (-2,-2,-2)..(2,2,2) list d => hh^* OO_X(d)
      chi(OO(x * t_0 + y * t_1 + z * t_2))
      RQ = ring oo
      assert(
          1/6 * sub(cubicForm X, RQ) + 1/12 * sub(c2Form X, RQ) 
          == 
          chi(OO(x * t_0 + y * t_1 + z * t_2))
          )

      -- Method #2 to use the package.      
      p1 = toricProjectiveSpace 1
      V0 = p1 ** p1 ** p1 ** p1
      isSimplicial V0
      max V0
      RZ = ZZ[a,b,c,d]
      Q1 = cyPolytope(rays V0, ID => 0)
      X = makeCY(Q1, ID => 0, PicardRing => RZ)
      V = normalToricVariety(X, CoefficientRing => ZZ/101)
      isWellDefined V
      L = OO_X(2,1,-1,-2)
      hh^* L      
      equations X

      describe X
      equations X
      X = calabiYauHypersurface V -- not written yet.
      X = toricCompleteIntersection(V, {-toricDivisor V}, Equations => generic)
        -- should detect it is CY?

      P2 = convexHull transpose matrix rays V
      isReflexive P2
      isSimplicial P2
      V1 = cyPolytope(rays V, ID => 0)
      degrees V1
      peek V1.cache
      X1 = cyData(V1, max V)
      X = variety(X1, ZZ/32003)
      hh^*(OO_X(1,1,1,1))
      cyData X

    topes = kreuzerSkarke(3, Limit => 100)      
    ks = topes_70
    V0 = normalToricVariety(ks, CoefficientRing => ZZ/101)
    Xs = makeCYs V0
    ambient Xs_0 -- a simplicial normal toric variety
    
    normalToricVariety KSEntry := opts -> ks -> (
        polytopeData := cyPolytope ks;
        normalToricVariety(rays polytopeData
        )

  V0 = normalToricVariety(ks, CoefficientRing => ZZ/101)
  Xs = findAllCYs V0 -- creates CalabiYauInToric's
  X = findOneCY V0 -- choose one triangulation
  ambient Xs_0 -- gives a simplicial toric variety (over same coefficient ring).
    
    V0 = cyPolytope ks
    for a in annotatedFaces V0 list if a#0 != dim V0 - 1 then continue else a#
    select(annotatedFaces V0, a -> a#0 == dim V0 - 1)
    PN = convexHull transpose matrix rays V0
    PN2 = polytope(ks, "N")


-- Usage #1.
  Q = cyPolytope(tope, ID => label)
  Xs = findAllCYs(Q, Ring => RZ) -- labels them
  X = makeCY(Q, ID => lab, PicardRing => RZ)
  -- given an X = Xs_0 say
  -- really want line bundles on X, cohomology on X.  But might also want equations.  Where to put those?
  
-- Usage #2. Start with a toric variety constructed elsewhere.
-- Case A: it is simplicial, from reflexive.
-- Case B: it is not simplicial, but is from reflexive.
  X = cyHypersurface(V, Ring => RZ, Equations => ...)
  V = ambient X
  
  Caveat
  SeeAlso
///



-- hodgeOfCYToricDivisor doc node moved to DocCalabiYauInToric.m2

///
  Key
    "computing sheaf cohomology of specific Calabi-Yau hypersurfaces in toric varieties"
  Headline
    facilities available for computing sheaf cohomology
  Description
    Text
    Example
  Caveat
  SeeAlso
    
///

-- "facilities available for working with triangulations" doc node moved to DocCYPolytope.m2


///
  Key
    vertexMatrix
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

-- "cohomology of line bundles on toric varieties" doc node moved to DocCI.m2

------------------------------------------------
-- Database creation and retrieval functions ---
------------------------------------------------
-- TODO: remove this node, as it has moved to DocCYPolytope.m2
///
  Key
    "Creating a CYDatabase file"
  Headline
    How to create a CYDatabase of reflexive polytopes and corresponding Calabi-Yau 3-folds
  Description
    Text
      A CYDatabase file is a @TO Database@, which contains pre-computed information about
      4-dimensional reflexive polytopes, and their triangulations.  The fact that it contains
      precomputed features allows us to spend much less time scanning over many polytopes.

      Here, we describe how to construct such a file.  We will construct the file for all examples
      with $h^{1,1} = 2$.  There are not very many of these, but the same method works for
      larger cases.
    Example
      topes = kreuzerSkarke(2, Limit => 1000);
      assert(#topes == 36)
      elapsedTime addToCYDatabase("can-delete-me-ntfe-h11-2.dbm", topes, NTFE => false)
    Text
      Let's test that this was created correctly.  We see that in particular the
      annotated faces (@TO (annotatedFaces, ReflexivePolytope)@) has been computed (this is
      one of the things that seems to take the longest.  That, and the list of triangulations.
      For higher $h^{1,1}(X)$, we must arrange to not compute these, as there are too many triangulations.
      This has not been done yet.

      The ring of a @TO CalabiYauInToric@ is the intersection ring: a ring over $\ZZ$ in
      $h^{1,1}(X)$ variables.  This is the same for all Calabi-Yaus with the same $h^{1,1}(X)$, so we
      create one such ring, and use it for all constructed Calabi-Yau's.
    Example
      RZ = ZZ[a,b]
      (Qs, Xs) = readCYDatabase("can-delete-me-ntfe-h11-2.dbm", Ring => RZ);
      assert(keys Qs === toList(0..35))
      sort keys Xs
      #oo == 36
    Text
      The line above confirms that there are 36 constructed Calabi-Yau 3-fold hypersurfaces of
      $h^{1,1} = 2$.  Note that each polytope only gives one triangulation.  As  $h^{1,1}$ increases,
      the number of triangulations increases, and at some point, becomes astronomical.

      For now, we take one, and look at some of its invariants.
    Example
      X = Xs#(0,0)
      label X
      hh^(1,1) X
      hh^(1,2) X
      cubicForm X
      c2Form X
      removeFile "can-delete-me-ntfe-h11-2.dbm"
    Text
      Now let's deal with a larger case, when we want to have multiple processes working on creating
      databases, and then we will merge them together.
    Pre
      debug needsPackage "StringTorics"
      DBNAME = "can-delete-me-cys-ntfe-h11-5.dbm"
      topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
      assert(#topes == 4990)
      "topes-h11-5.txt" << toExternalString topes << endl << close;
      topes2 = value get "topes-h11-5.txt";
      assert(topes === topes2)
      createM2Lines("cys-h11-5", "topes-h11-5.txt", 4990, 22)
      elapsedTime addToCYDatabase(DBNAME, topes) -- on Apple M4 Max, Jan 2025: 3.14 hours to create.
    Text
      Here we try to run lots of M2's to do this
    Pre
      needsPackage "StringTorics"
      R = ZZ[a..f]
      filename = elapsedTime processCYPolytopes("cys-h11-6", "topes-h11-6", 2, 3000)
      filename = elapsedTime addToCYDatabase("cys-h11-6", "topes-h11-6", 2, 3000)
      filename = elapsedTime addToCYDatabase("cys-h11-6", "topes-h11-6", 2999, 3000)
      (Qs, Xs) = readCYDatabase(filename, Ring => ZZ[a_0..a_5]);

       M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5' -e 'hi=8' -e 'processCYPolytopes("cys-h11-6", "topes-h11-6", (lo,hi))' -e 'exit 0'
       M2 --silent --stop -e 'needsPackage "StringTorics"' -e 'lohi = (5,8)' -e 'processCYPolytopes("cys-h11-6", "topes-h11-6", lohi)' -e 'exit 0'

      filename = elapsedTime processCYPolytopes("cys-h11-6", "topes-h11-6", (5,8))
      (Qs, Xs) = readCYDatabase(filename, Ring => ZZ[a_0..a_5]);
      
      topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
      assert(#topes == 4990)
      elapsedTime("topes-h11-5" << toExternalString topes << close)
      elapsedTime get "topes-h11-5";
      elapsedTime value oo;
      oo === topes

      topes = kreuzerSkarke(6, Limit => 20000); -- 17101 of these
      assert(#topes == 17101)
      elapsedTime("topes-h11-6" << toExternalString topes << close);
      topes2 = elapsedTime value get "topes-h11-6";


      
      
      DBNAME = "dbm-h11-6-range-5-20"
      elapsedTime addToCYDatabase(DBNAME, topes_{5..20})
      M2 --stop -e 'needsPackage "StringTorics" -e 'lo=5' -e 'hi=10' -e 'addToCYDatabase(DBNAME|"-range-"|lo|"-"|hi, topes

      cyDatabase(DBNAMEPREFIX, topesFile, lo, hi);
  SeeAlso
    addToCYDatabase
    readCYDatabase
///

-- addToCYDatabase doc nodes moved to DocDatabases.m2

///
  Key
  Headline
  Usage
  Inputs
  Outputs
  Consequences
    Item
  Description
    Text
    Example
  SeeAlso
///



///
-*
  restart
  needsPackage "StringTorics"
*-
  topes = kreuzerSkarke(5, Limit=>10, Access=>"wget")

  A = matrix topes_9
  P = convexHull A
  P2 = polar P
  LP = matrix{latticePoints P2}
  elems2 = regularSubdivision(LP, matrix{{1/2, 1/10, 10, 3, 2, 8, 12, 7/8, 9/11, 1}})
  hts = for i from 0 to numcols LP-1 list random 100
  elems3 = regularSubdivision(LP, matrix{hts})
  elems = regularFineTriangulation LP
  for e in elems list latticeVolume convexHull LP_e
  for e in elems2 list latticeVolume convexHull LP_e
  sum for e in elems3 list latticeVolume convexHull LP_e
  sort unique flatten elems == toList(0..numcols LP-1)
  sort unique flatten elems2 == toList(0..numcols LP-1)
  sort unique flatten elems3 == toList(0..numcols LP-1)
  -- I would like to check that these are well-defined triangulations
  latticeVolume P
  latticeVolume P2
  
  vertices P -- matrix over QQ
  vertexList P -- list (over ZZ or QQ?), different order than 'vertices'
  vertexMatrix P
  faceList
  faceDimensionHash
  minimalFace
  latticePoints P -- list of single column matrices, over ZZ
  latticePointList P
  latticePointHash
  interiorLatticePoints
  interiorLatticePointList
  faces(3,P)
  faceList(1,P)
  vertices P2
  matrix{latticePoints P2}
  
///

///
-*
  restart
  needsPackage "StringTorics"
*-
  topes = kreuzerSkarke(80, Limit=>10, Access=>"wget")
  A = matrix topes_9
  P = convexHull A
  P2 = polar P
  vertices P -- matrix over QQ
  vertexList P -- list (over ZZ or QQ?), different order than 'vertices'
  vertexMatrix P
  faceList
  faceDimensionHash
  minimalFace
  latticePoints P -- list of single column matrices, over ZZ
  latticePointList P
  latticePointHash
  interiorLatticePoints
  interiorLatticePointList
  faces(3,P)
  faceList(1,P)
  vertices P2
  matrix{latticePoints P2}
  
///

------------------------------------------------------------
-- Routines for complete intersections in toric varieties --
-- Some functionality is only for hypersurfaces! -----------
------------------------------------------------------------

-- CompleteIntersectionInToric doc node moved to DocCI.m2





-*
///
  Key
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
*-

