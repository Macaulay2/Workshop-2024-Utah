doc ///
    Key
        "creating and using Batyrev Calabi-Yau hypersurfaces"
    Headline
        information and basic constructors
    Description
        Text
          A @TO CalabiYauInToric@ is the Macaulay2 type for such an object.
          The basic data stored is the vertices of a reflexive polytope,
          and a triangulation (given as a list of lists of indices of the lattice points
          of the polytope).  Note that this is exactly the same data as that
          for a simplicial toric variety.
        Text
          It is generally easier to have triangulations computed for you.
          In many smaller cases ($h^{1,1}(X)$ less than 5 or 6 or so), the
          number of triangulations is small.  For larger cases, the number can be too
          many to enumerate...
    	Text
    	    @SUBSECTION "Constructors for Batryev Calabi Yaus"@
        Text
    	    @UL {
                TO (calabiYau, ReflexivePolytope, List),
                TO (makeCY, ReflexivePolytope),
                TO (makeCYs, ReflexivePolytope)
    	    }@
        Text
    	    @SUBSECTION "Basic information"@
        Text
    	    @UL {
                TO (dim, CalabiYauInToric),
                TO (rays, CalabiYauInToric),
                TO (max, CalabiYauInToric),
                TO (label, CalabiYauInToric),
                TO (isFavorable, CalabiYauInToric),
                TO (hh, Sequence, CalabiYauInToric),
                TO (describe, CalabiYauInToric)
    	    }@
        Text
    	    @SUBSECTION "Picard group"@
        Text
    	    @UL {
                TO (picardRing, CalabiYauInToric),
                TO (degrees, CalabiYauInToric),
                TO (basisIndices, CalabiYauInToric)
            }@
        Text
            @SUBSECTION "Interacting with NormalToricVarieties, Polyhedra"@
        Text
            @UL {
                TO (normalToricVariety, CalabiYauInToric),
                TO (ambient, CalabiYauInToric),
                TO (polytope, CalabiYauInToric),
                TO (abstractVariety, CalabiYauInToric)
            }@
///        

doc ///
    Key
        "working with complete intersections in toric varieties"
    Headline
        
    Description
        Text
          One way to construct a variety is as a complete intersection $X$ in a toric or projective
          variety $V$.  A subvariety of codimension $c$ is a complete intersection in $V$,
          whose ideal has exactly $c$ minimal generators (in the homogeneous coordinate ring, or Cox
          ring, of $V$).

          Given a set of equations cutting out $X$, there are a few caveats one must consider:
          (1) $X$ might not be smooth, and so some formulas or techniques might not work, and (2)
          the Picard group might not be the one induced by $V$.  For example, in the
          Batyrev hypersurface case, if $X$ is not favorable, then the Picard group is strictly
          larger than the induced Picard group from $V$.  We can handle that case, but for for more
          general complete intersections, it seems to be a very hard problem to determine its
          cohomology ring or Neron-Severi group.

          Still, we can learn much about $X$, especially if $X$ is smooth.
    	Text
    	    @SUBSECTION "Constructors for CompleteIntersectionInToric"@
        Text
    	    @UL {
                TO (completeIntersection, NormalToricVariety, List),
    	    }@
        Text
    	    @SUBSECTION "Basic information"@
        Text
    	    @UL {
                TO (dim, CompleteIntersectionInToric),
                TO (equations, CompleteIntersectionInToric),
    	    }@
        Text
    	    @SUBSECTION "Induced intersection ring"@
        Text
    	    @UL {
                TO (abstractVariety, CompleteIntersectionInToric),
                TO (intersectionRing, CompleteIntersectionInToric),
                TO (intersectionForm, CompleteIntersectionInToric),
                TO (c2Form, CompleteIntersectionInToric),
--                TO (cubicForm, CompleteIntersectionInToric)
            }@
        Text
    	    @SUBSECTION "Hodge numbers"@
        Text
    	    @UL {
                TO (hodgeDiamond, CompleteIntersectionInToric)
            }@
        Text
    	    @SUBSECTION "Cohomology of induced line bundles"@
        Text
    	    @UL {
                TO (symbol_, OO, CompleteIntersectionInToric),
                TO (lineBundle, CompleteIntersectionInToric, List),
                TO (cohomology, ZZ, CompleteIntersectionInToric, List, RingElement)
            }@
        Text
    	    @SUBSECTION "Hodge-Deligne polynomials, Danilov-Khovanskii algorithm"@
        Text
    	    @UL {
                "in an included package.  Needs alot of work!"
            }@
        Text
            There is alot we would like to understand still.
        Text
    	    @UL {
                "Is $X$ smooth?",
                "Many algorithms here only work in the hypersurface case"
            }@
        ///        

doc ///
  Key
    CalabiYauInToric
  Headline
    a Batyrev Calabi-Yau 3-fold hypersurface in a simplicial toric variety
  Description
    Text
      A Batyrev hypersurface in a simplicial toric 4-fold is constructed
      as follows: (1) Choose a reflexive 4D lattice polytope.
      (2) Choose a FRST triangulation of this polytope.  This is a
      triangulation of the polytope containing all lattice points
      not interior to facets, and each simplex contains the origin, and the
      triangulation is regular.
      (3) Define a complete (projective) fan $\Sigma$ by taking the cone generated
      by each simplex in the triangulation.  This defines a 4D simplicial
      toric variety $V = \PP_\Sigma$.
      (4) Choose a (smooth) anti-canonical divisor $X \subset V$ on $V$.
      Often, we take a "general" (i.e. random) such divisor.
    Text
      Batyrev proved in 1993 (XXX), that $X$ is a smooth Calabi-Yau 3-fold.
      He does this by analyzing the possible singularities, and showing that
      the singularities of $V$ must be canonical singularities.  This means points.
      $X$ is then taken to miss these finite number of points, and
      by base-point-free-ness, there can be no further singularities of $X$.
    Text
      This package defines two related types of objects: @TO ReflexivePolytope@,
      which encodes the reflexive polytope above, and also includes information
      that is contstant across triangulations (e.g. Picard group, basis for it,
      information about the polytope, etc).  The second is {\tt CalabiYauInToric},
      which includes the triangulation, equations, computation of cohomology,
      intersection rings, etc.
    Text
      As an example, let's consider the following polytope, and then its
      corresponding Calabi-Yau varieties.
    Example
      topes = kreuzerSkarke(3, Limit => 10);
      tope = topes_7;
      Q = reflexivePolytope(tope, ID => 7)
      transpose matrix vertices Q
      transpose matrix latticePoints Q
      faceDimensions Q
      X = makeCY(Q, PicardRing => ZZ[a,b,c], ID => 0)
      label X
    Text
      Here are some 
    Example
      V = normalToricVariety(X, CoefficientRing => ZZ/32003)
      ring V -- wrong ring... if findAllCYs is used...
      glsm' = transpose matrix degrees ring V
      glsm = transpose matrix degrees X
      glsm == glsm'
      isSimplicial V
      isSmooth V
      isProjective V
      dim V
      singularCones(0, V) -- V has 4 singular points
    Text
      
    Example
      --Xci = asCompleteIntersection X -- might be good?
      basisIndices X == {0,1,2}
      Xci = completeIntersection(V, {-toricDivisor V}, Basis => basisIndices X, Variables => {a,b,c})
      Xa = abstractVariety Xci
      use intersectionRing Xa
      use coefficientRing oo
      sub(integral((a*t_0 + b*t_1 + c*t_2)^3), QQ ** picardRing X)
      basisIndices X == {0,1,2}
      integral((a*t_0 + b*t_1 + c*t_2)^3)
      
      isFavorable X
      hh^(1,1) X
      hh^(1,2) X
      c2Form X
      cubicForm X
      -- hodgeDiamond X -- BUG: not available!
    Text
      {\bf Related rings and varieties}
      ideal prune intersectionRing Xa
    Text
    Example
      Xa = abstractVariety(X, base(a,b,c))
      ideal prune intersectionRing Xa
    Example
     dump X -- FAILS if ID isn't given
  SeeAlso
    ReflexivePolytope
///

doc ///
  Key
    picardRing
    (picardRing, CalabiYauInToric)
  Headline
    picard ring over the integers 
  Usage
    RZ = picardRing X
  Inputs
    X:CalabiYauInToric
  Outputs
    RZ:Ring
      This retrieves the stored Picard Ring, which is
      a polynomial ring with $h^{1,1}(X)$ variables
  Description
    Text
      This ring is used to create the $c_2(X)$ and intersection form for $X$,
      as well as to create the base rings the intersection ring
    Example
      Q = reflexivePolytope(
          {
              {-1, -1, -1, -1},
              {-1, -1, -1, 0},
              {-1, -1, 0, 2},
              {-1, 0, -1, -1},
              {0, -1, -1, -1},
              {1, -1, 0, -1},
              {1, 2, 2, 2}
              },
          ID => 7)
      isReflexive polytope Q
      R = ZZ[a,b,c]
      X = makeCY(Q, PicardRing => R, ID => 0)
      label X
      picardRing X === R
      c2 X
      c2Form X
      ring c2Form X === R
    Text
      The Picard Ring is also used when creating the abstract variety
      corresponding to $X$.
  Caveat
    This is not well thought out in the case when the class group of the ambient
    toric variety has torsion.
  SeeAlso
    c2Form
    cubicForm
    intersectionForm
    intersectionNumbers
    c2
    intersectionRing
///

doc ///
  Key
    (rays, CalabiYauInToric)
  Headline
    coordinates of the rays of the underlying simplicial fan
  Usage
    rays X
    transpose matrix rays X
  Inputs
    X:CalabiYauInToric
  Outputs
    :List
      the first lattice point on each ray of the underlying fan
      of $X$
  Description
    Text
      The fan defining the normal toric variety that $X$ sits in
      is determined by its rays and maximal cones.
    Example
      Q = reflexivePolytope(
          {
              {-1, -1, -1, -1},
              {-1, -1, -1, 0},
              {-1, -1, 0, 2},
              {-1, 0, -1, -1},
              {0, -1, -1, -1},
              {1, -1, 0, -1},
              {1, 2, 2, 2}
              },
          ID => 7)
      isReflexive polytope Q
      Ts = findAllFRSTs Q
      R = ZZ[a,b,c]
      X = calabiYau(Q, Ts#0, PicardRing => R, ID => 0)
      rays X
    Text
      It is often useful to create a matrix whose columns are these rays.
    Example
      transpose matrix rays X
  SeeAlso
    (rays, CalabiYauInToric)
    (normalToricVariety, CalabiYauInToric)
///

doc ///
  Key
    (label, CalabiYauInToric)
    (label, ReflexivePolytope)
  Headline
    label given during construction
  Usage
    label Q
    label X
  Inputs
    X:{CalabiYauInToric, ReflexivePolytope}
  Outputs
    :{ZZ, Sequence}
      the corresponding integer label of a polytope or a pair of
      integers for a Calabi Yau variety
  Description
    Text
      In order to process many @TO ReflexivePolytope@s or @TO CalabiYauInToric@s
      it is useful to label them.  It is also useful to have a labelling that
      is consistent with the Kreuzer-Skarke database in order to better inter-operate
      with other systems, such as {\tt CYTools}.
    Text
      By convention, polytopes have integers labels
      (often its index in the Kreuzer-Skarke database), and
      Calabi-Yau's have label $(polytope label, triangulation index)$.
    Example
      Q = reflexivePolytope(
          {
              {-1, -1, -1, -1},
              {-1, -1, -1, 0},
              {-1, -1, 0, 2},
              {-1, 0, -1, -1},
              {0, -1, -1, -1},
              {1, -1, 0, -1},
              {1, 2, 2, 2}
              },
          ID => 7)
      Ts = findAllFRSTs Q
      R = ZZ[a,b,c]
      X = calabiYau(Q, Ts#0, PicardRing => R, ID => 0)
      label Q -- this is the polytope at index 7 in the Kreuzer-Skarke list for h11=3
      label X
    Text
      These are often used as keys into a hash table of a bunch of polytopes
      or CY3 varieties.  For instance, can get all Batyrev Calabi-Yau hypersurfaces
      at h11=3 from the data base included with StringTorics.  The two hash tables
      returned have keys which are the labels of the resulting polytopes and
      CalabiYauInToric's.
    Example
      DB3 = databaseLOC | "/cy3-h11-3.dbm"
      (Qs, Xs) = readCYDatabase(DB3, Ring => R);
      keys Qs
      keys Xs
  SeeAlso
    reflexivePolytope
    (findAllFRSTs, ReflexivePolytope)
    calabiYau
    readCYDatabase
///

doc ///
  Key
    (max, CalabiYauInToric)
  Headline
    (indices of) the maximal cones of the underlying fan
  Usage
    tri = max X
  Inputs
    X:CalabiYauInToric
  Outputs
    tri:List
      of lists of integers, the indices into the rays of $X$
      describing the maximal cones of the corresponding
      simplicial fan
  Description
    Text
      The fan defining the normal toric variety that $X$ sits in
      is determined by its rays and maximal cones.
    Example
      Q = reflexivePolytope(
          {
              {-1, -1, -1, -1},
              {-1, -1, -1, 0},
              {-1, -1, 0, 2},
              {-1, 0, -1, -1},
              {0, -1, -1, -1},
              {1, -1, 0, -1},
              {1, 2, 2, 2}
              },
          ID => 7)
      -- TODO: isWellDefined Q
      -- TODO: isWellDefined X
      isReflexive polytope Q
      Ts = findAllFRSTs Q
      R = ZZ[a,b,c]
      X = calabiYau(Q, Ts#0, Ring => R, ID => 0)
      transpose matrix rays X
      max X
      # max X
    Text
      So there are 13 maximal cones, including
      the cone generated by the first 4 columns of the
      rays matrix.
  SeeAlso
    (rays, CalabiYauInToric)
    (normalToricVariety, CalabiYauInToric)
///

doc ///
   Key
     hodgeOfCYToricDivisor
     (hodgeOfCYToricDivisor,Polyhedron,List)
   Headline
     compute the cohomology vector of an (irreducible) toric divisor on a CY hypersurface
   Usage
     hodgeOfCYToricDivisor(P,pt)
   Inputs
     P:Polyhedron
       Any polytope will do, although so far it has only been tested on
       reflexive polytopes.
     pt:List
       a lattice point in the polar dual polytope {\tt polar P}
   Outputs
     :List
       The list of $\{ h^0(X,OO_D), h^1(X,OO_D), h^2(X,OO_D) \}$,
       where $D$ is the intersection of the toric divisor corresponding
       to the lattice point with a hypersurface $X$ corresponding to
       an anti-canonical divisor on the toric 4-fold $V$ corresoponding .
   Description
    Text
      We assume that $P$ is a reflexive 4-dimensional polytope, and let $V_0$ be the 4-dimensional
      toric variety corresponding to $P$ (i.e. corresponding to the normal fan of $P$).  Let $V$
      be a simplicial resolution of $V_0$, corresponding to a star triangulation of the polar dual
      $P^o$, which is fine, i.e. involves all of the lattice points of $P^o$.  Let $X \subset V$ be
      the inverse image of an anti-canonical divisor on $V_0$.  Note that from Batyrev, it turns
      out that $X$ is a smooth Calabi-Yau 3-fold.  Finally, a lattice point $pt$
      of $P^o$ corresponds to a toric divisor on $V$.  Let $D$ be the intersection of this divisor
      with $X$.

      This function computes the vector of cohomologies of the structure sheaf of the surface $D$.

      As an example, the following example is taken from the Kreuzer-Skarke database of 4D reflexive
      polytopes.

    Example
       polystr = "Kreuzer-Skarke: 4 12  M:24 12 N:16 11 H:11,19 [-16]
               1   0   0   0   0   1   2   1   0  -2   0  -2
               0   1   0   0   0   0  -2  -1   1   2  -1   0
               0   0   1   0   0  -1   0  -1  -1   1  -1   1
               0   0   0   1  -1   0   1   1  -1   0   1  -2
               "
      A = matrix first kreuzerSkarke polystr
      P = convexHull A
    Text
      This polytope has 12 vertices, 33 edges, 32 2-faces, and 11 facets.
    Example
      # faceList(0,P)
      # faceList(1,P)
      # faceList(2,P)
      # faceList(3,P)
      fVector P
    Text

      Note that only the faces which contain interior lattice points, or whose dual does, is included.
      So 6 of the 33 edges of the polytope have an interior vertex along that edge.

      There are 15 non-zero lattice points in the dual, meaning that the
      toric 4-fold has 15 toric divisors on it (each is a 3-fold, and also toric).
      It turn out that they all are rigid, in the sense that in each case,
      $h^0(OO_D) = 1$, $h^1(OO_D) = 0$, $h^2(OO_D) = 0$.
    Example
      # latticePointList polar P
      hodgeOfCYToricDivisors P
    Text
      The Hodge numbers $h^{1,1}(X)$ and $h^{2,1}(X)$ can be computed using
      information about $P$ only, not the specific triangulation used.
    Example
      isFavorable P
      h11OfCY P
      h21OfCY P
   Caveat
     This function currently only works for 4-d reflexive polytopes.  However, the
     formulas work for other dimensions, and these should be included.
   SeeAlso
     (hh, Sequence, CalabiYauInToric)
///

doc ///
  Key
    calabiYau
    (calabiYau, ReflexivePolytope, List)
  Headline
    construct a Calabi-Yau hypersurface from a reflexive polytope and triangulation
  Usage
    X = calabiYau(Q, T)
  Inputs
    Q:ReflexivePolytope
    T:List
      a triangulation of $Q$, given as a list of lists of indices into @TT "rays Q"@
    ID => ZZ
      an optional integer label for this CY
    PicardRing => Ring
      a polynomial ring with $h^{1,1}$ variables, used for intersection forms
  Outputs
    X:CalabiYauInToric
  Description
    Text
      Constructs a @TO CalabiYauInToric@ from a reflexive polytope $Q$ and a
      fine, regular, star triangulation $T$.  The triangulation is given as a list
      of maximal simplices, each a list of indices into @TT "rays Q"@.

      See @TO (makeCY, ReflexivePolytope)@ for a version that automatically computes
      a triangulation, and @TO (findAllCYs, ReflexivePolytope)@ for finding all CYs
      up to equivalence.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      Ts = findAllFRSTs Q
      R = ZZ[a,b,c]
      X = calabiYau(Q, Ts#0, PicardRing => R, ID => 0)
      label X
      hh^(1,1) X
      hh^(1,2) X
  SeeAlso
    (makeCY, ReflexivePolytope)
    (findAllCYs, ReflexivePolytope)
    CalabiYauInToric
///

doc ///
  Key
    makeCY
    (makeCY, ReflexivePolytope)
  Headline
    construct a Calabi-Yau using an automatically computed triangulation
  Usage
    X = makeCY Q
  Inputs
    Q:ReflexivePolytope
    ID => ZZ
      an optional integer label
    PicardRing => Ring
      a polynomial ring with $h^{1,1}$ variables
  Outputs
    X:CalabiYauInToric
  Description
    Text
      Constructs a @TO CalabiYauInToric@ from $Q$ by automatically computing
      one fine, regular, star triangulation.  This is a convenient shortcut when
      you just need one CY from a given polytope.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c]
      X = makeCY(Q, PicardRing => R, ID => 0)
      max X
      hh^(1,1) X
      cubicForm X
  SeeAlso
    calabiYau
    (findAllCYs, ReflexivePolytope)
///

doc ///
  Key
    makeCYs
    (makeCYs, ReflexivePolytope)
    findAllCYs
    (findAllCYs, ReflexivePolytope)
  Headline
    find all Calabi-Yau hypersurfaces up to equivalence
  Usage
    Xs = findAllCYs Q
    Xs = makeCYs Q
  Inputs
    Q:ReflexivePolytope
    PicardRing => Ring
      a polynomial ring with $h^{1,1}$ variables
    NTFE => Boolean
      whether to partition by NTFE (non-trivial face equivalence). Default: @TT "true"@
    Automorphisms => Boolean
      whether to use automorphisms to reduce the list. Default: @TT "true"@
  Outputs
    Xs:List
      of @TO CalabiYauInToric@ objects, one per equivalence class
  Description
    Text
      Computes all fine, regular, star triangulations of $Q$, partitions them
      by 2-face equivalence (modulo automorphisms of $Q$ if requested), and
      returns one @TO CalabiYauInToric@ per equivalence class.

      @TT "makeCYs"@ and @TT "findAllCYs"@ are synonyms.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c]
      Xs = findAllCYs(Q, PicardRing => R)
      #Xs
      Xs/cubicForm
  SeeAlso
    calabiYau
    (makeCY, ReflexivePolytope)
    partitionFRSTsByDFaceEquivalence
///

doc ///
  Key
    computeBasics
    (computeBasics, ReflexivePolytope)
  Headline
    compute and cache basic data for a reflexive polytope
  Usage
    computeBasics Q
  Inputs
    Q:ReflexivePolytope
  Description
    Text
      Computes and caches all basic data for $Q$: lattice points, face dimensions,
      annotated faces, GLSM charges, Hodge numbers ($h^{1,1}$, $h^{1,2}$),
      automorphisms, and all fine regular star triangulations.

      This is a convenience function that calls all the individual computation
      functions at once.  After calling @TT "computeBasics"@, subsequent calls
      to functions like @TT "hh^(1,1) Q"@, @TT "basisIndices Q"@, etc., will
      use the cached values.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      computeBasics Q
      peek Q.cache
      hh^(1,1) Q
      hh^(1,2) Q
      basisIndices Q
      # findAllFRSTs Q
      # automorphisms Q
  SeeAlso
    reflexivePolytope
    basisIndices
    (findAllFRSTs, ReflexivePolytope)
///

doc ///
  Key
    basisIndices
    (basisIndices, ReflexivePolytope)
    (basisIndices, CalabiYauInToric)
  Headline
    indices of the rays forming a basis for the Picard group
  Usage
    L = basisIndices Q
    L = basisIndices X
  Inputs
    Q:ReflexivePolytope
    X:CalabiYauInToric
  Outputs
    L:List
      of integers (for favorable polytopes) or of integers and pairs (for non-favorable)
  Description
    Text
      Returns the indices of the rays of $Q$ (or the underlying polytope of $X$)
      that form a basis for the Picard group of the corresponding Calabi-Yau.
      For favorable polytopes, these are simply integer indices into @TT "rays Q"@.
      For non-favorable polytopes, some entries may be pairs $(i, j)$ indicating
      a divisor coming from the interior of a 2-face.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      basisIndices Q
      isFavorable Q
    Example
      X = makeCY(Q, PicardRing => ZZ[a,b,c], ID => 0)
      basisIndices X
  SeeAlso
    picardRing
    isFavorable
    (degrees, ReflexivePolytope)
///

doc ///
  Key
    isFavorable
    (isFavorable, ReflexivePolytope)
    (isFavorable, CalabiYauInToric)
  Headline
    whether the polytope or CY is favorable
  Usage
    isFavorable Q
    isFavorable X
  Inputs
    Q:ReflexivePolytope
    X:CalabiYauInToric
  Outputs
    :Boolean
  Description
    Text
      A reflexive polytope is {\em favorable} if the Picard group of the corresponding
      Calabi-Yau hypersurface is generated by (restrictions of) toric divisors corresponding to
      rays on the 1-skeleton of the polytope.  Equivalently, no 2-face
      of the dual polytope has interior lattice points that contribute to $h^{1,1}$.

      Most polytopes in the Kreuzer-Skarke database are favorable, especially at
      small $h^{1,1}$.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      isFavorable Q
    Example
      X = makeCY(Q, PicardRing => ZZ[a,b,c], ID => 0)
      isFavorable X
  SeeAlso
    basisIndices
    annotatedFaces
///

doc ///
  Key
    (hh, Sequence, CalabiYauInToric)
  Headline
    Hodge numbers of a Calabi-Yau hypersurface
  Usage
    hh^(p,q) X
  Inputs
    (p,q):Sequence
      a pair of non-negative integers
    X:CalabiYauInToric
  Outputs
    :ZZ
  Description
    Text
      Returns the Hodge number $h^{p,q}(X)$.  Currently only implemented for
      Calabi-Yau 3-folds (from 4D reflexive polytopes), so $p,q \in \{0,1,2,3\}$.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      X = makeCY(Q, PicardRing => ZZ[a,b,c], ID => 0)
      hh^(1,1) X
      hh^(1,2) X
      hh^(0,0) X
  SeeAlso
    (hh, Sequence, ReflexivePolytope)
///

doc ///
  Key
    (hh, Sequence, ReflexivePolytope)
  Headline
    Hodge numbers of the corresponding Calabi-Yau 3-fold
  Usage
    hh^(p,q) Q
  Inputs
    (p,q):Sequence
      a pair of non-negative integers
    Q:ReflexivePolytope
  Outputs
    :ZZ
  Description
    Text
      Returns the Hodge number $h^{p,q}$ of the Calabi-Yau 3-fold corresponding
      to the reflexive polytope $Q$ (via the Batyrev construction).  These are
      independent of the choice of triangulation.

      Currently only implemented for 4D reflexive polytopes.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      hh^(1,1) Q
      hh^(1,2) Q
  SeeAlso
    (hh, Sequence, CalabiYauInToric)
    isFavorable
///

doc ///
  Key
    (dim, CalabiYauInToric)
  Headline
    dimension of a Calabi-Yau hypersurface
  Usage
    dim X
  Inputs
    X:CalabiYauInToric
  Outputs
    :ZZ
  Description
    Text
      Returns the dimension of $X$, which is one less than the dimension of
      the ambient toric variety.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      X = makeCY(Q, PicardRing => ZZ[a,b,c], ID => 0)
      dim X
  SeeAlso
    CalabiYauInToric
///

doc ///
  Key
    (reflexivePolytope, CalabiYauInToric)
    (cyPolytope, CalabiYauInToric)
  Headline
    the underlying reflexive polytope of a Calabi-Yau
  Usage
    Q = reflexivePolytope X
    Q = cyPolytope X
  Inputs
    X:CalabiYauInToric
  Outputs
    Q:ReflexivePolytope
  Description
    Text
      Returns the @TO ReflexivePolytope@ from which $X$ was constructed.
      @TT "cyPolytope"@ is a synonym.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      X = makeCY(Q, PicardRing => ZZ[a,b,c], ID => 0)
      reflexivePolytope X === Q
      cyPolytope X === Q
  SeeAlso
    CalabiYauInToric
    ReflexivePolytope
///

doc ///
  Key
    h11OfCY
    (h11OfCY, Polyhedron)
  Headline
    compute h11 of the CY 3-fold from a 4D reflexive polytope
  Usage
    n = h11OfCY P
  Inputs
    P:Polyhedron
      a 4-dimensional reflexive polytope
  Outputs
    n:ZZ
  Description
    Text
      Computes $h^{1,1}$ of the Calabi-Yau 3-fold corresponding to $P$
      using the Batyrev formula.  This works directly on a @TO Polyhedron@ object.
      For the @TO ReflexivePolytope@ version, use @TT "hh^(1,1) Q"@.
    Example
      P = convexHull transpose matrix {
          {-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
          {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}}
      h11OfCY P
  SeeAlso
    h21OfCY
    (hh, Sequence, ReflexivePolytope)
///

doc ///
  Key
    h21OfCY
    (h21OfCY, Polyhedron)
  Headline
    compute h21 of the CY 3-fold from a 4D reflexive polytope
  Usage
    n = h21OfCY P
  Inputs
    P:Polyhedron
      a 4-dimensional reflexive polytope
  Outputs
    n:ZZ
  Description
    Text
      Computes $h^{2,1}$ of the Calabi-Yau 3-fold corresponding to $P$
      using the Batyrev formula.  This works directly on a @TO Polyhedron@ object.
      For the @TO ReflexivePolytope@ version, use @TT "hh^(1,2) Q"@.
    Example
      P = convexHull transpose matrix {
          {-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
          {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}}
      h21OfCY P
  SeeAlso
    h11OfCY
    (hh, Sequence, ReflexivePolytope)
///

doc ///
  Key
    label
  Headline
    label given during construction
  Description
    Text
      Returns the label (ID) associated with a @TO ReflexivePolytope@ or
      @TO CalabiYauInToric@.  See @TO (label, CalabiYauInToric)@ and
      @TO (label, ReflexivePolytope)@ for details.
  SeeAlso
    (label, CalabiYauInToric)
    (label, ReflexivePolytope)
///

