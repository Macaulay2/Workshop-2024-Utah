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
      by base-point-free-ness, there can be no further sinularities of $X$.
    Text
      This package defines two related types of objects: @TO ReflexivePolytope@,
      which encodes the reflexive polytpoe above, and also includes information
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
    toric variety has torsion.  It might have some issues for
    the 
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
      DB3 = StringTorics#"auxiliary files" | "cy3-h11-3.dbm"
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

