---------------------------------------------------------------------------
-- Doc for complete intersections and hypersurfaces in a toric variety ----
---------------------------------------------------------------------------

doc ///
   Key
     "Example: (3,3) hypersurface in P2 x P2"
   Headline
     the bicubic threefold
   Description
    Text
      Let's analyze one particular toric variety, and the corresponding
      Calabi-Yau hypersurface.

      First, let's construct this toric variety, and the corresponding Calabi-Yau
      3-fold.  Let $V = \PP^2 \times \PP^2$, and let $X \subset V$ be defined
      by a random $(3,3)$ form in 6 variables.
    Example
      P2 = toricProjectiveSpace 2
      V0 = P2 ** P2
      isSmooth V0
      RZ = QQ[a,b]
      Q = reflexivePolytope(rays V0)
      isFavorable Q
      hh^(1,1) Q
      hh^(1,2) Q
    Text
      $Q \subset N \otimes \RR$ is the reflexive polytope in the $N = \ZZ^3$ lattice, and

      We now create the Calabi-Yau. The ring here should be in $h^(1,1)(X)$ variables (over the integers, or
      the rationals.
    Example
      X = makeCY(Q, PicardRing => RZ)
      normalToricVariety(X, CoefficientRing => ZZ/32003)
      dim X
      describe X
      cubicForm X
      c2Form X
    Text
      Hirzebruch-Riemann-Roch gives the following for the euler characteristic of OO(a,b).
    Example
      aX = abstractVariety(X, base(a,b))
      basisIndices X
      intersectionRing aX
      chi OO(a * t_0 + b*t_1)
      1/6 * cubicForm X + 1/12 * c2Form X
    Text
      Now let's investigate the cohomology of line bundles $\mathcal{O}_X(a,b)$.
      Very few cohomologies are non-zero on the ambient $V = \PP^2 \times \PP^2$.
    Example
      V = ambient X
      netList toricOrthants V
    Text
      For $a,b \ge 0$, $H^0(\mathcal{O}_X(a,b)) = H^0(\mathcal{O}_V(a,b)) = S_{ab}$, where $S$ is the Cox ring of $V$,
      and is zero outside of this range.
    Example
      S = ring V
      describe S
      cohomologyBasis(2, V, {-4,0})
      (lo,hi) = (-4, 4)
      elapsedTime cohoms = hashTable flatten for a from lo to hi list for b from lo to hi list elapsedTime (a,b) => (hh^*(OO_X(a,b)))
      matrix for a from lo to hi list for b from lo to hi list (cohoms#(a,b))_0
      matrix for a from lo to hi list for b from lo to hi list (cohoms#(a,b))_1
      matrix for a from lo to hi list for b from lo to hi list (cohoms#(a,b))_2
      matrix for a from lo to hi list for b from lo to hi list (cohoms#(a,b))_3
      for a from 0 to 10 list a => hh^*(OO_X(0,a))

      cohomologyBasis(2, V, {-6,3})
      cohomologyBasis(2, V, {-3,6})
      cohomologyMatrix(2, V, {-3, 6}, first equations X)
      matrix first oo
      rank oo
      for a from -3 to 3 list cohomologyBasis(2, V, {a,0})
      for a from -3 to 3 list cohomologyBasis(2, V, {-3,a-3})

///

doc ///
  Key
    "cohomology of line bundles on toric varieties"
  Headline
    introduction to computing line bundle cohomologies in Macaulay2
  Description
    Text
    Example
      topes = kreuzerSkarke(3, Access=>"wget");
      topes_9
      A = matrix topes_11
      P = convexHull A
      P2 = polar P
      V = reflexiveToSimplicialToricVariety P
      S = ring V -- the Cox ring
      GLSM = transpose matrix degrees S
      SR = dual monomialIdeal V
      assert isSimplicial V
      assert isSmooth V
      picardGroup V === ZZ^3
    Text
      As an example, let's compute the cohomology of the line bundle $OO_V(-2,3,-4)$.
    Example
      D = 2*V_1 + 3*V_0 - 4*V_6
      degree D
      for i from 0 to 4 list HH^i(V, OO(D))
      for i from 0 to 4 list rank HH^i(V, OO(D))
      cohomologyVector(V, {-2,3,-4})
      cohomologyVector(V, D)
      hashTable for i from 0 to # rays V - 1 list i => cohomologyVector(V, V_i)
    Text
      @SUBSECTION "Bases of cohomology groups"@
    Text
      We now delve a bit deeper into the bases of these cohomology groups.

      There are a number of pointed orthants in $\ZZ^n$, where $n$ is the number of rays,
      which support cohomology.
    Example
      netList toricOrthants V
      cohomologyBasis(1, V, {2, 1, 0})
      HH^1(V, OO_V(2,1,0))
      for i from 0 to 3 list cohomologyBasis(i, V, {2, 1, 0})
      for i from 0 to 3 list rank HH^i(V, OO_V(2, 1, 0))

      for i from 0 to 3 list cohomologyBasis(i, V, {2, 1, -2})
      for i from 0 to 3 list cohomologyBasis(i, V, {2, 1, -3})
    Text
      @SUBSECTION "Cohomology on Calabi-Yau hypersurfaces"@
    Text
  SeeAlso
///

doc ///
  Key
      CompleteIntersectionInToric
  Headline
      a complete intersection in a projective toric variety
  Description
    Text
      Functions for this type include allowing intersection theory on the induced
      intersection ring
    Text
      @SUBSECTION "Line bundles and cohomology on complete intersections in toric varieties"@
    Text
      @UL {
          -- {TO ""},
          -- {TO ""},
          -- {TO ""}
          }@
    Text
      Here is an example of using these facilities.  We consider a hypersurface in a normal toric
      variety.
    Example
      needsPackage "StringTorics"
      V = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/32003)
      rays V
      max V
      dual monomialIdeal V
      X = completeIntersection(V, {-toricDivisor V})
      dim X == 2
      -- TODO: would like to check smoothness (one way: saturate(ideal minors(1, jacobian ideal equations X), ideal V))
      -- TODO: codim
    Example
      pt = base(a,b,c)
      Xa = abstractVariety(X, pt)
      IX = intersectionRing Xa
      numgens IX
      L = OO_X(1,0,0)
      hh^* L
      hh^1 L
      hh^10 L
  Caveat
  SeeAlso
///
