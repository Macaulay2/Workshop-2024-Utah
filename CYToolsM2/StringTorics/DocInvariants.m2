------------------------------------------------------------
-- Invariants of cubic forms and c2 forms                 --
------------------------------------------------------------

doc ///
  Key
    "invariants of cubic forms and c2 forms"
  Headline
    topological invariants for distinguishing Calabi-Yau 3-folds
  Description
    Text
      This section describes various numerical invariants that can be computed
      from the cubic intersection form $F$ and the $c_2$ form $L$ of a
      Calabi-Yau 3-fold.  These invariants are preserved under integer
      changes of basis and so can be used to distinguish topologically
      inequivalent CYs.

      The invariants range from simple (GCD of coefficients) to more
      sophisticated (Aronhold invariants of plane cubics, Hessian factorization
      shapes, conductor invariants of singular loci).
    Text
      @SUBSECTION "Hubsch invariants"@
    Text
      @UL {
          TO (hubschInvariants, CalabiYauInToric)
      }@
    Text
      @SUBSECTION "Content and Hodge invariants"@
    Text
      @UL {
          TO (invariantsH11H12, CalabiYauInToric),
          "TO (invariantContents, CalabiYauInToric)"
      }@
    Text
      @SUBSECTION "Hessian invariants"@
    Text
      @UL {
          TO (hessianInvariants, CalabiYauInToric)
      }@
    Text
      @SUBSECTION "Singular locus invariants"@
    Text
      @UL {
          TO (singularContents, CalabiYauInToric),
          TO (cubicConductorInvariants, CalabiYauInToric),
          TO (cubicLinearConductorInvariants, CalabiYauInToric),
          TO (singularContentsQuartic, CalabiYauInToric)
      }@
    Text
      @SUBSECTION "Aronhold invariants"@
    Text
      @UL {
          TO (aronhold, RingElement)
      }@
    Text
      @SUBSECTION "Point counting"@
    Text
      @UL {
          TO PointCounter,
          TO (pointCounter, Ring),
          TO (pointCounts, PointCounter, RingElement, RingElement)
      }@
  SeeAlso
    "topological equivalence of Calabi-Yau 3-folds"
    "working with intersection numbers and intersection rings"
///

----------------------------------------------
-- Hubsch invariants -------------------------
----------------------------------------------

doc ///
  Key
    hubschInvariants
    (hubschInvariants, CalabiYauInToric)
  Headline
    compute Hubsch divisibility invariants of a Calabi-Yau 3-fold
  Usage
    L = hubschInvariants X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:Sequence
      $(d_1, d_2, d_3, d_4, d_5, d_6, d_7, g)$ where
      $d_1, \ldots, d_7$ are divisibility invariants and $g = \gcd(c_2)$
  Description
    Text
      Computes the Hubsch divisibility invariants of the intersection numbers
      and $c_2$ form of a Calabi-Yau 3-fold.  These are GCD-based invariants
      that are preserved under integer changes of basis.

      The first three invariants ($d_1, d_2, d_3$) come from the cubic
      intersection form alone:
      $d_1$ is the GCD of all intersection numbers,
      $d_2$ involves pairs of equal indices,
      $d_3$ involves triple-equal indices.
      The next four ($d_4, d_5, d_6, d_7$) come from a quadrilinear form
      derived from both the cubic form and $c_2$.
      The last value $g$ is the GCD of the $c_2$ form.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      hubschInvariants X
  SeeAlso
    invariants
    intersectionNumbers
    c2Form
///

----------------------------------------------
-- Content and Hodge invariants --------------
----------------------------------------------

doc ///
  Key
    invariantsH11H12
    (invariantsH11H12, CalabiYauInToric)
  Headline
    Hodge numbers as an invariant list
  Usage
    L = invariantsH11H12 X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:List
      $\{h^{1,1}, h^{1,2}\}$
  Description
    Text
      Returns the Hodge numbers $h^{1,1}$ and $h^{1,2}$ as a two-element list.
      This is a simple invariant used as a first filter when partitioning
      CYs by topology.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      invariantsH11H12 X
  SeeAlso
    (hh, Sequence, CalabiYauInToric)
///

doc ///
  Key
    invariantContents
    (invariantContents, CalabiYauInToric)
  Headline
    GCD of coefficients of the c2 and cubic forms
  Usage
    L = invariantContents X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:List
      $\{$content of $c_2$, content of cubic form$\}$
  Description
    Text
      Returns a list of two integers: the GCD of the coefficients
      (``polynomial content'') of the $c_2$ form and of the cubic form.
      These are simple invariants preserved under integer changes of basis.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      invariantContents X
  SeeAlso
    c2Form
    cubicForm
///

----------------------------------------------
-- Hessian invariants ------------------------
----------------------------------------------

doc ///
  Key
    hessianInvariants
    (hessianInvariants, CalabiYauInToric)
  Headline
    factorization shape of the Hessian determinant
  Usage
    S = hessianInvariants X
  Inputs
    X:CalabiYauInToric
  Outputs
    S:List
      the factorization shape of $\det(\mathrm{Hess}(F))$
  Description
    Text
      Computes the Hessian matrix of the cubic form $F$ (the matrix of
      second partial derivatives), takes its determinant, and returns the
      ``factorization shape'': a sorted list of pairs $\{$degree, multiplicity$\}$
      describing the irreducible factorization of the determinant over $\mathbb{Z}$.

      This is an invariant of the cubic form preserved under $GL(n, \mathbb{Z})$
      changes of basis.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      hessianInvariants X
  SeeAlso
    cubicForm
    hubschInvariants
///

----------------------------------------------
-- Singular locus invariants -----------------
----------------------------------------------

doc ///
  Key
    singularContents
    (singularContents, CalabiYauInToric)
    (singularContents, Ideal)
  Headline
    GCD structure of the singular locus of the cubic form
  Usage
    L = singularContents X
    L = singularContents I
  Inputs
    X:CalabiYauInToric
    I:Ideal
      a homogeneous ideal over $\mathbb{Z}$
  Outputs
    L:List
      cumulative GCDs by degree
  Description
    Text
      Computes the singular locus of the cubic form $F$ (the saturation of
      $\langle F \rangle + \mathrm{Jac}(F)$) and returns the cumulative GCDs
      of the generators organized by degree.

      For each degree $d$, the result gives $\gcd$ of the polynomial contents
      of all generators of degree $\leq d$.  This list is an invariant of the
      cubic form.

      The ideal version operates on any homogeneous ideal over $\mathbb{Z}$.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      singularContents X
  SeeAlso
    cubicConductorInvariants
    cubicForm
///

doc ///
  Key
    cubicConductorInvariants
    (cubicConductorInvariants, CalabiYauInToric)
  Headline
    conductor-type invariants of the singular locus of the cubic form
  Usage
    L = cubicConductorInvariants X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:List
      $\{$integer part, linear content$\}$ of the saturated Jacobian ideal
  Description
    Text
      Computes $\langle F \rangle + \mathrm{Jac}(F)$, saturates it,
      and returns two values: the integer part (the GCD of the constant
      generators) and the linear content (the GCD of the linear generators).

      These capture information about the singularities of the cubic surface
      defined by $F$ in projective space.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      cubicConductorInvariants X
  SeeAlso
    singularContents
    cubicLinearConductorInvariants
///

doc ///
  Key
    cubicLinearConductorInvariants
    (cubicLinearConductorInvariants, CalabiYauInToric)
  Headline
    conductor invariant of the intersection of the cubic and c2 forms
  Usage
    n = cubicLinearConductorInvariants X
  Inputs
    X:CalabiYauInToric
  Outputs
    n:ZZ
      the integer part of the saturated Jacobian ideal of $\langle L, F \rangle$
  Description
    Text
      Computes the ideal $\langle L, F \rangle$ (where $L$ is the $c_2$ form
      and $F$ the cubic form), adds the $2 \times 2$ minors of the Jacobian
      matrix, saturates, and returns the integer part.

      This invariant combines information from both the $c_2$ form and the
      cubic form.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      cubicLinearConductorInvariants X
  SeeAlso
    cubicConductorInvariants
    singularContents
///

doc ///
  Key
    singularContentsQuartic
    (singularContentsQuartic, CalabiYauInToric)
  Headline
    singular contents of the product of the c2 and cubic forms
  Usage
    L = singularContentsQuartic X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:List
      cumulative GCDs by degree of the singular locus of $L \cdot F$
  Description
    Text
      Computes @TO singularContents@ of the product $L \cdot F$ where
      $L$ is the $c_2$ form (divided by its content) and $F$ is the
      cubic form (divided by its content).  The product is a quartic,
      hence the name.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      singularContentsQuartic X
  SeeAlso
    singularContents
    cubicConductorInvariants
///

----------------------------------------------
-- Aronhold invariants -----------------------
----------------------------------------------

doc ///
  Key
    aronhold
    (aronhold, RingElement)
  Headline
    Aronhold invariants S and T of a ternary cubic form
  Usage
    L = aronhold F
  Inputs
    F:RingElement
      a homogeneous cubic polynomial in 3 variables
  Outputs
    L:List
      $\{S, T\}$ the two Aronhold invariants
  Description
    Text
      Computes the two Aronhold invariants $S$ (degree 4) and $T$ (degree 6)
      of a ternary cubic form.  These are classical invariants of plane cubic
      curves under the action of $SL(3)$.

      The $j$-invariant of the cubic curve can be recovered as
      $j = 1728 \cdot (4S)^3 / ((4S)^3 - T^2)$.

      For a Calabi-Yau 3-fold with $h^{1,1} = 3$, the cubic intersection
      form is a ternary cubic, so these invariants apply directly.
    Example
      R = ZZ[x,y,z];
      F = x^3 + y^3 + z^3 - 7*x*y*z
      aronhold F
    Text
      For a CY 3-fold with $h^{1,1} = 3$:
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      aronhold cubicForm X
  Caveat
    Only defined for cubic polynomials in exactly 3 variables.
  SeeAlso
    cubicForm
    hubschInvariants
///

----------------------------------------------
-- Point counting ----------------------------
----------------------------------------------

doc ///
  Key
    PointCounter
  Headline
    type for counting solutions of polynomial equations over finite fields
  Description
    Text
      A @TT "PointCounter"@ precomputes evaluation maps for a polynomial
      ring over a collection of finite fields.  This allows efficient
      repeated point counting for different polynomials in the same ring.

      Create one using @TO pointCounter@, then use @TO pointCounts@
      to count zeros.
  SeeAlso
    pointCounter
    pointCounts
///

doc ///
  Key
    pointCounter
    (pointCounter, Ring)
  Headline
    create a PointCounter for a polynomial ring
  Usage
    PC = pointCounter R
  Inputs
    R:Ring
      a polynomial ring over $\mathbb{Z}$
  Outputs
    PC:PointCounter
  Description
    Text
      Creates a @TO PointCounter@ for the ring $R$.  This precomputes
      all evaluation maps for the default set of finite fields
      (primes 2, 3, 5, 7, 11, 13 and some prime powers).

      The {\tt Projective} option (default {\tt true}) determines whether
      to count projective points (up to scaling) or affine points.
    Example
      R = ZZ[a,b,c];
      PC = pointCounter R
  SeeAlso
    PointCounter
    pointCounts
///

doc ///
  Key
    pointCounts
    (pointCounts, PointCounter, RingElement, RingElement)
    (pointCounts, PointCounter, CalabiYauInToric)
  Headline
    count zeros of the c2 and cubic forms over finite fields
  Usage
    L = pointCounts(PC, L, F)
    L = pointCounts(PC, X)
  Inputs
    PC:PointCounter
    L:RingElement
      the $c_2$ form
    F:RingElement
      the cubic form
    X:CalabiYauInToric
  Outputs
    L:List
      for each finite field, a triple $\{$zeros of $L$, zeros of $F$,
      zeros of both$\}$
  Description
    Text
      Counts the number of zeros of the $c_2$ form $L$ and cubic form $F$
      over each finite field in the @TO PointCounter@.  For each field,
      returns a list of three counts: zeros of $L$ alone, zeros of $F$ alone,
      and zeros of both $L$ and $F$ simultaneously.

      These point counts are invariants of the forms (preserved under
      invertible linear changes of variables over $\mathbb{Z}$) and can
      help distinguish topologically inequivalent CYs.

      The forms are first divided by their content before counting.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      PC = pointCounter R
      pointCounts(PC, X)
  SeeAlso
    PointCounter
    pointCounter
    invariants
///
