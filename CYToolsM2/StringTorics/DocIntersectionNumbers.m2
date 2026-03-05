--------------------------------------------------
-- Intersection numbers on Calabi-Yau varieties --
--------------------------------------------------

doc ///
  Key
    "working with intersection numbers and intersection rings"
  Headline
    facilities in StringTorics for intersection numbers
  Description
    Text
      Given a @TO CalabiYauInToric@ $X$, a Calabi-Yau 3-fold hypersurface in a
      simplicial toric 4-fold, we can compute triple intersection numbers
      of divisors on $X$.  These intersection numbers, together with
      the second Chern class $c_2(X)$, determine the intersection ring of $X$.
    Text
      @SUBSECTION "Intersection numbers"@
    Text
      @UL {
          TO (intersectionNumbers, CalabiYauInToric),
          TO (toricIntersectionNumbers, CalabiYauInToric),
          TO (intersectionNumbersOfCY, CalabiYauInToric)
      }@
    Text
      @SUBSECTION "Intersection forms and c2"@
    Text
      @UL {
          TO (cubicForm, CalabiYauInToric),
          TO (c2, CalabiYauInToric),
          TO (c2Form, CalabiYauInToric)
      }@
    Text
      @SUBSECTION "Mori cone"@
    Text
      @UL {
          TO (toricMoriCone, CalabiYauInToric),
          TO (toricMoriConeCap, CalabiYauInToric)
      }@
  SeeAlso
    "creating and using Batyrev Calabi-Yau hypersurfaces"
///

doc ///
  Key
    intersectionNumbers
    (intersectionNumbers, CalabiYauInToric)
  Headline
    triple intersection numbers of a Calabi-Yau 3-fold in the Picard basis
  Usage
    L = intersectionNumbers X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:List
      of pairs $\{i,j,k\} \Rightarrow n$, where $0 \le i \le j \le k$ are indices
      into the Picard basis and $n$ is the triple intersection number
  Description
    Text
      Returns the non-zero triple intersection numbers $D_i \cdot D_j \cdot D_k$
      of divisors on $X$, expressed in the Picard basis given by @TO basisIndices@.
      Only non-zero entries are returned.

      These are the intersection numbers in the basis for the Picard group, not the
      full toric intersection numbers (see @TO (toricIntersectionNumbers, CalabiYauInToric)@).
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      Xs = findAllCYs(Q, PicardRing => R)
      intersectionNumbers Xs#0
  SeeAlso
    cubicForm
    toricIntersectionNumbers
    basisIndices
///

doc ///
  Key
    toricIntersectionNumbers
    (toricIntersectionNumbers, CalabiYauInToric)
  Headline
    triple intersection numbers using all toric divisors
  Usage
    L = toricIntersectionNumbers X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:List
      of pairs $\{i,j,k\} \Rightarrow n$, where $i,j,k$ are indices into @TT "rays Q"@
      and $n$ is the triple intersection number
  Description
    Text
      Returns the non-zero triple intersection numbers $D_i \cdot D_j \cdot D_k$
      of toric divisors on $X$, using indices into the full set of rays of the
      underlying reflexive polytope, not just the Picard basis.

      Compare with @TO (intersectionNumbers, CalabiYauInToric)@, which returns
      the intersection numbers in the Picard basis.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      toricIntersectionNumbers X
  SeeAlso
    intersectionNumbers
    cubicForm
///

doc ///
  Key
    c2
    (c2, CalabiYauInToric)
  Headline
    second Chern class values on the Picard basis
  Usage
    L = c2 X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:List
      of integers, the values $c_2(X) \cdot D_i$ for each basis divisor $D_i$
  Description
    Text
      Returns the list of integers $\{c_2(X) \cdot D_0, c_2(X) \cdot D_1, \ldots\}$
      where $D_i$ are the divisors in the Picard basis given by @TO basisIndices@.

      Together with the triple intersection numbers, this determines the
      topological type of the Calabi-Yau (up to diffeomorphism equivalences).
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      c2 X
      c2Form X
  SeeAlso
    c2Form
    intersectionNumbers
    cubicForm
///

doc ///
  Key
    c2Form
    (c2Form, CalabiYauInToric)
  Headline
    second Chern class as a linear form in the Picard ring
  Usage
    f = c2Form X
  Inputs
    X:CalabiYauInToric
  Outputs
    f:RingElement
      a linear polynomial in the @TO picardRing@ of $X$
  Description
    Text
      Returns $c_2(X)$ expressed as a linear form $\sum c_i \, a_i$ in the
      @TO picardRing@ of $X$, where $c_i = c_2(X) \cdot D_i$.

      This is essentially the same data as @TO (c2, CalabiYauInToric)@,
      but expressed as a ring element.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      c2Form X
      ring c2Form X === R
  SeeAlso
    c2
    cubicForm
    picardRing
///

doc ///
  Key
    cubicForm
    (cubicForm, CalabiYauInToric)
  Headline
    cubic intersection form as a polynomial in the Picard ring
  Usage
    f = cubicForm X
  Inputs
    X:CalabiYauInToric
  Outputs
    f:RingElement
      a cubic polynomial in the @TO picardRing@ of $X$
  Description
    Text
      Returns the cubic intersection form $\sum d_{ijk} \, a_i \, a_j \, a_k$
      in the @TO picardRing@ of $X$.  The coefficient of the monomial
      $a_i \, a_j \, a_k$ (with $i \le j \le k$) is the multinomial coefficient
      times the triple intersection number $D_i \cdot D_j \cdot D_k$.

      This encodes the same data as @TO (intersectionNumbers, CalabiYauInToric)@
      but as a polynomial.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      Xs = findAllCYs(Q, PicardRing => R)
      Xs/cubicForm
      Xs/(X -> c2Form X)
  SeeAlso
    intersectionNumbers
    c2Form
    picardRing
///

doc ///
  Key
    intersectionNumbersOfCY
    (intersectionNumbersOfCY, CalabiYauInToric)
    (intersectionNumbersOfCY, NormalToricVariety, List)
    (intersectionNumbersOfCY, Ring, List)
  Headline
    intersection numbers via the abstract intersection ring (alternate method)
  Usage
    L = intersectionNumbersOfCY X
    L = intersectionNumbersOfCY(V, basisIndices)
    L = intersectionNumbersOfCY(IX, basisIndices)
  Inputs
    X:CalabiYauInToric
    V:NormalToricVariety
    IX:Ring
      the intersection ring from Schubert2
    basisIndices:List
      indices of the basis divisors
  Outputs
    L:List
      of pairs $\{i,j,k\} \Rightarrow n$
  Description
    Text
      An alternate, slower method for computing triple intersection numbers,
      using the abstract intersection ring from {\tt Schubert2}.  This is useful
      as a check on the primary method @TO (intersectionNumbers, CalabiYauInToric)@.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      intersectionNumbers X
      intersectionNumbersOfCY X
  SeeAlso
    intersectionNumbers
///

doc ///
  Key
    toricMoriCone
    (toricMoriCone, CalabiYauInToric)
    (toricMoriCone, NormalToricVariety, List)
  Headline
    the Mori cone of the ambient toric variety restricted to the Picard basis
  Usage
    C = toricMoriCone X
    C = toricMoriCone(V, basisIndices)
  Inputs
    X:CalabiYauInToric
    V:NormalToricVariety
    basisIndices:List
  Outputs
    C:Cone
  Description
    Text
      Computes the Mori cone (cone of effective curves) of the ambient toric
      variety, restricted to the Picard basis of the Calabi-Yau.  This is a
      polyhedral cone in $\RR^{h^{1,1}}$.

      Currently only works for favorable polytopes.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      C = toricMoriCone X
      rays C
  SeeAlso
    toricMoriConeCap
///

doc ///
  Key
    toricMoriConeCap
    (toricMoriConeCap, CalabiYauInToric)
  Headline
    intersection of Mori cones over all triangulations
  Usage
    L = toricMoriConeCap X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:List
      of lists of integers, the rays of the intersection cone
  Description
    Text
      For a given reflexive polytope, different triangulations give different
      toric Mori cones.  This function computes the intersection (``cap'') of
      the Mori cones over all triangulations that are equivalent to $X$
      (i.e., have the same 2-face restriction).  The actual Mori cone of $X$ is
      contained in this one.

      The result is returned as a list of ray vectors (lists of integers).
      Returns  @TT "null"@ if $X$ is not favorable.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      toricMoriConeCap X
  SeeAlso
    toricMoriCone
///
