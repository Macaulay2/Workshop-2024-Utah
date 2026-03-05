------------------------------------------------------------
-- Line bundle cohomology on toric varieties              --
-- and complete intersections                             --
------------------------------------------------------------

doc ///
  Key
    "working with cohomology of line bundles and sheaves in torics"
  Headline
    computing line bundle cohomology on toric varieties and their subvarieties
  Description
    Text
      StringTorics provides several methods for computing cohomology of
      line bundles on normal toric varieties and on complete intersections
      (including Calabi-Yau hypersurfaces) in toric varieties.
    Text
      @SUBSECTION "Cohomology on toric varieties"@
    Text
      The basic approach represents cohomology groups $H^i(V, \mathcal{O}_V(D))$ via
      Laurent monomials in the fraction field of the Cox ring.  The key functions are:

      @UL {
          TO (toricCohomologySetup, NormalToricVariety),
          TO (cohomologyBasis, ZZ, NormalToricVariety, List),
          TO (toricOrthants, NormalToricVariety),
          TO (cohomologyFractions, NormalToricVariety, List, List)
      }@
    Text
      @SUBSECTION "Cohomology matrices"@
    Text
      Given a line bundle on $V$ and a section $F$ of $\mathcal{O}_V(-K_V)$,
      the multiplication map $H^i(V, \mathcal{O}_V(D+K_V)) \to H^i(V, \mathcal{O}_V(D))$
      can be expressed as a matrix.  These matrices are used to compute
      cohomology on hypersurfaces via long exact sequences.

      @UL {
          TO (cohomologyMatrix, ZZ, NormalToricVariety, List, RingElement),
          TO (cohomologyMatrixRank, ZZ, NormalToricVariety, List, RingElement),
          TO (genericCohomologyMatrix, ZZ, NormalToricVariety, List)
      }@
    Text
      @SUBSECTION "Cohomology on complete intersections"@
    Text
      For Calabi-Yau hypersurfaces and complete intersections in toric varieties:

      @UL {
          TO (cohomologyVector, NormalToricVariety, List),
          TO (cohomologyVector, CompleteIntersectionInToric, List),
          TO (cohomologyVector, CalabiYauInToric, List, RingElement),
          TO (cohomologyOmega1, CompleteIntersectionInToric),
          TO (hodgeDiamond, CompleteIntersectionInToric),
          TO (hodgeVector, List)
      }@
    Text
      @SUBSECTION "Degree normalization utilities"@
    Text
      @UL {
          TO (normalDegrees, NormalToricVariety),
          TO (normalDegree, NormalToricVariety, List),
          TO (setOfColumns, NormalToricVariety, List),
          TO (getFraction, NormalToricVariety, List),
          TO (findTope, Matrix, List, List)
      }@
  SeeAlso
    "working with complete intersections in toric varieties"
    "working with intersection numbers and intersection rings"
///

----------------------------------------------
-- Cohomology setup --------------------------
----------------------------------------------

doc ///
  Key
    toricCohomologySetup
    (toricCohomologySetup, NormalToricVariety)
    CohomologySetup
  Headline
    set up cohomology data for a toric variety
  Usage
    H = toricCohomologySetup V
  Inputs
    V:NormalToricVariety
  Outputs
    H:HashTable
      mapping cohomological degree $i$ to lists of (orthant, ring) pairs
  Description
    Text
      Initializes and caches the data needed to compute cohomology of
      line bundles on $V$ using Laurent monomials.  The result is a hash
      table indexed by cohomological degree $i$, where each value is a list
      of pairs describing which ``orthants'' in $\mathbb{Z}^n$ support
      $H^i$.

      This function is called automatically by @TO cohomologyBasis@
      and does not usually need to be called directly.
    Example
      V = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/32003)
      H = toricCohomologySetup V
      keys H
  SeeAlso
    cohomologyBasis
    toricOrthants
///

doc ///
  Key
    cohomologyBasis
    (cohomologyBasis, ZZ, NormalToricVariety, List)
  Headline
    basis of a line bundle cohomology group as Laurent monomials
  Usage
    B = cohomologyBasis(i, V, D)
  Inputs
    i:ZZ
      the cohomological degree
    V:NormalToricVariety
    D:List
      the multidegree of the line bundle $\mathcal{O}_V(D)$
  Outputs
    B:List
      of Laurent monomials in the fraction field of the Cox ring
  Description
    Text
      Returns a basis of $H^i(V, \mathcal{O}_V(D))$ represented as
      Laurent monomials in the fraction field of the Cox ring of $V$.
    Example
      V = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/32003)
      cohomologyBasis(0, V, {1, 0, 0})
      cohomologyBasis(1, V, {2, 1, 0})
  SeeAlso
    toricCohomologySetup
    cohomologyMatrix
///

----------------------------------------------
-- Toric orthants ----------------------------
----------------------------------------------

doc ///
  Key
    toricOrthants
    (toricOrthants, NormalToricVariety)
  Headline
    orthants in the degree lattice supporting cohomology
  Usage
    L = toricOrthants V
  Inputs
    V:NormalToricVariety
  Outputs
    L:List
      of pairs $(i, S)$ where $i$ is the cohomological degree and
      $S$ is a subset of $\{0, \ldots, n-1\}$
  Description
    Text
      Returns a list describing which ``orthants'' in $\mathbb{Z}^n$ (where $n$
      is the number of rays) can support nonzero cohomology.  Each entry
      $(i, S)$ means that $H^i$ can be nonzero when the degree has negative
      entries in the positions indexed by $S$.
    Example
      V = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/32003)
      netList toricOrthants V
  SeeAlso
    toricCohomologySetup
    H1orthants
///

doc ///
  Key
    H1orthants
    (H1orthants, NormalToricVariety)
  Headline
    orthants supporting H1 cohomology
  Usage
    L = H1orthants V
  Inputs
    V:NormalToricVariety
  Outputs
    L:List
      of subsets of $\{0, \ldots, n-1\}$
  Description
    Text
      Returns the subsets of ray indices whose corresponding orthants
      can support $H^1$ cohomology.  This is a faster computation than
      @TO toricOrthants@ when only $H^1$ information is needed.
    Example
      V = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/32003)
      H1orthants V
  SeeAlso
    toricOrthants
///

----------------------------------------------
-- Cohomology matrices -----------------------
----------------------------------------------

doc ///
  Key
    cohomologyMatrix
    (cohomologyMatrix, ZZ, NormalToricVariety, List, RingElement)
  Headline
    matrix of a multiplication map on cohomology
  Usage
    (M, B2, B1) = cohomologyMatrix(i, V, D, F)
  Inputs
    i:ZZ
      the cohomological degree
    V:NormalToricVariety
    D:List
      a multidegree
    F:RingElement
      an element of the Cox ring
  Outputs
    M:MutableMatrix
      the multiplication map
    B2:List
      basis of the target
    B1:List
      basis of the source
  Description
    Text
      Computes the matrix of the multiplication-by-$F$ map
      $H^i(V, \mathcal{O}_V(D - \deg F)) \to H^i(V, \mathcal{O}_V(D))$
      with respect to the Laurent monomial bases of source and target.
    Example
      P2 = toricProjectiveSpace 2
      V = P2 ** P2
      S = ring V
      F = random({3,3}, S)
      cohomologyMatrix(2, V, {-3, 6}, F)
  SeeAlso
    cohomologyBasis
    cohomologyMatrixRank
///

doc ///
  Key
    cohomologyMatrixRank
    (cohomologyMatrixRank, ZZ, NormalToricVariety, List, RingElement)
  Headline
    dimensions and rank of a cohomology multiplication map
  Usage
    (nr, nc, rk) = cohomologyMatrixRank(i, V, D, F)
  Inputs
    i:ZZ
    V:NormalToricVariety
    D:List
    F:RingElement
  Outputs
    nr:ZZ
      number of rows (dimension of target)
    nc:ZZ
      number of columns (dimension of source)
    rk:ZZ
      rank of the multiplication matrix
  Description
    Text
      Returns the number of rows, columns, and rank of the
      multiplication matrix computed by @TO cohomologyMatrix@.
      This is the key computation for determining line bundle
      cohomology on hypersurfaces via long exact sequences.
  SeeAlso
    cohomologyMatrix
    cohomologyVector
///

doc ///
  Key
    genericCohomologyMatrix
    (genericCohomologyMatrix, ZZ, NormalToricVariety, List)
  Headline
    cohomology multiplication matrix for a random anticanonical section
  Usage
    M = genericCohomologyMatrix(i, V, D)
  Inputs
    i:ZZ
    V:NormalToricVariety
    D:List
  Outputs
    M:MutableMatrix
  Description
    Text
      Computes the cohomology multiplication matrix for a random section
      of $\mathcal{O}_V(-K_V)$.  This gives the generic rank of the
      multiplication map, which determines the cohomology of a generic
      Calabi-Yau hypersurface.
  Caveat
    The matrix is computed over a finite field ($\mathbb{Z}/32003$),
    so there is a very small probability that the rank could be
    lower than the true generic rank.
  SeeAlso
    cohomologyMatrix
    cohomologyMatrixRank
///

----------------------------------------------
-- Degree normalization ----------------------
----------------------------------------------

doc ///
  Key
    normalDegrees
    (normalDegrees, NormalToricVariety)
  Headline
    normalized degree matrix of a toric variety
  Usage
    A = normalDegrees V
  Inputs
    V:NormalToricVariety
  Outputs
    A:Matrix
      the degree matrix in normal form $[-I \mid C]$
  Description
    Text
      Computes and caches a normal form for the degree matrix of $V$.
      A maximal cone with determinant $\pm 1$ is found, and the degree
      matrix is transformed so that the corresponding columns become
      $-I$ (the negative identity).

      This normalization is used internally by @TO cohomologyFractions@
      and related functions.
  SeeAlso
    normalDegree
    setOfColumns
    getFraction
///

doc ///
  Key
    normalDegree
    (normalDegree, NormalToricVariety, List)
  Headline
    translate a degree to the normalized basis
  Usage
    d = normalDegree(V, deg)
  Inputs
    V:NormalToricVariety
    deg:List
      a degree in the original basis
  Outputs
    d:List
      the degree in the normalized basis
  SeeAlso
    normalDegrees
///

doc ///
  Key
    setOfColumns
    (setOfColumns, NormalToricVariety, List)
  Headline
    translate column indices to the normalized ordering
  Usage
    S = setOfColumns(V, cols)
  Inputs
    V:NormalToricVariety
    cols:List
      column indices in the original ordering
  Outputs
    S:List
      the corresponding indices in the normalized ordering
  SeeAlso
    normalDegrees
///

doc ///
  Key
    getFraction
    (getFraction, NormalToricVariety, List)
  Headline
    convert normalized exponent vector to a Laurent monomial
  Usage
    f = getFraction(V, mon)
  Inputs
    V:NormalToricVariety
    mon:List
      an exponent vector in normalized coordinates
  Outputs
    f:RingElement
      a Laurent monomial in the fraction field of the Cox ring
  SeeAlso
    normalDegrees
    cohomologyFractions
///

doc ///
  Key
    findTope
    (findTope, Matrix, List, List)
  Headline
    find the lattice polytope of Laurent monomials in a cohomology group
  Usage
    (P, H, R) = findTope(A, neg, deg)
  Inputs
    A:Matrix
      the normalized degree matrix
    neg:List
      indices of the ``negative'' directions (the orthant)
    deg:List
      the degree in the normalized basis
  Outputs
    P:Polyhedron
      the polytope of valid exponent vectors
  Description
    Text
      Given the normalized degree matrix, a set of negative directions
      (specifying which orthant we are in), and a degree, returns the
      polytope whose lattice points give the exponent vectors of
      Laurent monomials in the corresponding cohomology group.

      This is used internally by @TO cohomologyFractions@.
  SeeAlso
    cohomologyFractions
    normalDegrees
///

doc ///
  Key
    cohomologyFractions
    (cohomologyFractions, NormalToricVariety, List, List)
  Headline
    Laurent monomial basis of a cohomology group via normalized degrees
  Usage
    L = cohomologyFractions(V, negSet, deg)
  Inputs
    V:NormalToricVariety
    negSet:List
      the negative set (subset of ray indices specifying the orthant)
    deg:List
      the multidegree
  Outputs
    L:List
      of Laurent monomials in the Cox ring
  Description
    Text
      Computes a basis of the cohomology group corresponding to the
      given orthant and degree, using the normalized degree matrix.
      This is an alternative to @TO cohomologyBasis@ that works via
      lattice point enumeration.
  SeeAlso
    cohomologyBasis
    normalDegrees
    findTope
///

----------------------------------------------
-- Cohomology on complete intersections ------
----------------------------------------------

doc ///
  Key
    cohomologyVector
    (cohomologyVector, NormalToricVariety, List)
    (cohomologyVector, NormalToricVariety, ToricDivisor)
    (cohomologyVector, CompleteIntersectionInToric, List)
    (cohomologyVector, CompleteIntersectionInToric, ToricDivisor)
    (cohomologyVector, CompleteIntersectionInToric)
    (cohomologyVector, CalabiYauInToric, List, RingElement)
    (cohomologyVector, CompleteIntersectionInToric, List, RingElement)
    (cohomologyVector, LineBundle)
  Headline
    cohomology dimensions of a line bundle
  Usage
    v = cohomologyVector(V, D)
    v = cohomologyVector(X, D)
    v = cohomologyVector(X, D, F)
    v = cohomologyVector L
  Inputs
    V:NormalToricVariety
    X:CompleteIntersectionInToric
      or a @TO CalabiYauInToric@
    D:List
      a multidegree (or a @TO ToricDivisor@)
    F:RingElement
      a defining equation (element of the Cox ring)
    L:LineBundle
  Outputs
    v:List
      of cohomology dimensions $(h^0, h^1, \ldots, h^{\dim})$
  Description
    Text
      Computes the cohomology vector $(h^0, h^1, \ldots, h^d)$ of the
      line bundle $\mathcal{O}(D)$.

      On a toric variety $V$, this uses CohomCalg.

      On a complete intersection or CY hypersurface, this uses
      multiplication matrices and long exact sequences.  The version
      taking a @TO RingElement@ $F$ uses the specific polynomial $F$
      rather than a generic one.
    Example
      V = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/32003)
      cohomologyVector(V, {1, 0, 0})
      cohomologyVector(V, {-2, 3, -4})
  SeeAlso
    cohomologyBasis
    cohomologyMatrix
    hodgeDiamond
///

doc ///
  Key
    cohomologyFromLES
    (cohomologyFromLES, Matrix)
  Headline
    compute cohomology from a long exact sequence
  Usage
    M' = cohomologyFromLES M
  Inputs
    M:Matrix
      a matrix whose entries form a long exact sequence
  Outputs
    M':Matrix
      the matrix reduced modulo the LES relations
  Description
    Text
      Given a matrix whose entries represent the terms of a long exact
      sequence (with zeros indicating where the sequence splits),
      computes the alternating sum relations and reduces the matrix
      modulo these relations.

      This is used internally by the cohomology computation for
      complete intersections.
  SeeAlso
    cohomologyVector
///

doc ///
  Key
    collectLineBundles
    (collectLineBundles, NormalToricVariety, List, List)
  Headline
    collect line bundles needed for complete intersection cohomology
  Usage
    L = collectLineBundles(V, CI, Ds)
  Inputs
    V:NormalToricVariety
    CI:List
      multidegrees of the complete intersection equations
    Ds:List
      multidegrees of the line bundles to compute
  Outputs
    L:List
      of multidegrees of all line bundles on $V$ whose cohomology
      is needed
  Description
    Text
      Given a complete intersection defined by equations of degrees
      in the list {\tt CI}, and a set of line bundle degrees {\tt Ds}
      on the complete intersection, determines all line bundle degrees
      on the ambient toric variety $V$ that must be computed via the
      Koszul complex.
  SeeAlso
    cohomologyVector
    basicCohomologies
///

doc ///
  Key
    basicCohomologies
    (basicCohomologies, CompleteIntersectionInToric)
  Headline
    compute basic cohomologies of a complete intersection
  Usage
    H = basicCohomologies X
  Inputs
    X:CompleteIntersectionInToric
  Outputs
    H:MutableHashTable
      cached cohomology data
  Description
    Text
      Computes cohomology vectors for the structure sheaf, all toric
      divisors restricted to $X$, and the anticanonical divisors.
      These are the basic cohomologies needed for further computations
      (e.g., $\Omega^1_X$ and the Hodge diamond).

      The results are cached in {\tt X.cache.cohom}.
  SeeAlso
    cohomologyVector
    cohomologyOmega1
    hodgeDiamond
///

doc ///
  Key
    cohomologyOmega1
    (cohomologyOmega1, CompleteIntersectionInToric)
  Headline
    cohomology of the cotangent sheaf of a complete intersection
  Usage
    v = cohomologyOmega1 X
  Inputs
    X:CompleteIntersectionInToric
  Outputs
    v:List
      cohomology vector $(h^0(\Omega^1_X), h^1(\Omega^1_X), \ldots)$
  Description
    Text
      Computes the cohomology of $\Omega^1_X$ using the conormal
      sequence and long exact sequences in cohomology.  This is
      needed for computing the Hodge diamond of $X$.
    Example
      V = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/32003)
      X = completeIntersection(V, {-toricDivisor V})
      cohomologyOmega1 X
  SeeAlso
    hodgeDiamond
    basicCohomologies
///

doc ///
  Key
    hodgeDiamond
    (hodgeDiamond, CompleteIntersectionInToric)
  Headline
    Hodge diamond of a complete intersection
  Usage
    M = hodgeDiamond X
  Inputs
    X:CompleteIntersectionInToric
  Outputs
    M:Matrix
      the Hodge diamond $h^{p,q}$
  Description
    Text
      Computes the Hodge diamond of a smooth complete intersection $X$
      in a toric variety, using cohomology of $\mathcal{O}_X$ and
      $\Omega^1_X$.  For $\dim X \leq 3$, all Hodge numbers are
      determined.  For $\dim X \geq 4$, entries requiring $\Omega^2_X$
      are returned as $-1$.
    Example
      V = smoothFanoToricVariety(3, 10, CoefficientRing => ZZ/32003)
      X = completeIntersection(V, {-toricDivisor V})
      hodgeDiamond X
  SeeAlso
    cohomologyOmega1
    hodgeVector
///

doc ///
  Key
    hodgeVector
    (hodgeVector, List)
    (hodgeVector, List, List)
  Headline
    cohomology vector of a complete intersection of toric divisors
  Usage
    v = hodgeVector Ds
    v = hodgeVector(Ds, Fs)
  Inputs
    Ds:List
      of @TO ToricDivisor@s on a toric variety
    Fs:List
      of polynomials in the Cox ring (for the codimension 2 version)
  Outputs
    v:List
      the cohomology vector of the complete intersection
  Description
    Text
      Computes the cohomology vector of the complete intersection
      defined by the given toric divisors.

      The first form uses generic cohomology (via CohomCalg and
      long exact sequences).  The second form takes specific
      polynomials defining the intersection (currently limited
      to codimension 2).
  SeeAlso
    hodgeDiamond
    cohomologyVector
///

doc ///
  Key
    removeUnusedVariables
    (removeUnusedVariables, List)
  Headline
    remove unused variables from a list of matrices
  Usage
    L' = removeUnusedVariables L
  Inputs
    L:List
      of matrices over a polynomial ring
  Outputs
    L':List
      the matrices over a ring with only the used variables
  Description
    Text
      Takes a list of matrices and returns them over a new ring
      containing only the variables that actually appear.  This is
      a utility function used in cohomology computations.
  SeeAlso
    cohomologyFromLES
///

doc ///
  Key
    nextBundle
    (nextBundle, CompleteIntersectionInToric)
  Headline
    allocate variables for the next unknown cohomology vector
  Usage
    v = nextBundle X
  Inputs
    X:CompleteIntersectionInToric
  Outputs
    v:List
      of fresh variables representing $(h^0, h^1, \ldots, h^{\dim X})$
  Description
    Text
      Allocates a new set of variables in the cohomology computation ring
      to represent an unknown cohomology vector.  This is used internally
      when building long exact sequences to solve for unknown cohomologies.
  SeeAlso
    basicCohomologies
    cohomologyOmega1
///
