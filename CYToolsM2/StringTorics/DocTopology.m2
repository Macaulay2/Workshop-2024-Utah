----------------------------------------------
-- Topology of Calabi-Yau 3-folds           --
----------------------------------------------

doc ///
  Key
    "topological equivalence of Calabi-Yau 3-folds"
  Headline
    facilities in StringTorics for topological equivalence of CY 3-folds
  Description
    Text
      The topological type of a Calabi-Yau 3-fold $X$ (up to diffeomorphism)
      is largely determined by the Hodge numbers $h^{1,1}$ and $h^{1,2}$,
      the cubic intersection form, and the second Chern class $c_2(X)$.
      Two Calabi-Yau 3-folds $X_1$ and $X_2$ with the same Hodge numbers are
      topologically equivalent if there exists an integer change-of-basis matrix $A$
      (with $\det A = \pm 1$) that maps the $c_2$ form and cubic form of $X_1$ to
      those of $X_2$.
    Text
      @SUBSECTION "Topological data"@
    Text
      @UL {
          TO TopologicalDataOfCY3,
          TO (topologicalData, CalabiYauInToric),
          TO (c2Form, TopologicalDataOfCY3),
          TO (cubicForm, TopologicalDataOfCY3),
          TO (hh, Sequence, TopologicalDataOfCY3)
      }@
    Text
      @SUBSECTION "Equivalence testing"@
    Text
      @UL {
          TO (isEquivalent, CalabiYauInToric, CalabiYauInToric, Matrix),
          TO (mapIsIsomorphism, Matrix, CalabiYauInToric, CalabiYauInToric),
          TO (invariants, CalabiYauInToric)
      }@
    Text
      @SUBSECTION "Partitioning into equivalence classes"@
    Text
      @UL {
          TO TopologySet,
          TO (topologySet, List, HashTable),
          TO (representatives, TopologySet),
          TO (equivalences, TopologySet),
          TO (separateIfDifferent, TopologySet, Function),
          TO (partitionByTopology, List)
      }@
  SeeAlso
    "working with intersection numbers and intersection rings"
    "creating and using Batyrev Calabi-Yau hypersurfaces"
///

doc ///
  Key
    TopologicalDataOfCY3
  Headline
    type storing the topological data of a Calabi-Yau 3-fold
  Description
    Text
      A @TT "TopologicalDataOfCY3"@ is a list containing
      the four pieces of data that (conjecturally) determine
      the diffeomorphism type of a Calabi-Yau 3-fold:
      the $c_2$ form (a linear polynomial),
      the cubic intersection form (a cubic polynomial),
      $h^{1,1}$, and $h^{1,2}$.

      These are accessed via @TO (c2Form, TopologicalDataOfCY3)@,
      @TO (cubicForm, TopologicalDataOfCY3)@, and
      @TO (hh, Sequence, TopologicalDataOfCY3)@.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      T = topologicalData X
      c2Form T
      cubicForm T
      hh^(1,1) T
      hh^(1,2) T
  SeeAlso
    topologicalData
    c2Form
    cubicForm
///

doc ///
  Key
    topologicalData
    (topologicalData, CalabiYauInToric)
  Headline
    extract topological data from a Calabi-Yau 3-fold
  Usage
    T = topologicalData X
  Inputs
    X:CalabiYauInToric
  Outputs
    T:TopologicalDataOfCY3
  Description
    Text
      Returns a @TO TopologicalDataOfCY3@ containing the $c_2$ form,
      cubic intersection form, $h^{1,1}$, and $h^{1,2}$ of $X$.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      T = topologicalData X
      c2Form T
      cubicForm T
  SeeAlso
    TopologicalDataOfCY3
    c2Form
    cubicForm
    isEquivalent
///

doc ///
  Key
    (c2Form, TopologicalDataOfCY3)
  Headline
    the c2 form stored in topological data
  Usage
    f = c2Form T
  Inputs
    T:TopologicalDataOfCY3
  Outputs
    f:RingElement
      a linear polynomial in the Picard ring
  Description
    Text
      Returns the $c_2$ form (a linear polynomial) from the topological data.
      This is the same polynomial as @TO (c2Form, CalabiYauInToric)@.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      T = topologicalData X
      c2Form T == c2Form X
  SeeAlso
    TopologicalDataOfCY3
    (c2Form, CalabiYauInToric)
///

doc ///
  Key
    (cubicForm, TopologicalDataOfCY3)
  Headline
    the cubic form stored in topological data
  Usage
    f = cubicForm T
  Inputs
    T:TopologicalDataOfCY3
  Outputs
    f:RingElement
      a cubic polynomial in the Picard ring
  Description
    Text
      Returns the cubic intersection form from the topological data.
      This is the same polynomial as @TO (cubicForm, CalabiYauInToric)@.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      T = topologicalData X
      cubicForm T == cubicForm X
  SeeAlso
    TopologicalDataOfCY3
    (cubicForm, CalabiYauInToric)
///

doc ///
  Key
    (hh, Sequence, TopologicalDataOfCY3)
  Headline
    Hodge numbers of topological data
  Usage
    n = hh^(p,q) T
  Inputs
    (p,q):Sequence
    T:TopologicalDataOfCY3
  Outputs
    n:ZZ
  Description
    Text
      Returns the Hodge number $h^{p,q}$ stored in the topological data.
      For a Calabi-Yau 3-fold, the interesting Hodge numbers are
      $h^{1,1}$ and $h^{1,2}$ (the rest are determined by these).
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      T = topologicalData X
      hh^(1,1) T
      hh^(1,2) T
      hh^(1,1) T == hh^(1,1) X
  SeeAlso
    TopologicalDataOfCY3
    (hh, Sequence, CalabiYauInToric)
///

doc ///
  Key
    isEquivalent
    (isEquivalent, CalabiYauInToric, CalabiYauInToric, Matrix)
    (isEquivalent, Sequence, Sequence, Matrix)
  Headline
    test whether two Calabi-Yau 3-folds are topologically equivalent via a given matrix
  Usage
    b = isEquivalent(X1, X2, A)
    b = isEquivalent((L1,F1), (L2,F2), A)
  Inputs
    X1:CalabiYauInToric
    X2:CalabiYauInToric
    A:Matrix
      an $h^{1,1} \times h^{1,1}$ integer matrix with $\det A = \pm 1$
  Outputs
    b:Boolean
      whether $A$ maps the topological data of $X_1$ to that of $X_2$
  Description
    Text
      Two Calabi-Yau 3-folds are topologically equivalent if there is an
      integer change-of-basis matrix $A$ (with $\det A = \pm 1$) such that
      applying the substitution $A$ to the $c_2$ form and cubic form of $X_1$
      gives those of $X_2$.

      This function tests whether a specific matrix $A$ gives such an equivalence.
      If $A$ is the identity matrix, this tests whether $X_1$ and $X_2$ have
      literally the same topological data (in the same basis).

      The second form takes pairs $(L_i, F_i)$ where $L_i$ is the $c_2$ form and
      $F_i$ is the cubic form, both in the same polynomial ring.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      A = id_(ZZ^3)
      isEquivalent(X, X, A)
  SeeAlso
    mapIsIsomorphism
    topologicalData
///

doc ///
  Key
    mapIsIsomorphism
    (mapIsIsomorphism, Matrix, CalabiYauInToric, CalabiYauInToric)
  Headline
    test whether a matrix gives a topological isomorphism between two Calabi-Yau 3-folds
  Usage
    b = mapIsIsomorphism(M, X1, X2)
  Inputs
    M:Matrix
      an $h^{1,1} \times h^{1,1}$ matrix
    X1:CalabiYauInToric
    X2:CalabiYauInToric
  Outputs
    b:Boolean
  Description
    Text
      Tests whether the matrix $M$, when used as a substitution on the
      generators of the Picard ring, maps the $c_2$ form and cubic form of
      $X_1$ to those of $X_2$.  This is equivalent to
      @TO (isEquivalent, CalabiYauInToric, CalabiYauInToric, Matrix)@
      but also works for @TT "CYToolsCY3"@ objects.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      M = id_(ZZ^3)
      mapIsIsomorphism(M, X, X)
  SeeAlso
    isEquivalent
///

doc ///
  Key
    invariants
    (invariants, CalabiYauInToric)
  Headline
    compute a list of topological invariants of a Calabi-Yau 3-fold
  Usage
    L = invariants X
  Inputs
    X:CalabiYauInToric
  Outputs
    L:List
      a list of invariants: $\{h^{1,1}, h^{1,2}$, point counts, bad primes,
      content of $c_2$, content of cubic form, number of factors, dimension, number of components$\}$
  Description
    Text
      Computes a list of numerical invariants of the topological data of $X$.
      Two topologically equivalent Calabi-Yau 3-folds must have the same invariants,
      so this can be used as a quick test to distinguish non-equivalent CYs.

      The invariants include Hodge numbers, point counts of the cubic form modulo
      small primes (2, 3, 5, 7, 11), the GCD of coefficients (``content'') of the
      $c_2$ and cubic forms, and algebraic-geometric properties of the cubic surface
      defined by the cubic form.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      invariants X
  SeeAlso
    isEquivalent
    separateIfDifferent
///

doc ///
  Key
    TopologySet
  Headline
    type for tracking topological equivalence classes
  Description
    Text
      A @TT "TopologySet"@ is a mutable hash table that organizes a collection
      of Calabi-Yau 3-folds into buckets of potentially equivalent topologies.
      It has two fields:

      @TT "\"Sets\""@: a list of buckets, where each bucket is a list of lists.
      Two CYs in different buckets are definitely topologically distinct.
      Within a bucket, CYs in the same inner list are known to be equivalent;
      CYs in different inner lists within the same bucket may or may not be equivalent.

      @TT "\"CYHash\""@: a hash table mapping labels to @TO CalabiYauInToric@ objects.

      The basic workflow is to create a @TO TopologySet@ from a list of labels and a hash table
      of CY objects, then refine the buckets using @TO separateIfDifferent@ with various
      invariant functions.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      labels = {label X}
      Xs = hashTable {label X => X}
      T = topologySet(labels, Xs)
      representatives T
  SeeAlso
    topologySet
    representatives
    equivalences
    separateIfDifferent
///

doc ///
  Key
    topologySet
    (topologySet, List, HashTable)
  Headline
    create a TopologySet from labels and CY objects
  Usage
    T = topologySet(labels, Xs)
  Inputs
    labels:List
      of labels (typically strings or sequences)
    Xs:HashTable
      mapping each label to a @TO CalabiYauInToric@ object
  Outputs
    T:TopologySet
  Description
    Text
      Creates a new @TO TopologySet@ with all labels in a single bucket
      (no information about equivalences yet).  The set can then be refined
      using @TO separateIfDifferent@ with various invariant functions.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      labels = {label X}
      Xs = hashTable {label X => X}
      T = topologySet(labels, Xs)
  SeeAlso
    TopologySet
    separateIfDifferent
///

doc ///
  Key
    representatives
    (representatives, TopologySet)
    [representatives, IgnoreSingles]
  Headline
    representatives of topological equivalence classes
  Usage
    L = representatives T
  Inputs
    T:TopologySet
  Outputs
    L:List
      of lists of labels, one per equivalence class
  Description
    Text
      Returns a list of lists, where each inner list gives the labels
      of CYs in one bucket of the @TO TopologySet@.  By default,
      singleton buckets (those with only one class) are omitted; use
      @TT "IgnoreSingles => false"@ to include them.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      labels = {label X}
      Xs = hashTable {label X => X}
      T = topologySet(labels, Xs)
      representatives(T, IgnoreSingles => false)
  SeeAlso
    TopologySet
    equivalences
    IgnoreSingles
///

doc ///
  Key
    equivalences
    (equivalences, TopologySet)
    [equivalences, IgnoreSingles]
  Headline
    equivalence maps stored in a TopologySet
  Usage
    H = equivalences T
  Inputs
    T:TopologySet
  Outputs
    H:HashTable
      mapping each representative label to a list of $\{$label, matrix$\}$ pairs
  Description
    Text
      Returns a hash table recording the known equivalences between CYs
      in the @TO TopologySet@.  Each key is a label, and its value is a
      (possibly empty) list of pairs $\{$label, matrix$\}$ indicating other
      CYs equivalent to it and the change-of-basis matrix giving the equivalence.

      By default, singleton classes are omitted; use
      @TT "IgnoreSingles => false"@ to include them.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      labels = {label X}
      Xs = hashTable {label X => X}
      T = topologySet(labels, Xs)
      equivalences(T, IgnoreSingles => false)
  SeeAlso
    TopologySet
    representatives
    IgnoreSingles
///

doc ///
  Key
    IgnoreSingles
  Headline
    option to ignore singleton equivalence classes
  Description
    Text
      An option for @TO representatives@ and @TO equivalences@.
      When @TT "IgnoreSingles => true"@ (the default), equivalence classes
      containing only a single element are omitted from the output.
      Set to @TT "false"@ to include all classes.
  SeeAlso
    representatives
    equivalences
///

doc ///
  Key
    separateIfDifferent
    (separateIfDifferent, TopologySet, Function)
  Headline
    refine a TopologySet by separating CYs with different invariant values
  Usage
    T' = separateIfDifferent(T, fun)
  Inputs
    T:TopologySet
    fun:Function
      a function from @TO CalabiYauInToric@ to some value
  Outputs
    T':TopologySet
      with finer buckets
  Description
    Text
      Refines the partition in a @TO TopologySet@ by splitting each bucket
      according to the values of the function @TT "fun"@.  If @TT "fun X1 =!= fun X2"@,
      then $X_1$ and $X_2$ are placed in different buckets (since they cannot be
      topologically equivalent).  If @TT "fun X1 === fun X2"@, no conclusion is drawn.

      This is the basic operation for narrowing down the topology partition:
      start with a @TO TopologySet@ and apply @TT "separateIfDifferent"@ with
      progressively finer invariants.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      labels = {label X}
      Xs = hashTable {label X => X}
      T = topologySet(labels, Xs)
      T' = separateIfDifferent(T, invariants)
      representatives(T', IgnoreSingles => false)
  SeeAlso
    TopologySet
    topologySet
    invariants
///

doc ///
  Key
    partitionByTopology
    (partitionByTopology, List)
    (partitionByTopology, List, HashTable, ZZ)
  Headline
    partition Calabi-Yau 3-folds by topological equivalence using GV invariants
  Usage
    H = partitionByTopology LGVs
    H = partitionByTopology(labels, Xs, degreeLimit)
  Inputs
    LGVs:List
      of pairs $X \Rightarrow \text{gvPartition}$
    labels:List
      of labels
    Xs:HashTable
      mapping labels to @TO CalabiYauInToric@ objects
    degreeLimit:ZZ
      degree limit for Gopakumar-Vafa invariant computation
  Outputs
    H:HashTable
      mapping representative labels to lists of $\{$label, matrix$\}$ pairs
  Description
    Text
      Partitions a list of Calabi-Yau 3-folds into topological equivalence classes
      using Gopakumar-Vafa (GV) invariants to find candidate change-of-basis matrices.

      The GV invariants determine linear maps between Picard groups that are
      compatible with the intersection pairing.  These candidate maps are then
      checked to see if they actually give a topological isomorphism (i.e., map
      $c_2$ and cubic forms correctly).

      The first form takes a pre-computed list of pairs $X \Rightarrow \text{gvPartition}$.
      The second form takes labels, a hash table of CYs, and a degree limit,
      and computes GV invariants internally.

      The result is a hash table whose keys are representative labels and
      whose values are lists of $\{$label, matrix$\}$ for the equivalent CYs.
  SeeAlso
    mapIsIsomorphism
    isEquivalent
    TopologySet
///
