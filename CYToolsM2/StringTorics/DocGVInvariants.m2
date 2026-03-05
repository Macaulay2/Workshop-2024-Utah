------------------------------------------------------------
-- Gopakumar-Vafa invariants and birational geometry      --
------------------------------------------------------------

doc ///
  Key
    "Gopakumar-Vafa invariants of (general in complex moduli) Calabi-Yau 3-folds"
  Headline
    computing and using GV invariants of CY 3-folds
  Description
    Text
      Gopakumar-Vafa (GV) invariants are integer-valued invariants of a
      Calabi-Yau 3-fold $X$ associated to curve classes.  For a generic
      complex structure on $X$, they count (in a BPS sense) the number of
      curves in each homology class.

      The GV invariants are computed by the external C++ program {\tt computeGV},
      which must be installed separately and whose path is specified
      in the package @TO Configuration@.
      Currently, GV computation is only supported for favorable CY 3-folds.
    Text
      @SUBSECTION "Computing GV invariants"@
    Text
      @UL {
          TO (gvInvariants, CalabiYauInToric),
          TO (gvCone, CalabiYauInToric),
          TO (gvInvariantsAndCone, CalabiYauInToric, ZZ),
          TO (partitionGVConeByGV, CalabiYauInToric)
      }@
    Text
      @SUBSECTION "GV data organized by rays"@
    Text
      @UL {
          TO GVTable,
          TO (gvTable, CalabiYauInToric),
          TO (gvTable, HashTable),
          TO (gvRays, GVTable),
          TO (gvRay, GVTable, List)
      }@
    Text
      @SUBSECTION "Extremal curve classification"@
    Text
      @UL {
          TO (extremalRayGVs, CalabiYauInToric, List),
          TO (classifyExtremalCurve, List),
          TO (classifyExtremalCurves, CalabiYauInToric),
          TO (isNilpotent, GVTable, List)
      }@
    Text
      @SUBSECTION "Finding linear maps between GV cones"@
    Text
      @UL {
          TO (findLinearMaps, HashTable, HashTable)
      }@
  SeeAlso
    "birational geometry of Calabi-Yaus"
    "working with intersection numbers and intersection rings"
///

doc ///
  Key
    "birational geometry of Calabi-Yaus"
  Headline
    flopping curves and the CY3 type
  Description
    Text
      The @TO CY3@ type represents a Calabi-Yau 3-fold together with
      its GV invariants and a record of which curves have been flopped.
      Starting from a @TO CalabiYauInToric@, one creates a @TO CY3@
      using @TO makeCY3@ and then performs flops using @TO performFlop@.

      A flop along a curve class $C$ with GV invariant $n$ changes the
      topology as follows: $c_2 \mapsto c_2 + 2n\ell$ and
      $F \mapsto F - n\ell^3$, where $\ell$ is the linear form
      corresponding to $C$.
    Text
      @UL {
          TO CY3,
          TO (makeCY3, CalabiYauInToric),
          TO (makeCY3, ZZ, ZZ, RingElement, RingElement),
          TO (performFlop, CY3, List),
          TO (negatedCurves, CY3),
          TO (moriCone, CY3),
          TO (c2Form, CY3),
          TO (cubicForm, CY3)
      }@
  SeeAlso
    "Gopakumar-Vafa invariants of (general in complex moduli) Calabi-Yau 3-folds"
///

----------------------------------------------
-- GV invariant computation ------------------
----------------------------------------------

doc ///
  Key
    gvInvariants
    (gvInvariants, CalabiYauInToric)
    (gvInvariants, NormalToricVariety, List)
    [gvInvariants, Mori]
    [gvInvariants, DegreeLimit]
    [gvInvariants, FilePrefix]
    [gvInvariants, Executable]
  Headline
    compute Gopakumar-Vafa invariants of a Calabi-Yau 3-fold
  Usage
    GV = gvInvariants X
    GV = gvInvariants(V, basisIndices)
  Inputs
    X:CalabiYauInToric
    V:NormalToricVariety
    basisIndices:List
      indices of divisors forming a basis for the Picard group
    Mori => List
      Hilbert basis generators of the Mori cone (computed if not given)
    DegreeLimit => ZZ
      degree bound for the computation (default: infinity)
    FilePrefix => String
      prefix for temporary input/output files
    Executable => String
      path to the {\tt computeGV} executable
  Outputs
    GV:HashTable
      mapping curve class sequences to GV invariants (integers)
  Description
    Text
      Computes Gopakumar-Vafa invariants by calling the external C++ program
      {\tt computeGV}.  The result is a hash table whose keys are curve classes
      (as sequences of integers) and whose values are the corresponding GV
      invariants.

      This function requires a favorable CY 3-fold and returns {\tt null}
      for non-favorable ones.

      The {\tt DegreeLimit} option controls how far the computation goes;
      larger values give more curve classes but take longer.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      GV = gvInvariants(X, DegreeLimit => 10)
      sort pairs GV
  Caveat
    Requires the external program {\tt computeGV} to be installed and
    configured.  Only works for favorable CY 3-folds.
  SeeAlso
    gvTable
    gvCone
    partitionGVConeByGV
///

doc ///
  Key
    gvCone
    (gvCone, CalabiYauInToric)
    [gvCone, Mori]
    [gvCone, DegreeLimit]
    [gvCone, FilePrefix]
    [gvCone, Executable]
  Headline
    compute the cone generated by curve classes with nonzero GV invariants
  Usage
    C = gvCone X
  Inputs
    X:CalabiYauInToric
    DegreeLimit => ZZ
  Outputs
    C:Cone
  Description
    Text
      Computes GV invariants of $X$ and returns the cone generated by all
      curve classes that have nonzero GV invariants.  For a generic complex
      structure, this is the Mori cone.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      C = gvCone(X, DegreeLimit => 10)
      rays C
  Caveat
    Requires {\tt computeGV}.  Returns {\tt null} for non-favorable CYs.
  SeeAlso
    gvInvariants
    toricMoriCone
///

doc ///
  Key
    gvInvariantsAndCone
    (gvInvariantsAndCone, CalabiYauInToric, ZZ)
    [gvInvariantsAndCone, Mori]
    [gvInvariantsAndCone, DegreeLimit]
    [gvInvariantsAndCone, FilePrefix]
    [gvInvariantsAndCone, Executable]
  Headline
    compute GV invariants and the GV cone in two stages
  Usage
    (GV, C) = gvInvariantsAndCone(X, D)
  Inputs
    X:CalabiYauInToric
    D:ZZ
      initial degree bound for cone generation
    DegreeLimit => ZZ
      overall degree limit for GV computation
  Outputs
    GV:HashTable
      the GV invariants
    C:Cone
      the GV cone, possibly refined using higher-degree data
  Description
    Text
      Computes GV invariants and constructs the GV cone in two stages.
      First, curves up to degree $D$ are used to generate an initial cone.
      Then, higher-degree curves that lie outside this initial cone are
      added to refine it.

      This two-stage approach can give a better approximation of the
      true Mori cone when the {\tt DegreeLimit} is not large enough
      to capture all extremal rays at low degree.
  Caveat
    Requires {\tt computeGV}.  Returns {\tt null} for non-favorable CYs.
  SeeAlso
    gvInvariants
    gvCone
///

doc ///
  Key
    partitionGVConeByGV
    (partitionGVConeByGV, CalabiYauInToric)
    [partitionGVConeByGV, Mori]
    [partitionGVConeByGV, DegreeLimit]
    [partitionGVConeByGV, FilePrefix]
    [partitionGVConeByGV, Executable]
  Headline
    partition extremal rays of the GV cone by their GV invariant
  Usage
    H = partitionGVConeByGV X
  Inputs
    X:CalabiYauInToric
    DegreeLimit => ZZ
  Outputs
    H:HashTable
      mapping GV invariants to lists of extremal ray generators
  Description
    Text
      Computes GV invariants, finds the GV cone, and partitions the
      extremal rays of the cone according to their GV invariant value.
      This is useful for finding candidate linear maps between the
      Picard groups of two CY 3-folds via @TO findLinearMaps@.
  Caveat
    Requires {\tt computeGV}.  Returns {\tt null} for non-favorable CYs.
  SeeAlso
    gvInvariants
    gvCone
    findLinearMaps
///

----------------------------------------------
-- GVTable and ray-organized data ------------
----------------------------------------------

doc ///
  Key
    GVTable
  Headline
    type for GV invariants organized by primitive curve rays
  Description
    Text
      A @TT "GVTable"@ stores GV invariants organized by primitive curve
      classes.  Instead of a flat hash table mapping curve classes to
      invariants, a @TT "GVTable"@ maps each primitive curve $C$ to a
      list $\{n_C, n_{2C}, n_{3C}, \ldots\}$ of GV invariants along
      multiples of $C$, up to a degree bound.

      A @TT "GVTable"@ also stores the degree limit and heft vector used
      in the computation.

      Create one using @TO gvTable@.  Access the ray data via @TO gvRays@
      and @TO gvRay@.
  SeeAlso
    gvTable
    gvRays
    gvRay
///

doc ///
  Key
    gvTable
    (gvTable, CalabiYauInToric)
    (gvTable, HashTable)
    [gvTable, DegreeLimit]
    [gvTable, Heft]
  Headline
    create a GVTable from GV invariant data
  Usage
    G = gvTable X
    G = gvTable GVhash
  Inputs
    X:CalabiYauInToric
    GVhash:HashTable
      mapping curve class lists to GV invariants
    DegreeLimit => ZZ
      required: the degree bound used
    Heft => List
      required (for HashTable input): the heft vector
  Outputs
    G:GVTable
  Description
    Text
      Creates a @TO GVTable@ from either a @TO CalabiYauInToric@ (by computing
      GV invariants internally) or from a pre-computed hash table of GV invariants.

      For the hash table form, both {\tt DegreeLimit} and {\tt Heft}
      must be provided.  For the CalabiYauInToric form, the heft vector
      is computed automatically.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      G = gvTable(X, DegreeLimit => 10)
      gvRays G
  Caveat
    Requires {\tt computeGV} when called on a CalabiYauInToric.
  SeeAlso
    GVTable
    gvRays
    gvRay
    gvInvariants
///

doc ///
  Key
    gvRays
    (gvRays, GVTable)
    (gvRays, CY3)
  Headline
    the primitive curve rays and their GV invariants
  Usage
    H = gvRays G
    H = gvRays X3
  Inputs
    G:GVTable
    X3:CY3
  Outputs
    H:HashTable
      mapping primitive curve classes (lists) to lists of GV invariants
  Description
    Text
      Returns the hash table of GV invariants organized by primitive curve
      class.  Each key is a primitive curve class (a list of integers),
      and the corresponding value is a list $\{n_C, n_{2C}, n_{3C}, \ldots\}$
      of GV invariants for multiples of that curve.
  SeeAlso
    GVTable
    gvTable
    gvRay
///

doc ///
  Key
    gvRay
    (gvRay, GVTable, List)
    (gvRay, CY3, List)
    (gvRay, CalabiYauInToric, List)
    (gvRay, HashTable, List, ZZ, List)
  Headline
    GV invariants along multiples of a curve class
  Usage
    L = gvRay(G, C)
    L = gvRay(X3, C)
  Inputs
    G:GVTable
    X3:CY3
    C:List
      a primitive curve class
  Outputs
    L:List
      GV invariants $\{n_C, n_{2C}, n_{3C}, \ldots\}$
  Description
    Text
      Returns the list of GV invariants along multiples of the curve class $C$.
      If $C$ is not a primitive curve in the table, an empty list is returned.

      For a @TO CY3@ object, if $C$ has been negated by a flop, the curve
      is automatically negated before lookup.
  SeeAlso
    GVTable
    gvRays
///

doc ///
  Key
    gvByRay
    (gvByRay, HashTable, ZZ, List)
  Headline
    reorganize a GV hash table by primitive curve classes
  Usage
    H = gvByRay(GVhash, degreeLimit, heftVector)
  Inputs
    GVhash:HashTable
      mapping curve class lists to GV invariants
    degreeLimit:ZZ
    heftVector:List
  Outputs
    H:HashTable
      mapping primitive curves to lists of GV invariants along the ray
  Description
    Text
      Takes a flat hash table of GV invariants (as returned by @TO gvInvariants@)
      and reorganizes it so that each key is a primitive curve class and
      the value is a list of GV invariants for multiples of that curve.

      This is the internal function used by @TO gvTable@.
  SeeAlso
    gvTable
    gvRays
///

----------------------------------------------
-- Extremal curve classification -------------
----------------------------------------------

doc ///
  Key
    extremalRayGVs
    (extremalRayGVs, CalabiYauInToric, List)
  Headline
    GV invariants along an extremal curve ray
  Usage
    L = extremalRayGVs(X, C)
  Inputs
    X:CalabiYauInToric
    C:List
      a curve class
  Outputs
    L:List
      GV invariants along multiples of $C$
  Description
    Text
      Computes GV invariants specifically along the ray determined by the
      curve class $C$, using a targeted degree limit based on the heft of $C$.
      This is more efficient than computing all GV invariants when only
      one ray is needed.
  Caveat
    Requires {\tt computeGV}.
  SeeAlso
    gvRay
    classifyExtremalCurve
///

doc ///
  Key
    classifyExtremalCurve
    (classifyExtremalCurve, List)
    (classifyExtremalCurve, HashTable, List, ZZ, List)
  Headline
    classify an extremal curve by its GV ray
  Usage
    type = classifyExtremalCurve rayGVs
  Inputs
    rayGVs:List
      GV invariants along multiples of a curve
  Outputs
    type:List
      a pair $\{$type name, data$\}$
  Description
    Text
      Classifies an extremal curve based on the GV invariants along its ray.
      The possible types are:

      @UL {
          {"\"ZERO\" -- all GV invariants are zero (curve not effective)"},
          {"\"FLOP\" -- a floppable curve (first two GV values nonnegative, rest zero)"},
          {"\"TYPEIII0\" -- a type III contraction (one of first two GV values is $-2$)"},
          {"\"TYPEIIIg\" -- a type III contraction with genus (one of first two values is negative)"},
          {"\"TYPEII\" -- a type II contraction (nonzero GV values beyond degree 2)"},
          {"\"OTHER\" -- too few values to classify"}
      }@
  SeeAlso
    classifyExtremalCurves
    extremalRayGVs
///

doc ///
  Key
    classifyExtremalCurves
    (classifyExtremalCurves, CalabiYauInToric)
  Headline
    classify all extremal curves of a Calabi-Yau 3-fold
  Usage
    H = classifyExtremalCurves X
  Inputs
    X:CalabiYauInToric
  Outputs
    H:HashTable
      mapping classification types to lists of curve classes
  Description
    Text
      Classifies all extremal curves in the toric Mori cone cap of $X$
      by computing GV invariants along each ray and applying
      @TO classifyExtremalCurve@.  The result is a hash table partitioning
      the extremal curves by type.
  Caveat
    Requires {\tt computeGV}.  Only works for favorable CYs.
  SeeAlso
    classifyExtremalCurve
    toricMoriConeCap
///

doc ///
  Key
    isNilpotent
    (isNilpotent, GVTable, List)
    (isNilpotent, CY3, List)
  Headline
    test whether a curve class has eventually zero GV invariants
  Usage
    b = isNilpotent(G, C)
    b = isNilpotent(X3, C)
  Inputs
    G:GVTable
    X3:CY3
    C:List
      a curve class
  Outputs
    b:Boolean
  Description
    Text
      Tests whether the GV invariants along the ray of $C$ appear to
      be eventually zero (specifically, whether the last one or two computed
      values are zero).  A nilpotent curve can potentially be flopped.
  SeeAlso
    classifyExtremalCurve
    performFlop
///

----------------------------------------------
-- findLinearMaps ----------------------------
----------------------------------------------

doc ///
  Key
    findLinearMaps
    (findLinearMaps, HashTable, HashTable)
  Headline
    find linear maps between two GV cone partitions
  Usage
    Ms = findLinearMaps(gv1, gv2)
  Inputs
    gv1:HashTable
      result of @TO partitionGVConeByGV@ for the first CY
    gv2:HashTable
      result of @TO partitionGVConeByGV@ for the second CY
  Outputs
    Ms:List
      of matrices over QQ that map the generators of one to the other
  Description
    Text
      Given two GV cone partitions (as returned by @TO partitionGVConeByGV@),
      this function finds all linear maps that send the extremal rays of the
      first cone to those of the second, respecting the partition by GV
      invariant values.

      These candidate maps can then be checked using @TO mapIsIsomorphism@
      to determine if they give actual topological equivalences.

      The function returns an empty list if the partitions are incompatible
      (different keys or different numbers of rays per key).
  SeeAlso
    partitionGVConeByGV
    mapIsIsomorphism
    partitionByTopology
///

----------------------------------------------
-- CY3 type and flops -----------------------
----------------------------------------------

doc ///
  Key
    CY3
  Headline
    type for a Calabi-Yau 3-fold with GV data and flop history
  Description
    Text
      A @TT "CY3"@ object represents a Calabi-Yau 3-fold together with
      its GV invariant data (as a @TO GVTable@), a list of negated
      (flopped) curves, and derived quantities like the Mori cone and
      nef cone.

      Starting from a @TO CalabiYauInToric@, create a @TT "CY3"@ via
      @TO makeCY3@.  Flops are performed with @TO performFlop@, producing
      a new @TT "CY3"@ with updated topology and flop history.

      The following methods are available for @TT "CY3"@ objects:
      @TO (c2Form, CY3)@, @TO (cubicForm, CY3)@,
      @TO (gvTable, CY3)@, @TO (gvRays, CY3)@, @TO (gvRay, CY3, List)@,
      @TO (negatedCurves, CY3)@, @TO (moriCone, CY3)@,
      @TO (performFlop, CY3, List)@.
  SeeAlso
    makeCY3
    performFlop
    CalabiYauInToric
///

doc ///
  Key
    makeCY3
    (makeCY3, CalabiYauInToric)
    (makeCY3, ZZ, ZZ, RingElement, RingElement)
    [makeCY3, GVTable]
    [makeCY3, DegreeLimit]
    [makeCY3, NegatedCurves]
    [makeCY3, Label]
  Headline
    create a CY3 object from a CalabiYauInToric or from raw data
  Usage
    X3 = makeCY3 X
    X3 = makeCY3(h11, h12, c2, cubic)
  Inputs
    X:CalabiYauInToric
    h11:ZZ
    h12:ZZ
    c2:RingElement
      the $c_2$ form (linear polynomial)
    cubic:RingElement
      the cubic intersection form
    DegreeLimit => ZZ
      degree limit for GV computation (required for CalabiYauInToric input)
    GVTable => GVTable
      a pre-computed GV table (for raw data input)
    NegatedCurves => List
      list of negated curve classes (default: empty)
    Label => Thing
      an optional label
  Outputs
    X3:CY3
  Description
    Text
      Creates a @TO CY3@ object.  When given a @TO CalabiYauInToric@,
      GV invariants are computed automatically up to the given {\tt DegreeLimit}.
      When given raw topological data ($h^{1,1}$, $h^{1,2}$, $c_2$ form,
      cubic form), a @TO GVTable@ must be provided separately.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      X3 = makeCY3(X, DegreeLimit => 10)
      c2Form X3
      cubicForm X3
      negatedCurves X3
  Caveat
    The CalabiYauInToric form requires {\tt computeGV} and a favorable CY.
  SeeAlso
    CY3
    performFlop
    gvTable
///

doc ///
  Key
    (c2Form, CY3)
  Headline
    the c2 form of a CY3
  Usage
    L = c2Form X3
  Inputs
    X3:CY3
  Outputs
    L:RingElement
  Description
    Text
      Returns the $c_2$ form of the CY3 object.  After a flop,
      this will differ from the original by $2n\ell$ where $n$ is the
      GV invariant of the flopped curve and $\ell$ is its linear form.
  SeeAlso
    CY3
    (c2Form, CalabiYauInToric)
    performFlop
///

doc ///
  Key
    (cubicForm, CY3)
  Headline
    the cubic form of a CY3
  Usage
    F = cubicForm X3
  Inputs
    X3:CY3
  Outputs
    F:RingElement
  Description
    Text
      Returns the cubic intersection form of the CY3 object.  After a flop,
      this will differ from the original by $-n\ell^3$ where $n$ is the
      GV invariant of the flopped curve and $\ell$ is its linear form.
  SeeAlso
    CY3
    (cubicForm, CalabiYauInToric)
    performFlop
///

doc ///
  Key
    negatedCurves
    (negatedCurves, CY3)
  Headline
    list of curve classes that have been negated by flops
  Usage
    Cs = negatedCurves X3
  Inputs
    X3:CY3
  Outputs
    Cs:List
      of curve classes (lists of integers)
  Description
    Text
      Returns the list of primitive curve classes that have been
      negated by performing flops on the CY3 object.  These curves
      are used to adjust the Mori cone: a negated curve $C$ contributes
      $-C$ to the Mori cone instead of $C$.
  SeeAlso
    CY3
    performFlop
    moriCone
///

doc ///
  Key
    performFlop
    (performFlop, CY3, List)
  Headline
    perform a flop along a curve class
  Usage
    X3' = performFlop(X3, C)
  Inputs
    X3:CY3
    C:List
      a primitive curve class (with $\gcd = 1$)
  Outputs
    X3':CY3
      the flopped CY3
  Description
    Text
      Performs a flop of the Calabi-Yau 3-fold along the curve class $C$.
      If $n$ is the GV invariant of $C$ and $\ell = \sum C_i x_i$, then
      the new topological data is:
      $$c_2' = c_2 + 2n\ell, \qquad F' = F - n\ell^3.$$
      The curve $C$ is added to the list of negated curves (or removed
      if it was already negated, i.e., flopping back).

      The resulting @TO CY3@ shares the same @TO GVTable@ as the original.
    Example
      Q = reflexivePolytope(
          {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
           {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
          ID => 7)
      R = ZZ[a,b,c];
      X = makeCY(Q, PicardRing => R, ID => 0)
      X3 = makeCY3(X, DegreeLimit => 10)
      c2Form X3
      cubicForm X3
      -- Find a nilpotent (floppable) curve
      nilps = select(keys gvRays X3, c -> isNilpotent(X3, c))
      if #nilps > 0 then (
          C := first nilps;
          X3' := performFlop(X3, C);
          c2Form X3'
          )
  Caveat
    The curve class $C$ must be primitive (content 1).
  SeeAlso
    CY3
    makeCY3
    negatedCurves
    isNilpotent
///

doc ///
  Key
    moriCone
    (moriCone, CY3)
    (moriCone, GVTable, List)
  Headline
    the Mori cone of a CY3, accounting for flopped curves
  Usage
    C = moriCone X3
  Inputs
    X3:CY3
  Outputs
    C:Cone
  Description
    Text
      Returns the Mori cone of the @TO CY3@ object, computed from the
      primitive curves in the @TO GVTable@ with negated curves reversed.
      After a flop along curve $C$, the ray $C$ is replaced by $-C$
      in the cone.
  SeeAlso
    CY3
    negatedCurves
    performFlop
    toricMoriCone
///
