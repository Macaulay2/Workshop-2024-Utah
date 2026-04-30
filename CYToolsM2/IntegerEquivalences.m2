newPackage(
    "IntegerEquivalences",
    Version => "0.1",
    Date => "13 Jan 2024",
    Headline => "finding invertible integral matrices preserving points, linear forms and ideals",
    Authors => {{ Name => "", Email => "", HomePage => ""}},
    AuxiliaryFiles => false,
    DebuggingMode => true
    )

export {
    -- easy interface
    "findEquivalence",
    -- by hand interface
    "equivalenceIdeal",
    "factorsByType",
    "factorType",
    "idealsByBetti",
    "genericLinearMap", -- genericLinearMap(R).  Constructs two new rings, T, U, a matrix A over T nxn, n = numgens R, and phi = map(U, U, transpose A).
    "invertibleMatrixOverZZ",
    "hessianMatches",
    "selectLinear",
    "matchingData",
    "allSigns",
    "signedPermutations",
    "matches",
    "singularPoints",
    "singularPointMatches",
    "singularMatches",
    "tryEquivalences",
    -- types defined
    "MatchingData",
    -- utilities
    "extendToMatrix", -- extendToMatrix(List of integers) ==> Matrix (over ZZ).
    "cartesian",
    "hessian", -- place in Core?
    -- optional arguments and symbols usd in matching
    "RowVector",
    "ColumnVector",
    "Unknown",
    "SignedPermutations",
    "Permutations",
    -- symbols returnd
    "CONSISTENT",
    "INCONSISTENT",
    "INDETERMINATE"
    }

importFrom_"LLLBases"{"gcdLLL"};

extendToMatrix = method()
extendToMatrix List := Matrix => L -> (
    if not all(L, a -> instance(a, ZZ))
    then error "expected a list of integers";
    (g, A) := gcdLLL L; -- coming from LLLBases.
    transpose A
    -- TODO: should this insure that the matrix has determinant 1 (not -1)?
    -- I don't really need that...
    )

genericLinearMap = method(Options => {Variable => null})
genericLinearMap Ring := Sequence => opts -> R -> (
    -- R should be a polynomial ring in n >= 1 variables.
    n := numgens R;
    if n == 0 then (
        kk := coefficientRing R; -- TODO: if none, this should give a better error message
        A := map(kk^0, kk^0, {});
        return (A, id_R);
        );
    K := coefficientRing R;
    t := if opts.Variable === null then getSymbol "t" else opts.Variable;
    T := K[t_(1,1)..t_(n,n)];
    U := T [gens R, Join => false];
    A = map(T^n,,transpose genericMatrix(T, T_0, n, n));
    phi := map(U, U, transpose A);
    (A, phi)
    )

-- MatchingData: is a list of elements each of the form
--  L => {M0, M1, ..., Ms}
-- where L is a list of:
--   RingElement: a polynomial in the original ring RZ or RQ.
--   Ideal: an ideal in the original ring.
--   Matrix: either a row vector or column vector, over ZZ or QQ (or RZ or RQ)
-- and each list Mi has the same length as L, and the same types of its elements.
MatchingData = new Type of List

matchingData = method()
matchingData List := LMs -> (
    if not all(LMs, x -> instance(x, Option) or (instance(x, List) and #x === 3))
    then error "expected each element to be Item => Item, or a list {type, ItemList, ItemList}";
    ans := new MatchingData from LMs;
    if not isWellDefined ans then error "expected matching data to match.  Set `debugLevel=1` to investigate";
    ans
    )

MatchingData | MatchingData := MatchingData => (md1, md2) -> (
    matchingData join(toList md1, toList md2)
    )
MatchingData | Nothing := (MD, nothing) -> null
Nothing | MatchingData := (nothing, MD) -> null
Nothing | Nothing := (nothing1, nothing2) -> null

-- helper function for validMatchingItem.  
--   Input: either source or one of the targets of the matching data.
--   Output: a list of types.
itemType = L -> (
    -- L is a list or a single element of the following form
    -- returns a list of RingElement, Ideal, RowVector, ColumnVector.
    if not instance(L, List) then L = {L};
    for elem in L list (
        if instance(elem, RingElement) then RingElement
        else if instance(elem, Ideal) then Ideal
        else (
            if instance(elem, Matrix) then (
                if numrows elem === 1 then RowVector
                else if numcols elem === 1 then ColumnVector
                else Unknown
            ) else 
                Unknown
        ))
    )

-- helper function for (isWellDefined, MatchingData)
--   Input: one element (an Option, source and possible targets) of the matching data
--          which index this data sits at (used for error messages if debugLevel > 0)
--   Output: Boolean, whether this item is valid.
validMatchingItem = (LM,i) -> (
    -- TODO: if typ is null, then L should be a RingElement, Ideal or row or column Matrix, and M should be the same type
    if instance(LM, Option) then (
        if not (instance(LM#0, RingElement) or instance(LM#0, Ideal) or instance(LM#0, Matrix))
          or class LM#0 =!= class LM#1
          then (
              if debugLevel > 0 then 
                << "excepted each item of Option to be the same type: an Ideal, RingElement or row or column matrix" << endl;
                return false;
              );
          -- todo: make sure matrices are row or column matrices.
          return true;
        );
    (typ, L, M) := (LM#0, LM#1, LM#2);
    if typ =!= SignedPermutations and typ =!= Permutations then (
        if debugLevel > 0 then 
            << "expected first entry of list to be SignedPermutations or Permutations, instead, received "
                << toString typ << endl;
        return false;                
        );

    if not instance(L, List) or not instance(M, List) then (
        if debugLevel > 0 then 
          << "expected elements in source " << i << " to be lists" << endl;
          return false;
        );
    -- M should be a list, same length and types as L.
    Ltype := itemType L;
    Mtype := itemType M;
    if any(Ltype, x -> x === Unknown) then (
        if debugLevel > 0 then << "one element in source is unknown" << endl;
        return false;
        );
    if Mtype =!= Ltype then (
        if debugLevel > 0 then << "in source " << i << ", target doesn't match source = " << Ltype << " obtaining instead " << Mtype << endl;
        return false;
        );
    true
    )

isWellDefined MatchingData := Boolean => LMs -> (
    for i from 0 to #LMs-1 do (
        LM := LMs#i;
        if not instance(LM, Option) and not (instance(LM, List) and #LM === 3) then (
            if debugLevel > 0 then << "elements of list must be of the form {type, List, List}" << endl;
            return false;
            );
        if not validMatchingItem(LM,i) then return false;
        );
    true
    )

-- Used in creating MatchingData: use permutations or signedPermutations.
allSigns = method()
allSigns List := L -> (
    if #L <= 0 then return {{}};
    if #L == 1 then return {{L#0}, {-L#0}};
    flatten for q in allSigns(drop(L, 1)) list {prepend(L#0, q), prepend(-L#0, q)}
    )

signedPermutations = method()
signedPermutations List := List =>  L -> (
    flatten for p in permutations L list allSigns p
    )

cartesian = method()
cartesian List := (Ls) -> (
    -- cartesian product of Ls: one element from each
    -- so the result is a list of lists.
    if #Ls == 1 then return for p in Ls#0 list {p};
    Ls1 := cartesian drop(Ls,1);
    flatten for p in Ls#0 list for q in Ls1 list prepend(p, q)
    )

matches = method()
matches MatchingData := List => (MD) -> (
    -- for each element of the MatchingData, we make the list of all possible targets
    src := flatten for elem in MD list
        if instance(elem, Option) then elem#0 else elem#1;
    -- now create targets.  This means making all permutations, signed or not, and taking cartesian products.
    targs := cartesian for elem in MD list (
        if instance(elem, Option) then {elem#1}
        else if elem#0 === SignedPermutations then signedPermutations elem#2
        else if elem#0 === Permutations then permutations elem#2
        );
    targs = targs/flatten;
    (src, targs)
    )

selectLinear = method()
selectLinear MatchingData := MatchingData => (MD) -> (
    -- ASSUMPTION: the base ring for Ideals and RingElement's is standard graded polynomial ring.
    -- we keep only the ring elements that are linear.
    -- for each ideal, we take only the linear elements.
    -- We keep all row and column matrices.
    matchingData for elem in MD list (
        if instance(elem, Option) then (
            if instance(elem#0, RingElement) then (
                if degree elem#0 === {1} then elem else continue
            ) else if instance(elem#0, Ideal) then (
                elem0 := select(elem#0_*, f -> degree f === {1});
                if #elem0 === 0 then continue;
                elem1 := select(elem#1_*, f -> degree f === {1});
                ideal elem0 => ideal elem1
            ) else if instance(elem#0, Matrix) then elem
        ) else (
            if instance(elem#1#0, RingElement) then (
                if degree elem#1#0 === {1} then elem else continue
                )
            else if instance(elem#1#0, Ideal) then (
                elem1 = for x in elem#1 list
                    select(x_*, f -> degree f === {1});
                if #elem1#0 === 0 then continue;
                elem2 := for x in elem#2 list
                    select(x_*, f -> degree f === {1});
                {elem#0, elem1/ideal, elem2/ideal}
                )
            else elem
        ))
    )
selectLinear Nothing := nothing -> null

invertibleMatrixOverZZ = method()
invertibleMatrixOverZZ(Matrix, Ideal) := Sequence => (A, J) -> (
    -- returns (determinacy, A0), or (INCONSISTENT, null) or ...?
    if J == 1 then 
        (INCONSISTENT, null)
    else (
        A0 := A % J;
        detA0 := (det A0) % J;
        suppA0 := support A0;
        if suppA0 === {} then (
            -- In this case we either have an integer matrix, or a rational matrix.
            if liftable(detA0, ZZ) and all(flatten entries A0, f -> liftable(f, ZZ)) then
                return (CONSISTENT, sub(A0, ZZ));
            return (INCONSISTENT, sub(A0, QQ));
            );
        jc := decompose J;
        if isPrime J then return (INDETERMINATE, jc);
        possibles := for j in jc list (
            ans := invertibleMatrixOverZZ(A0, j);
            if ans#0 == CONSISTENT then return ans else ans
            );
        if all(possibles, a -> a#0 === INCONSISTENT) then (
            ans := select(1, possibles, a -> instance(a#1, Matrix));
            if #ans > 0 then return ans#0 else return (INCONSISTENT, null);
            )
        else
            return (INDETERMINATE, jc);
        )
    )

equivalenceIdeal = method()
equivalenceIdeal(List, List, Ring, Sequence) := Ideal => (List1, List2, RQ, Aphi) -> (
    if itemType List1 =!= itemType List2 then 
        error("expected two lists to have the same list of types, they are: " 
            | toString itemType List1 | " and " | toString itemType List2);
    (A,phi) := Aphi;
    U := target phi;
    if U =!= source phi then error "expected ring map with same source and target";
    if #List1 == 0 then return ideal(0_U);
    B := coefficientRing U;
    toU := map(U, RQ, vars U);
    toB := map(B, U);
    List1 = List1/(I -> if ring I =!= RQ then toU (sub(I, RQ)) else toU I);
    List2 = List2/(I -> if ring I =!= RQ then toU (sub(I, RQ)) else toU I);
    -- List1 and List2 are lists with the same length, consisting of RingElement's, Ideal's, Matrices.
    -- List1 and List2 should each have RingElement's and Ideal's in the same spot.
    ids := for i from 0 to #List1-1 list (
        if instance(List1#i, RingElement) then (
            ideal toB (last coefficients(phi List1#i - List2#i))
            )
        else if instance(List1#i, Ideal) then (
            ideal toB (last coefficients((gens phi List1#i) % List2#i))
            )
        else if instance(List1#i, Matrix) then (
            rowvec := (numrows List1#i === 1);
            -- if rowvec is false, then this must be a column vector.
            if rowvec then
                ideal toB last coefficients sub(List2#i * (transpose A) - List1#i, U)
            else
                ideal toB last coefficients sub((transpose A) * List1#i - List2#i, U)
            )
        );
    sum ids
    )

equivalenceIdeals = method()
equivalenceIdeals(MatchingData, Ring, Sequence) := List => (MD, RQ, Aphi) -> (
    -- returns the list of equivalence ideals.
    (A,phi) := Aphi;
    (src, tar) := matches MD;
    for i from 0 to #tar-1 list (
        trim equivalenceIdeal(src, tar#i, RQ, Aphi)
        )
    )
equivalenceIdeals(Nothing, Ring, Sequence) := (MDnull, RQ, Aphi) -> {}

tryEquivalences = method()
tryEquivalences(MatchingData, Ring, Sequence) := (MD, RQ, Aphi) -> (
    (A,phi) := Aphi;
    badJs := {};
    inconsistentMatrix := null;
    (src, tar) := matches MD;
    for i from 0 to #tar-1 do (
        << "doing " << i << endl;
        J := trim equivalenceIdeal(src, tar#i, RQ, Aphi);
        ans := invertibleMatrixOverZZ(A, J);
        if ans#0 == CONSISTENT then return ans;
        if ans#0 == INCONSISTENT and instance(ans#1, Matrix) then inconsistentMatrix = ans#1;
        if ans#0 == INDETERMINATE then (
            badJs = join(badJs, ans#1);
            );
        );
    if #badJs > 0 then return (INDETERMINATE, badJs);
    (INCONSISTENT, inconsistentMatrix)
    )
tryEquivalences(Nothing, Ring, Sequence) := (MD, RQ, Aphi) -> (INCONSISTENT, null)

-------------------------------
-- Finding matching data of (L1,F1), (L2,F2)
-- using singular loci and hessian factorizations
-------------------------------
factors = method()
factors RingElement := (F) -> (
     facs := factor F;
     facs//toList/toList/reverse
     )

factorsByType = method()
factorsByType RingElement := HashTable => F -> (
    facs := factors F;
    faclist := for fx in facs list (fx#0, sum first exponents fx#1, fx#1);
    H := partition(x -> {x#0, x#1}, faclist);
    hashTable for k in keys H list k => for x in H#k list x_2
    )   

factorType = method()
factorType RingElement := HashTable => F -> (
    fac1 := factorsByType F;
    keys1 := sort select(keys fac1, k -> k =!= {1,0});
    hashTable for k in keys1 list k => #fac1#k
    )

idealsByBetti = method()
idealsByBetti(List, List) := MatchingData => (J1s, J2s) -> (
    H1 := partition(J -> betti gens J, J1s);
    H2 := partition(J -> betti gens J, J2s);
    if sort keys H1 =!= sort keys H2 then return null;
    matchingData for k in sort keys H1 list (
        if #H1#k === 1 then H1#k#0 => H2#k#0 else {Permutations, H1#k, H2#k}
        )
    )

-- TODO: place in M2 Core?
hessian = method()
hessian RingElement := F -> diff(vars ring F, diff(transpose vars ring F, F))

hessianMatches = method()
hessianMatches(RingElement, RingElement) := MatchingData => (F1, F2) -> (
    fac1 := factorsByType(det hessian F1);
    fac2 := factorsByType(det hessian F2);
    keys1 := sort select(keys fac1, k -> k =!= {1,0}); -- remove constant
    keys2 := sort select(keys fac2, k -> k =!= {1,0}); -- remove constant
    if keys1 =!= keys2 then return null; -- no matches.
    matchingData for k in keys1 list {SignedPermutations, fac1#k, fac2#k}
    )

-- Note: the returned value over a finite field is fine, and over QQ, all denominators and lcms have been cleared
-- TODO: maybe this should be used if the ideal is generated by linears, and of codim = number of vars - 1?
-- IE: modify a MatchingData to change this.
singularPoints = method()
singularPoints RingElement := List => F -> (
    if not isHomogeneous F then error "expected homogeneous polynomial";
    R := ring F;
    kk := coefficientRing R;
    n := numgens R;
    singlocus := trim saturate(ideal F + ideal jacobian F);
    if singlocus == 1 then return {};
    comps := (decompose singlocus);
    comps0 := select(comps, c -> codim c == n-1 and degree c === 1); -- zero-dimensional rational points
    comps1 := select(comps, c -> not(codim c == n-1 and degree c === 1)); -- the rest
    comps1 = comps1/trim;
    if #comps1 != 0 then 
        << "other singular components: " << netList comps1 << endl;
    for c in comps0 list (
        pt := (vars R) % c;
        vs := support pt;
        if #vs != 1 then error "internal error: number of variables is not 1";
        pt = sub(pt, vs#0 => 1_kk);
        pt = flatten entries pt;
        if kk === QQ then (
            g := gcd pt;
            pt = 1/g * pt
            );
        transpose matrix{pt}
        )
    )
singularPointMatches = method()
singularPointMatches(RingElement, RingElement) := MatchingData => (F1, F2) -> (
    pts1 := singularPoints F1;
    pts2 := singularPoints F2;
    if #pts1 =!= #pts2 then return null; -- no matches.
    matchingData{{SignedPermutations, pts1, pts2}}
    )

-- NOT DONE YET!!
singularMatches = method()
singularMatches(RingElement, RingElement) := MatchingData => (F1, F2) -> (
    sing1 := trim saturate(ideal F1 + ideal jacobian F1);
    sing2 := trim saturate(ideal F2 + ideal jacobian F2);
    if sing1 == 1 then return matchingData{};
    comps1 := (decompose sing1)/trim;
    comps2 := (decompose sing2)/trim;
    idealsByBetti(comps1, comps2)
    )

findEquivalence = method()
findEquivalence(List, List) := (LF1, LF2) -> (
    (L1, F1) := toSequence LF1;
    (L2, F2) := toSequence LF2;
    R := ring L1;
    -- TODO: check that R is the ring of all 4 of these.
    -- TODO: check that coefficient ring is ZZ, QQ, finite field, or what else is allowed?
    RQ := R;
    toRQ := identity;
    if coefficientRing R === ZZ then (
        RQ = QQ (monoid R); -- change ZZ to QQ, leave finite fields alone.
        toRQ = map(RQ, R, vars RQ);
        );
    L1 = toRQ L1;
    L2 = toRQ L2;
    F1 = toRQ F1;
    F2 = toRQ F2;
    (A, phi) := genericLinearMap RQ;
    md := hessianMatches(F1, F2) |
          singularMatches(F1, F2) |
          matchingData {L1 => L2, F1 => F2};
    if false then (
        H1 := det hessian F1;
        H2 := det hessian F2;
        singH1 := ideal H1 + ideal jacobian H1;
        singH2 := ideal H2 + ideal jacobian H2;
        comps1 := (decompose singH1)/trim;
        comps2 := (decompose singH2)/trim;
        MDh := idealsByBetti(comps1, comps2);
        md = md | MDh;
        );
    if false then (
        linmd := (selectLinear md) | matchingData{F1 => F2};
        result := tryEquivalences(linmd, RQ, (A,phi));
        return result;
        );
    -- if result is INDETERMINATE, try the entire matching data
    -- TODO: if we get a consistent match, try that first!
    -- only if that fails should we move on to this.
    tryEquivalences(md, RQ, (A, phi))
    -- if result#0 =!= INCONSISTENT then (
    --     result2 := tryEquivalences(md, RQ, (A,phi));
    --     (result, result2)
    --     )
    -- else result
    )

findEquivalenceHessianSingularities = method()
findEquivalenceHessianSingularities(List, List) := (LF1, LF2) -> (
    (L1, F1) := toSequence LF1;
    (L2, F2) := toSequence LF2;
    R := ring L1;
    -- TODO: check that R is the ring of all 4 of these.
    -- TODO: check that coefficient ring is ZZ, QQ, finite field, or what else is allowed?
    RQ := R;
    toRQ := identity;
    if coefficientRing R === ZZ then (
        RQ = QQ (monoid R); -- change ZZ to QQ, leave finite fields alone.
        toRQ = map(RQ, R, vars RQ);
        );
    L1 = toRQ L1;
    L2 = toRQ L2;
    F1 = toRQ F1;
    F2 = toRQ F2;
    (A, phi) := genericLinearMap RQ;
    md := hessianMatches(F1, F2);
    if true then (
        H1 := det hessian F1;
        H2 := det hessian F2;
        singH1 := ideal H1 + ideal jacobian H1;
        singH2 := ideal H2 + ideal jacobian H2;
        comps1 := (decompose singH1)/trim;
        comps2 := (decompose singH2)/trim;
        MDh := idealsByBetti(comps1, comps2);
        md = md | MDh;
        );
    md = md | matchingData {L1 => L2};
    if true then (
        linmd := (selectLinear md) | matchingData{F1 => F2};
        result := tryEquivalences(linmd, RQ, (A,phi));
        return result;
        );
    )

beginDocumentation()

doc ///
  Key
    IntegerEquivalences
  Headline
    finding invertible integral matrices preserving points, linear forms and ideals
  Description
    Text
      This package provides tools for determining whether two sets of polynomials,
      linear forms, ideals, and lattice points are related by an invertible integer
      change of coordinates (i.e., an element of $GL(n, \ZZ)$).

      The main use case is determining whether two Calabi-Yau 3-folds have equivalent
      topological data (cubic forms, second Chern classes) up to a change of basis
      in $H^2(X, \ZZ)$.

      @SUBSECTION "Overview of the approach"@
    Text
      Given source data (polynomials $F_1$, linear forms $L_1$, ideals, lattice points)
      and target data $(F_2, L_2, \ldots)$, the package:

      (1) Sets up a generic $n \times n$ matrix $A$ of unknowns via @TO genericLinearMap@.

      (2) Creates @TO MatchingData@ specifying which source items map to which target items,
      and how to enumerate possible matchings (fixed, permutations, or signed permutations).

      (3) For each matching, computes the @TO equivalenceIdeal@ — the ideal of polynomial
      constraints on the entries of $A$.

      (4) Checks via @TO invertibleMatrixOverZZ@ whether the constraints are satisfied by
      an integer matrix with determinant $\pm 1$.

      The function @TO findEquivalence@ automates this pipeline for the common case of
      matching a linear form and a cubic form.
    Text
      @SUBSECTION "A simple example"@
    Text
      We check whether two pairs $(L_1, F_1)$ and $(L_2, F_2)$ of a linear form
      and cubic form (representing $c_2$ and cubic intersection form of CY 3-folds)
      are related by a $GL(3, \ZZ)$ change of coordinates.
    Example
      RZ = ZZ[a,b,c]
      L1 = 8*a-4*b+36*c
      F1 = 2*a^3-3*a^2*b-3*a*b^2+8*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*a*c^2
      L2 = -4*a+8*b+36*c
      F2 = 8*a^3-3*a^2*b-3*a*b^2+2*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*b*c^2
      result = findEquivalence({L1, F1}, {L2, F2})
      result#0
      result#1
    Text
      @SUBSECTION "A sample use of the pipeline"@
    Text
      We check whether two pairs $(L_1, F_1)$ and $(L_2, F_2)$ of a linear form
      and cubic form (representing $c_2$ and cubic intersection form of CY 3-folds)
      are related by a $GL(3, \ZZ)$ change of coordinates.  However, this time
      the two cubic forms (in 4 variables), using the default pipeline doesn't work well.
    Example
      RZ = ZZ[a,b,c,d]
      RQ = QQ[a,b,c,d]
      L1 = 2*a+26*b+8*c+18*d
      F1 = 5*a^3-9*a^2*b+3*a*b^2-b^3-6*a^2*c+12*a*b*c-6*b*c^2+2*c^3+3*b^2*d+
            12*b*c*d-6*c^2*d-3*b*d^2+6*c*d^2-3*d^3
      L2 = 2*a+8*b+32*c+26*d
      F2 = 5*a^3-6*a^2*b+2*b^3-9*a^2*c+12*a*b*c-12*b^2*c+3*a*c^2+18*b*c^2-
          10*c^3-9*a^2*d+12*a*b*d-6*b^2*d+6*a*c*d+
          12*b*c*d-6*c^2*d+3*a*d^2-d^3
    Text
      Currently, the following call would not terminate quickly:
    Pre
      result = findEquivalence({L1, F1}, {L2, F2})
    Text
      Instead we run through the pipeline by hand, using
    Example
      (A, phi) = genericLinearMap RQ
      H1 = det hessian F1
      H2 = det hessian F2
      sing1 = trim saturate(ideal H1 + ideal jacobian H1)
      sing2 = trim saturate(ideal H2 + ideal jacobian H2)
      comps1 = (decompose sing1)/trim;
      comps2 = (decompose sing2)/trim;
      MD = idealsByBetti(comps1, comps2)
      MD1 = (selectLinear MD) | matchingData{L1 => L2, F1 => F2}
      MDall = MD | matchingData{L1 => L2, F1 => F2}
      tryEquivalences(MD1, RQ, (A,phi))
      tryEquivalences(MDall, RQ, (A,phi))
      MD = (selectLinear singularMatches(H1, H2)) | matchingData{L1 => L2, F1 => F2}
      tryEquivalences(MD, RQ, (A,phi))
      Js = equivalenceIdeals(MD, RQ, (A,phi))

      (A, phi) = genericLinearMap RQ
      H1 = det hessian F1
      H2 = det hessian F2
      MD = (selectLinear singularMatches(H1, H2)) | matchingData{L1 => L2, F1 => F2}
      tryEquivalences(MD, RQ, (A,phi))
      Js = equivalenceIdeals(MD, RQ, (A,phi))
  SeeAlso
    findEquivalence
    MatchingData
    equivalenceIdeal
    genericLinearMap
    invertibleMatrixOverZZ
///

doc ///
  Key
    findEquivalence
    (findEquivalence, List, List)
  Headline
    find an integer change of basis relating two pairs of forms
  Usage
    result = findEquivalence({L1, F1}, {L2, F2})
  Inputs
    :List
      of the form {\{L1, F1\}} where $L_1$ is a linear form and $F_1$ is a polynomial
    :List
      of the form {\{L2, F2\}} where $L_2$ is a linear form and $F_2$ is a polynomial
  Outputs
    result:Sequence
      a pair {\tt (status, A)} where status is @TO CONSISTENT@, @TO INCONSISTENT@,
      or @TO INDETERMINATE@, and $A$ is the change-of-basis matrix (if consistent)
  Description
    Text
      This is the main high-level function.  Given two pairs of forms $(L_1, F_1)$
      and $(L_2, F_2)$ in a polynomial ring, it searches for an invertible integer
      matrix $A$ such that the corresponding change of variables sends $L_1 \mapsto L_2$
      and $F_1 \mapsto F_2$.

      The function automatically uses hessian and singular locus data to constrain
      the search.
    Example
      RZ = ZZ[a,b,c]
      L1 = 8*a-4*b+36*c
      F1 = 2*a^3-3*a^2*b-3*a*b^2+8*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*a*c^2
      L2 = -4*a+8*b+36*c
      F2 = 8*a^3-3*a^2*b-3*a*b^2+2*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*b*c^2
      (stat, A) = findEquivalence({L1, F1}, {L2, F2})
      stat == CONSISTENT
      A
  SeeAlso
    MatchingData
    tryEquivalences
    equivalenceIdeal
    CONSISTENT
    INCONSISTENT
    INDETERMINATE
///

doc ///
  Key
    MatchingData
  Headline
    type specifying how source and target data should be matched
  Description
    Text
      A {\tt MatchingData} object is a list of matching specifications.  Each
      element is either:

      (1) {\tt source => target} — a fixed matching (ring element to ring element,
      ideal to ideal, or row/column matrix to row/column matrix), or

      (2) {\tt \{type, sourceList, targetList\}} — where {\tt type} is @TO Permutations@
      or @TO SignedPermutations@, specifying that all permutations (or signed permutations)
      of the target list should be tried.

      Use @TO matchingData@ to construct validated instances.
    Example
      R = QQ[a,b,c]
      L1 = 10*a + 28*b + 26*c
      L2 = 16*a + 10*b + 26*c
      F1 = a^3 - 3*a^2*b
      F2 = -2*a^3 - 3*a^2*b
      md = matchingData {L1 => L2, F1 => F2}
  SeeAlso
    matchingData
    matches
    Permutations
    SignedPermutations
///

doc ///
  Key
    matchingData
    (matchingData, List)
  Headline
    create validated matching data
  Usage
    md = matchingData L
  Inputs
    L:List
      of matching specifications (see @TO MatchingData@)
  Outputs
    md:MatchingData
  Description
    Text
      Creates a @TO MatchingData@ object from a list of matching specifications.
      Each element of the list should be either {\tt source => target} (for fixed matchings)
      or {\tt \{type, sourceList, targetList\}} where type is @TO Permutations@ or
      @TO SignedPermutations@.

      The constructor validates that sources and targets have matching types and lengths.
    Example
      R = ZZ[a,b,c]
      md = matchingData {
          a^2+b => a^2-b,
          {SignedPermutations, {a, b, c}, {a+b, a+c, 2*b+c}},
          matrix{{1,2,3}} => matrix{{1,-2,1}}
          }
    Text
      Two MatchingData objects can be combined using @TT "|"@.
    Example
      md1 = matchingData {a => b}
      md2 = matchingData {ideal(a,b) => ideal(b,c)}
      md1 | md2
  SeeAlso
    MatchingData
    matches
    Permutations
    SignedPermutations
///

doc ///
  Key
    matches
    (matches, MatchingData)
  Headline
    enumerate all source-target matchings from matching data
  Usage
    (src, targs) = matches md
  Inputs
    md:MatchingData
  Outputs
    src:List
      the flattened source data
    targs:List
      a list of all possible flattened target lists
  Description
    Text
      Given @TO MatchingData@, this function computes the source list (which is fixed)
      and all possible target lists obtained by applying the specified permutations
      and signed permutations.

      For fixed matchings ({\tt source => target}), the target is always the same.
      For {\tt Permutations} entries, all permutations of the target list are generated.
      For {\tt SignedPermutations} entries, all signed permutations are generated.
      The Cartesian product of all these choices gives the full list of targets.
    Example
      R = ZZ[a,b,c]
      md = matchingData {
          a => b,
          {Permutations, {a, b}, {a+b, a+c}}
          }
      (src, targs) = matches md
      src
      #targs
  SeeAlso
    MatchingData
    equivalenceIdeal
///

doc ///
  Key
    genericLinearMap
    (genericLinearMap, Ring)
    [genericLinearMap, Variable]
  Headline
    create a generic linear change of coordinates
  Usage
    (A, phi) = genericLinearMap R
    (A, phi) = genericLinearMap(R, Variable => t)
  Inputs
    R:Ring
      a polynomial ring in $n$ variables
    Variable => Symbol
      the variable name to use for entries of $A$ (default: {\tt t})
  Outputs
    A:Matrix
      an $n \times n$ generic matrix over a new ring $T = K[t_{1,1}, \ldots, t_{n,n}]$
    phi:RingMap
      a ring map $U \to U$ where $U = T[x_1,\ldots,x_n]$ sending each variable
      to the corresponding linear combination given by $A$
  Description
    Text
      Creates a generic $n \times n$ matrix $A$ of new indeterminates and the corresponding
      ring map $\phi$ that acts on polynomials by the linear change of variables defined by $A$.
      This is the setup step for computing @TO equivalenceIdeal@.
    Example
      R = QQ[a,b,c]
      (A, phi) = genericLinearMap R
      A
      U = target phi
      phi(U_0)
    Text
      The {\tt Variable} option allows choosing the variable name.
    Example
      (A, phi) = genericLinearMap(R, Variable => symbol s)
      A
  SeeAlso
    equivalenceIdeal
///

doc ///
  Key
    equivalenceIdeal
    (equivalenceIdeal, List, List, Ring, Sequence)
  Headline
    compute the ideal of constraints for a linear map to match source to target
  Usage
    J = equivalenceIdeal(src, tar, RQ, (A, phi))
  Inputs
    src:List
      source data (ring elements, ideals, row/column matrices)
    tar:List
      target data (same types as source)
    RQ:Ring
      the polynomial ring (typically over QQ)
    :Sequence
      the pair {\tt (A, phi)} from @TO genericLinearMap@
  Outputs
    J:Ideal
      an ideal in the entries of $A$ whose solutions give change-of-basis matrices
  Description
    Text
      For each source-target pair, this function imposes constraints:

      $\bullet$ Ring elements: $\phi(F_1) = F_2$ (coefficient matching).

      $\bullet$ Ideals: $\phi(I_1) \subseteq I_2$ (generators of $\phi(I_1)$ reduce to zero modulo $I_2$).

      $\bullet$ Row vectors: $v_2 A^T = v_1$ (row vectors transform contravariantly).

      $\bullet$ Column vectors: $A^T v_1 = v_2$ (column vectors transform covariantly).
    Example
      RQ = QQ[a,b,c]
      (A, phi) = genericLinearMap RQ
      L1 = 10*a + 28*b + 26*c
      L2 = 16*a + 10*b + 26*c
      J = equivalenceIdeal({L1}, {L2}, RQ, (A, phi))
  SeeAlso
    genericLinearMap
    invertibleMatrixOverZZ
///

doc ///
  Key
    invertibleMatrixOverZZ
    (invertibleMatrixOverZZ, Matrix, Ideal)
  Headline
    check whether an equivalence ideal has a solution over the integers
  Usage
    (status, result) = invertibleMatrixOverZZ(A, J)
  Inputs
    A:Matrix
      the generic matrix from @TO genericLinearMap@
    J:Ideal
      the equivalence ideal from @TO equivalenceIdeal@
  Outputs
    status:Symbol
      one of @TO CONSISTENT@, @TO INCONSISTENT@, or @TO INDETERMINATE@
    result:Thing
      a @TO Matrix@ over $\ZZ$ if consistent, or {\tt null}/other info if not
  Description
    Text
      Reduces the generic matrix $A$ modulo the ideal $J$ and checks whether
      the result is an integer matrix with determinant $\pm 1$.

      If $J = (1)$, the system is inconsistent (no solutions at all).
      If the reduced matrix is over $\ZZ$ with $\det = \pm 1$, it is consistent.
      Otherwise, the function decomposes $J$ and checks each component.
    Example
      RQ = QQ[a,b,c]
      (A, phi) = genericLinearMap RQ
      L1 = 8*a-4*b+36*c
      F1 = 2*a^3-3*a^2*b-3*a*b^2+8*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*a*c^2
      L2 = -4*a+8*b+36*c
      F2 = 8*a^3-3*a^2*b-3*a*b^2+2*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*b*c^2
      J = equivalenceIdeal({L1, F1}, {L2, F2}, RQ, (A, phi))
      invertibleMatrixOverZZ(A, J)
  SeeAlso
    equivalenceIdeal
    CONSISTENT
    INCONSISTENT
    INDETERMINATE
///

doc ///
  Key
    tryEquivalences
    (tryEquivalences, MatchingData, Ring, Sequence)
  Headline
    try all matchings and search for an integer equivalence
  Usage
    result = tryEquivalences(md, RQ, (A, phi))
  Inputs
    md:MatchingData
    RQ:Ring
      the polynomial ring (typically over QQ)
    :Sequence
      the pair {\tt (A, phi)} from @TO genericLinearMap@
  Outputs
    result:Sequence
      a pair {\tt (status, matrix or info)}
  Description
    Text
      Iterates over all matchings generated by @TO matches@ from the @TO MatchingData@,
      computes the @TO equivalenceIdeal@ for each, and checks for integer solutions via
      @TO invertibleMatrixOverZZ@.  Returns as soon as a @TO CONSISTENT@ solution is found.
    Example
      RQ = QQ[a,b,c]
      (A, phi) = genericLinearMap RQ
      L1 = 10*a + 28*b + 26*c
      L2 = 16*a + 10*b + 26*c
      F1 = a^3-3*a^2*b-3*a*b^2-2*b^3-3*a^2*c+6*a*b*c+6*b^2*c+3*a*c^2+6*b*c^2-c^3
      F2 = -2*a^3-3*a^2*b-3*a*b^2+b^3+6*a*b*c-3*b^2*c+6*a*c^2+3*b*c^2-c^3
      md = matchingData{L1 => L2, F1 => F2}
      result = tryEquivalences(md, RQ, (A, phi))
      result#0
  SeeAlso
    MatchingData
    matches
    equivalenceIdeal
    invertibleMatrixOverZZ
///

doc ///
  Key
    CONSISTENT
  Headline
    result code indicating an integer equivalence was found
  Description
    Text
      Returned by @TO invertibleMatrixOverZZ@ and @TO tryEquivalences@ when an
      integer matrix with determinant $\pm 1$ satisfying all constraints has been found.
  SeeAlso
    INCONSISTENT
    INDETERMINATE
    invertibleMatrixOverZZ
///

doc ///
  Key
    INCONSISTENT
  Headline
    result code indicating no equivalence exists
  Description
    Text
      Returned by @TO invertibleMatrixOverZZ@ and @TO tryEquivalences@ when the
      constraints have no solution, or no solution that gives an integer matrix
      with determinant $\pm 1$.
  SeeAlso
    CONSISTENT
    INDETERMINATE
    invertibleMatrixOverZZ
///

doc ///
  Key
    INDETERMINATE
  Headline
    result code indicating the equivalence could not be determined
  Description
    Text
      Returned by @TO invertibleMatrixOverZZ@ and @TO tryEquivalences@ when the
      system could not be fully resolved — e.g., a prime ideal component
      could not be solved uniquely.
  SeeAlso
    CONSISTENT
    INCONSISTENT
    invertibleMatrixOverZZ
///

doc ///
  Key
    Permutations
  Headline
    matching type: try all permutations
  Description
    Text
      Used in @TO MatchingData@ to indicate that all permutations of the target list
      should be tried.  Compare with @TO SignedPermutations@.
    Example
      R = ZZ[a,b,c]
      md = matchingData {{Permutations, {a, b}, {a+b, a+c}}}
      (src, targs) = matches md
      #targs == 2
  SeeAlso
    SignedPermutations
    MatchingData
///

doc ///
  Key
    SignedPermutations
  Headline
    matching type: try all signed permutations
  Description
    Text
      Used in @TO MatchingData@ to indicate that all signed permutations (permutations
      combined with sign changes) of the target list should be tried.  For a list of
      length $k$, this produces $k! \cdot 2^k$ matchings.
    Example
      R = ZZ[a,b,c]
      md = matchingData {{SignedPermutations, {a, b}, {a+b, a+c}}}
      (src, targs) = matches md
      #targs == 8  -- 2! * 2^2
  SeeAlso
    Permutations
    MatchingData
///

doc ///
  Key
    RowVector
  Headline
    type indicator for row matrices in matching data
  Description
    Text
      Used internally to classify a $1 \times n$ matrix in @TO MatchingData@ as a row vector.
      Row vectors transform as $v_2 A^T = v_1$.
  SeeAlso
    ColumnVector
    MatchingData
///

doc ///
  Key
    ColumnVector
  Headline
    type indicator for column matrices in matching data
  Description
    Text
      Used internally to classify an $n \times 1$ matrix in @TO MatchingData@ as a column vector.
      Column vectors transform as $A^T v_1 = v_2$.
  SeeAlso
    RowVector
    MatchingData
///

doc ///
  Key
    Unknown
  Headline
    type indicator for unrecognized items in matching data
  Description
    Text
      Returned by internal type classification when an item in @TO MatchingData@
      is not a recognized type (ring element, ideal, row vector, or column vector).
      This causes validation to fail.
  SeeAlso
    MatchingData
    RowVector
    ColumnVector
///

doc ///
  Key
    extendToMatrix
    (extendToMatrix, List)
  Headline
    extend an integer vector to an invertible integer matrix
  Usage
    A = extendToMatrix v
  Inputs
    v:List
      a list of integers
  Outputs
    A:Matrix
      an invertible integer matrix whose last row times the column vector $v$
      gives $(0, \ldots, 0, \gcd(v))$
  Description
    Text
      Given a list of integers, uses the LLL algorithm to find an invertible
      integer matrix $A$ such that $A \cdot v^T = (0, \ldots, 0, g)^T$ where
      $g = \gcd(v)$.
    Example
      A = extendToMatrix{10,15,6}
      A * transpose matrix{{10,15,6}}
      det A
    Example
      A = extendToMatrix{10,15,1}
      A * transpose matrix{{10,15,1}}
      det A
  SeeAlso
    invertibleMatrixOverZZ
///

doc ///
  Key
    hessian
    (hessian, RingElement)
  Headline
    compute the Hessian matrix of a polynomial
  Usage
    H = hessian F
  Inputs
    F:RingElement
  Outputs
    H:Matrix
      the matrix of second partial derivatives of $F$
  Description
    Text
      Computes the Hessian matrix $(\partial^2 F / \partial x_i \partial x_j)$.
    Example
      R = QQ[a,b,c]
      F = a^3 + b^3 + c^3
      hessian F
///

doc ///
  Key
    hessianMatches
    (hessianMatches, RingElement, RingElement)
  Headline
    create matching data from Hessian determinant factorizations
  Usage
    md = hessianMatches(F1, F2)
  Inputs
    F1:RingElement
    F2:RingElement
  Outputs
    md:MatchingData
  Description
    Text
      Computes the determinant of the Hessian of each polynomial, factors it,
      groups factors by type (multiplicity and degree), and creates @TO MatchingData@
      with @TO SignedPermutations@ matching for each group.

      This constrains the search space for @TO tryEquivalences@ by requiring that
      factors of the Hessian determinant are mapped to corresponding factors.
    Example
      R = QQ[a,b,c]
      F1 = 2*a^3-3*a^2*b-3*a*b^2+8*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*a*c^2
      F2 = 8*a^3-3*a^2*b-3*a*b^2+2*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*b*c^2
      md = hessianMatches(F1, F2)
  SeeAlso
    hessian
    factorsByType
    MatchingData
///

doc ///
  Key
    factorsByType
    (factorsByType, RingElement)
  Headline
    factor a polynomial and group factors by multiplicity and degree
  Usage
    H = factorsByType F
  Inputs
    F:RingElement
  Outputs
    H:HashTable
      keys are {\tt \{multiplicity, degree\}}, values are lists of irreducible factors
  Description
    Text
      Factors a polynomial and groups the irreducible factors by their multiplicity
      and total degree.  Used by @TO hessianMatches@ to match factors of Hessian
      determinants between two polynomials.
    Example
      R = QQ[a,b,c]
      F = a^2*b*(a+b)^2*(a-c)
      factorsByType F
  SeeAlso
    hessianMatches
///

doc ///
  Key
    singularPoints
    (singularPoints, RingElement)
  Headline
    find rational singular points of a projective hypersurface
  Usage
    pts = singularPoints F
  Inputs
    F:RingElement
      a homogeneous polynomial
  Outputs
    pts:List
      a list of column vectors representing projective singular points
  Description
    Text
      Computes the singular locus of the projective hypersurface $V(F)$ and
      returns the rational singular points as column vectors (with integer
      entries if the coefficient ring is $\QQ$).
    Example
      R = QQ[a,b,c]
      F = a^2*b - a*b^2
      singularPoints F
  SeeAlso
    singularPointMatches
    singularMatches
///

doc ///
  Key
    singularPointMatches
    (singularPointMatches, RingElement, RingElement)
  Headline
    create matching data from singular points of two hypersurfaces
  Usage
    md = singularPointMatches(F1, F2)
  Inputs
    F1:RingElement
    F2:RingElement
  Outputs
    md:MatchingData
      or @TO null@ if the number of singular points differs
  Description
    Text
      Computes the singular points of the projective hypersurfaces $V(F_1)$ and $V(F_2)$
      and creates @TO MatchingData@ with @TO SignedPermutations@ matching.
      Returns @TO null@ if the two hypersurfaces have different numbers of singular points.
  SeeAlso
    singularPoints
    singularMatches
    MatchingData
///

doc ///
  Key
    singularMatches
    (singularMatches, RingElement, RingElement)
  Headline
    create matching data from singular locus components of two hypersurfaces
  Usage
    md = singularMatches(F1, F2)
  Inputs
    F1:RingElement
    F2:RingElement
  Outputs
    md:MatchingData
  Description
    Text
      Computes the decomposition of the singular loci of $V(F_1)$ and $V(F_2)$
      and matches components using @TO idealsByBetti@ (grouping by Betti numbers).
  SeeAlso
    singularPoints
    singularPointMatches
    idealsByBetti
///

doc ///
  Key
    idealsByBetti
    (idealsByBetti, List, List)
  Headline
    create matching data from two lists of ideals grouped by Betti numbers
  Usage
    md = idealsByBetti(J1s, J2s)
  Inputs
    J1s:List
      a list of ideals
    J2s:List
      a list of ideals
  Outputs
    md:MatchingData
  Description
    Text
      Groups ideals by their Betti numbers.  Ideals in groups of size 1
      are matched directly; groups with more than one ideal are matched
      using @TO Permutations@.
  SeeAlso
    singularMatches
    MatchingData
///

doc ///
  Key
    selectLinear
    (selectLinear, MatchingData)
  Headline
    extract linear constraints from matching data
  Usage
    md1 = selectLinear md
  Inputs
    md:MatchingData
  Outputs
    md1:MatchingData
      matching data restricted to linear forms and row/column matrices
  Description
    Text
      Filters @TO MatchingData@ to keep only the linear ring elements, the linear
      generators of ideals, and all row/column matrices.  This can be used to
      first try a faster solve using only linear constraints.
  SeeAlso
    MatchingData
    tryEquivalences
///

doc ///
  Key
    allSigns
    (allSigns, List)
  Headline
    generate all sign combinations of a list
  Usage
    result = allSigns L
  Inputs
    L:List
  Outputs
    result:List
      a list of lists, each with the same elements as $L$ but with
      all possible sign combinations
  Description
    Example
      allSigns {1,2}
  SeeAlso
    signedPermutations
///

doc ///
  Key
    signedPermutations
    (signedPermutations, List)
  Headline
    generate all signed permutations of a list
  Usage
    result = signedPermutations L
  Inputs
    L:List
  Outputs
    result:List
      all permutations of $L$ combined with all sign changes
  Description
    Text
      For a list of length $k$, produces $k! \cdot 2^k$ signed permutations.
    Example
      signedPermutations {1,2}
      #signedPermutations {1,2,3}
  SeeAlso
    allSigns
    SignedPermutations
///

doc ///
  Key
    cartesian
    (cartesian, List)
  Headline
    Cartesian product of a list of lists
  Usage
    result = cartesian Ls
  Inputs
    Ls:List
      a list of lists
  Outputs
    result:List
      the Cartesian product, as a list of lists
  Description
    Example
      cartesian {{1,2}, {3,4}}
      cartesian {{1,2}, {3,4}, {5,6}}
      cartesian {{1,2}, {3}, {5,6}}
  SeeAlso
    matches
///

TEST ///
-- These 3 forms were generated from h11=5 database, with:
-- {(2249, 0), (2255, 0), (2270, 0)}
-- L1 = c2Form Xs#(2249,0)
-- L2 = c2Form Xs#(2255,0)
-- L3 = c2Form Xs#(2270,0)
-- F1 = cubicForm Xs#(2249,0)
-- F2 = cubicForm Xs#(2255,0)
-- F3 = cubicForm Xs#(2270,0)

-*
  restart
  needsPackage "IntegerEquivalences"
*-
-- This is supposed to test formation of matchingData
  debug needsPackage "IntegerEquivalences"
  RZ = ZZ[a,b,c,d,e]
  RQ = QQ (monoid RZ)
  (A,phi) = genericLinearMap RQ
  use RQ
  L1 = 12*a+10*b+34*c-4*d+4*e
  L2 = -4*a+12*b+18*c+28*d+4*e
  L3 = -4*a+12*b+26*c+10*d+20*e
  F1 = b^3-6*a^2*c-3*b^2*c+3*b*c^2+c^3+6*a^2*d+12*a*c*d+6*c^2*d-12*a*d^2-
      12*c*d^2+8*d^3-6*a^2*e+12*a*c*e+12*a*e^2-12*c*e^2-8*e^3
  F2 = 8*a^3-12*a^2*b+6*a*b^2-12*a^2*c+12*a*b*c-6*b^2*c+6*a*c^2-3*c^3-
      12*a^2*d+12*a*b*d-6*b^2*d+12*a*c*d-6*c^2*d+6*a*d^2-6*c*d^2-
      2*d^3+12*c*d*e+6*d^2*e-12*c*e^2-8*e^3
  F3 = 8*a^3-12*a^2*b+6*a*b^2-12*a^2*c+12*a*b*c-6*b^2*c+6*a*c^2-c^3+
      3*c^2*d-3*c*d^2+d^3-6*c*e^2-4*e^3

  findEquivalence({L1, F1}, {L2, F2})
  findEquivalence({L3, F3}, {L2, F2})

  findEquivalence({L1, F1}, {L3, F3})
  
  FT1 = factorsByType det hessian F1
  FT2 = factorsByType det hessian F2
  FT3 = factorsByType det hessian F3

  md = hessianMatches(F1, F2)
  md = md | matchingData{L1 => L2}
  tryEquivalences(md, RQ, (A,phi))

  md = hessianMatches(F1, F3)
  md = md | matchingData{L1 => L3} | singularPointMatches(F1, F3)
  tryEquivalences(md, RQ, (A,phi))

  md = hessianMatches(F2, F3)
  md = md | matchingData{L2 => L3} 
  tryEquivalences(md, RQ, (A,phi))



  srcs = flatten {FT1#{1,1}, FT1#{2,1}, L1}
  targs = flatten \ (cartesian {signedPermutations FT2#{1,1}, signedPermutations FT2#{2,1}, {L2}})
  for t in targs list (
    A % trim equivalenceIdeal(srcs, t, RQ, (A,phi))
    )

  srcs = flatten {FT1#{1,1}, FT1#{2,1}, L1}
  targs = flatten \ (cartesian {signedPermutations FT3#{1,1}, signedPermutations FT3#{2,1}, {L3}})
  for t in targs list (
    A0 = A % trim equivalenceIdeal(srcs, t, RQ, (A,phi));
    A0 = try lift(A0, QQ) else continue;
    A0 = try lift(A0, ZZ) else continue;
    f := map(RQ, RQ, (transpose A0)**QQ);
    {A0, f(F1) - F3, f(L1) - L3}
    )

  srcs = flatten {FT2#{1,1}, FT2#{2,1}, L2}
  targs = flatten \ (cartesian {signedPermutations FT3#{1,1}, signedPermutations FT3#{2,1}, {L3}})
  for t in targs list (
    A0 = A % trim equivalenceIdeal(srcs, t, RQ, (A,phi));
    A0 = try lift(A0, QQ) else continue;
    A0 = try lift(A0, ZZ) else continue;
    f := map(RQ, RQ, (transpose A0)**QQ);
    {A0, f(F2) - F3, f(L2) - L3}
    )

  singularPointMatches(F1, F2)
  singularMatches(F1, F2)
  equivalenceIdeal({transpose matrix{{2,0,0,1,1}}}, {transpose matrix{{0,0,-2,2,1}}}, RQ, (A,phi))
  equivalenceIdeal({transpose matrix{{2,0,0,1,1}}}, {transpose matrix{{0,0,2,-2,-1}}}, RQ, (A,phi))
  
  md = hessianMatches(F1, F2) | singularMatches(F1, F2) | matchingData {L1 => L2, F1 => F2}
  matches md
  selectLinear md
  tryEquivalences(oo, RQ, (A,phi))
  tryEquivalences(md, RQ, (A,phi))
    -- XXX
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  A = extendToMatrix{10,15,6}
  assert(A * transpose matrix{{10,15,6}} == transpose matrix{{0,0,1}})
  assert(abs det A == 1)

  A = extendToMatrix{10,15,1}
  assert(A * transpose matrix{{10,15,1}} == transpose matrix{{0,0,1}})
  assert(abs det A == 1)

  A = extendToMatrix{10,15,0}
  assert(A * transpose matrix{{10,15,0}} == transpose matrix{{0,0,5}}) -- 5 is gcd!
  assert(abs det A == 1)
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  -- trivial case
  R = ZZ[]
  (A, phi) = genericLinearMap R
  assert(numRows A == 0 and numcols A == 0)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)

  R = ZZ[a..d]
  (A, phi) = genericLinearMap R
  assert(numRows A == 4 and numcols A == 4)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)

  (A, phi) = genericLinearMap(R, Variable => symbol s)
  assert(numRows A == 4 and numcols A == 4)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)
  use ring A
  assert(A_(0,0) == s_(1,1))  
  
  RQ = QQ (monoid R)
  (A, phi) = genericLinearMap RQ
  assert(numRows A == 4 and numcols A == 4)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)
  
  R = ZZ[a]  
  (A, phi) = genericLinearMap R
  assert(numRows A == 1 and numcols A == 1)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)

  R = ZZ/101[a..d]
  (A, phi) = genericLinearMap R
  U = target phi
  assert(source phi === U)
  assert(ring A === coefficientRing U)
  for i from 0 to 3 do 
    assert(phi U_i == (A^{i} * (transpose vars U))_(0,0))

  R = ZZ[a..d]
  (A, phi) = genericLinearMap R
  U = target phi
  assert(source phi === U)
  assert(ring A === coefficientRing U)
  for i from 0 to 3 do 
    assert(phi U_i == (A^{i} * (transpose vars U))_(0,0))

  R = QQ[a..e]
  (A, phi) = genericLinearMap R
  U = target phi
  assert(source phi === U)
  assert(ring A === coefficientRing U)
  for i from 0 to numgens R - 1 do 
    assert(phi U_i == (A^{i} * (transpose vars U))_(0,0))
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  RZ = ZZ[a,b,c]
  F1 = 3*a^2*b-3*a*b^2+b^3+6*a^2*c-6*a*c^2+2*c^3
  L1 = 24*a+10*b+8*c
  F2 = F1
  L2 = L1

  use ring L1
  findEquivalence({L1, F1}, {L1 - 2*a, F1+2*a^3})

  M1 = a
  M2 = b
  M3 = c
  N1 = a+b
  N2 = a+c
  N3 = 2*b+c

  md = matchingData {
    F1 => F2, 
    L1 => L2,
    {Permutations, {ideal M1, ideal M2}, {ideal N1, ideal N2}},
    {SignedPermutations, {M1,M2,M3}, {N1,N2,N3}},
    matrix{{1,2,3}} => matrix{{1,-2,1}}
    }
  assert isWellDefined md
  matches md
  netList first matches md
  netList last matches md

  md = matchingData {
    F1 => F2, 
    L1 => L2,
    {Permutations, {ideal(M1, M2^2), ideal(M2,M3^2)}, {ideal(N1,N2^2), ideal(N2,N3^2)}},
    {SignedPermutations, {M1,M2,M3}, {N1,N2,N3}},
    matrix{{1,2,3}} => matrix{{1,-2,1}}
    }
  assert isWellDefined md
  matches md
  netList first matches md
  netList last matches md

  md1 = selectLinear md
  
  assert try (matchingData {
    {F1, F2},
    {Permutations, {L1}, {L2}},
    {SignedPermutations, {M1,M2,M3}, {N1,N2,N3}},
    {matrix{{1,2,3}}, matrix{{1,-2,1}}}
    }; false
  ) else true

  assert try (
      matchingData {
          F1 => F2,
          L1 => L2,
          {M1,M2,M3} => {{N1,N2,N3,N1}},
          matrix{{1,2,3}} => matrix{{1,-2,1}}
          };
      false
  ) else true;

///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  assert(# signedPermutations{1,2,3} === 48)
  assert(# permutations{1,2,3} === 6)
  assert(cartesian{{1,2}} == {{1}, {2}})
  assert(cartesian{{1,2}, {3,4}} === {{1, 3}, {1, 4}, {2, 3}, {2, 4}})
  assert(cartesian{{1,2},{3,4},{5,6}} === {{1, 3, 5}, {1, 3, 6}, {1, 4, 5}, {1, 4, 6}, {2, 3, 5}, {2, 3, 6}, {2, 4, 5}, {2, 4, 6}})
  assert(cartesian{{1,2},{3},{5,6}} === {{1, 3, 5}, {1, 3, 6}, {2, 3, 5}, {2, 3, 6}})
  assert(cartesian{{1,2},{},{5,6}} === {})
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ)
  (A,phi) = genericLinearMap RQ

  -- Here is how we construct these polynomials (using h11=3 database).  
  --(L1, F1) = (c2Form Xs#(26,0), cubicForm Xs#(26,0))
  --(L2, F2) = (c2Form Xs#(37,0), cubicForm Xs#(37,0))
  use RQ  
  (L1, F1) = (10*a+28*b+26*c,a^3-3*a^2*b-3*a*b^2-2*b^3-3*a^2*c+6*a*b*c+6*b^2*c+3*a*c^2+6*b*c^2-c^3)
  (L2, F2) = (16*a+10*b+26*c,-2*a^3-3*a^2*b-3*a*b^2+b^3+6*a*b*c-3*b^2*c+6*a*c^2+3*b*c^2-c^3)
  
  MD = matchingData {
      L1 => L2,
      F1 => F2
      }
  (src, tars) = matches MD
  J = equivalenceIdeal(src, tars#0, RQ, (A,phi))
  A % J
  assert(first invertibleMatrixOverZZ(A, J) == INCONSISTENT)
  
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ)
  (A,phi) = genericLinearMap RQ

  -- Here is how we construct these polynomials (using h11=3 database).  
  --(L1, F1) = (c2Form Xs#(26,0), cubicForm Xs#(26,0))
  --(L2, F2) = (c2Form Xs#(37,0), cubicForm Xs#(37,0))

  use RZ
  (L1, F1) = (8*a-4*b+36*c,2*a^3-3*a^2*b-3*a*b^2+8*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*a*c^2)
  (L2, F2) = (-4*a+8*b+36*c,8*a^3-3*a^2*b-3*a*b^2+2*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*b*c^2)
  findEquivalence({L1, F1}, {L2, F2})

  --factorsByType det hessian F1
  
  MD = matchingData {
      L1 => L2,
      F1 => F2
      }
  (src, tars) = matches MD
  J = equivalenceIdeal(src, tars#0, RQ, (A,phi))
  A % J
  assert(first invertibleMatrixOverZZ(A, J) == CONSISTENT)

///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
///

end--

-* Development section *-
restart
needsPackage "IntegerEquivalences"
check "IntegerEquivalences"

uninstallPackage "IntegerEquivalences"
restart
installPackage "IntegerEquivalences"
viewHelp "IntegerEquivalences"
