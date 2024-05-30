-------------------------------------------------------
-- OLD: replaced with IntegerEquivalences package -----
-- TO BE REMOVED, do not use --------------------------
-------------------------------------------------------

debug needsPackage "StringTorics" -- let's arrange this so it doesn't need "debug"...
 -- debug needed for (at least): genericLinearMap.

-- MatchingData: is a list of elements each f the form
--  L => {M0, M1, ..., Ms}
-- where L is a list of:
--   RingElement: a polynomial in the original ring RZ or RQ.
--   Ideal: an ideal in the original ring.
--   Matrix: either a row vector or column vector, over ZZ or QQ (or RZ or RQ)
-- and each list Mi has the same length as L, and the same types of its elements.
MatchingData = new Type of List

matchingData = method()
matchingData List := LMs -> (
    ans := new MatchingData from LMs;
    if not isWellDefined ans then error "expected matching data to match.  Set `debugLevel=1` to investigate";
    ans
    )

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
                UNKNOWN
        ))
    )

-- helper function for (isWellDefined, MatchingData)
--   Input: one element (an Option, source and possible targets) of the matching data
--          which index this data sits at (used for error messages if debugLevel > 0)
--   Output: Boolean, whether this item is valid.
validMatchingItem = (LM,i) -> (
    L := LM#0;
    M := LM#1;
    -- M should be a list of lists, all same length and type as L.
    Ltype := itemType L;
    M' := if instance(M, List) then M else {M};
    Mtypes := for m in M' list itemType m;
    if any(Ltype, x -> x === UNKNOWN) then (
        if debugLevel > 0 then << "one element in source is unknown" << endl;
        return false;
        );
    for x in Mtypes do (
        if x =!= Ltype then (
            if debugLevel > 0 then << "in source " << i << ", target doesn't match source = " << Ltype << " obtaining instead " << x << endl;
            error "debug me";
            return false;
            );
        );
    true
    )
    
isWellDefined MatchingData := Boolean => LMs -> (
    for i from 0 to #LMs-1 do (
        LM := LMs#i;
        if not instance(LM, Option) then (
            if debugLevel > 0 then << "elements of list must be of the form L => M" << endl;
            return false;
            );
        if not validMatchingItem(LM,i) then return false;
        );
    true
    )

matches = method()
matches MatchingData := List => (MD) -> (
    targets := cartesian (MD/(x -> if instance(x#1, List) then x#1 else {x#1}));
    src := MD/first//flatten//toList;
    (src, targets/flatten//toList)
    -- XXX this is being tested now.
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
    TR := target phi; -- TODO: check: this is also source phi.
    if TR =!= source phi then error "expected ring map with same source and target";
    if #List1 == 0 then return ideal(0_TR);
    B := coefficientRing TR;
    toTR := map(TR, RQ, vars TR);
    toB := map(B, TR);
    List1 = List1/(I -> if ring I =!= RQ then toTR (sub(I, RQ)) else toTR I);
    List2 = List2/(I -> if ring I =!= RQ then toTR (sub(I, RQ)) else toTR I);
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
                ideal toB last coefficients sub(List2#i * (transpose A) - List1#i, TR)
            else
                ideal toB last coefficients sub((transpose A) * List1#i - List2#i, TR)
            )
        );
    sum ids
    )

tryEquivalencesV2 = method()
tryEquivalencesV2(MatchingData, Ring, Sequence) := (MD, RQ, Aphi) -> (
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

factorsByType = method()
factorsByType RingElement := HashTable => F -> (
    facs := factors F;
    faclist := for fx in facs list (fx#0, sum first exponents fx#1, fx#1);
    H := partition(x -> {x#0, x#1}, faclist);
    hashTable for k in keys H list k => for x in H#k list x_2
    )   

idealsByBettiV2 = method()
idealsByBettiV2(List, List) := List => (J1s, J2s) -> (
    H1 := partition(J -> betti gens J, J1s);
    H2 := partition(J -> betti gens J, J2s);
    if sort keys H1 =!= sort keys H2 then return {};
    for k in sort keys H1 list (
        H1#k => permutations H2#k
        )
    )

-*
hessianMatchesV2 = method()
hessianMatchesV2(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
    -- XXXX WORK ON THIS
    (A,phi) := Aphi;
    X1 := Xs#lab1;
    X2 := Xs#lab2;
    (L1, F1) := (c2Form X1, cubicForm X1);
    (L2, F2) := (c2Form X2, cubicForm X2);
    fac1 := factorsByType(det hessian F1);
    fac2 := factorsByType(det hessian F2);
    keyset1 := sort for k in sort keys fac1 list if k === {1,0} then (cH1 = fac1#k; continue) else k;
    keyset2 := sort for k in sort keys fac2 list if k === {1,0} then (cH2 = fac2#k; continue) else k;
    if keyset1 =!= keyset2 then error "expected same factors and their degrees";
    ans := for k in keyset1 list
        fac1#k => signedPermutations fac2#k
    checkFormat ans;
    ans
    )

singularLocusMatchesV2 = method()
singularLocusMatchesV2(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
    -- XXXX WORK ON THIS    
    -- returns a list of pairs:
    -- L1 => {L2a, L2b, ...}
    -- where L1 and each L2i have the same length,
    -- and every element is a RingElement, or an ideal
    -- possible todo: allow a point in H^2 as well.
    (A,phi) := Aphi;
    X1 := Xs#lab1;
    X2 := Xs#lab2;
    (L1, F1) := (c2Form X1, cubicForm X1);
    (L2, F2) := (c2Form X2, cubicForm X2);
    RQ := QQ (monoid ring L1);
    FQ1 = sub(F1, RQ);
    FQ2 = sub(F2, RQ);
    sing1 := trim saturate(ideal FQ1 + ideal jacobian FQ1);
    sing2 := trim saturate(ideal FQ2 + ideal jacobian FQ2);
    if sing1 == 1 then return {};
    comps1 := (decompose sing1)/trim;
    comps2 := (decompose sing2)/trim;
    ans := append(idealsByBettiV2(comps1, comps2), {sing1} => {sing2});
    checkFormat ans;
    ans
    )

equivalenceByHessian(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
    << "Calling new code" << endl;
    (A,phi) := Aphi;
    X1 := Xs#lab1;
    X2 := Xs#lab2;
    (L1, F1) := (c2Form X1, cubicForm X1);
    (L2, F2) := (c2Form X2, cubicForm X2);
    (list1, list2perms) := hessianMatches(lab1, lab2, Xs, Aphi);
    list1 = join(list1, {L1, F1});
    list2perms = for list2 in list2perms list join(list2, {L2, F2});
    tryEquivalences(list1, list2perms, Aphi)
    )
*-

end--

restart
load "FindEquivalencesV2.m2"

DBNAME = "./Databases/cys-ntfe-h11-3.dbm"
RZ = ZZ[a,b,c]
RQ = QQ (monoid RZ);
elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);

F1 = 3*a^2*b-3*a*b^2+b^3+6*a^2*c-6*a*c^2+2*c^3
L1 = 24*a+10*b+8*c
F2 = F1
L2 = L1
M1 = a
M2 = b
M3 = c
N1 = a+b
N2 = a+c
N3 = 2*b+c

debugLevel = 1
matchingData {
    F1 => F2, 
    L1 => L2,
    {M1,M2,M3} => {{N1,N2,N3}},
    matrix{{1,2,3}} => matrix{{1,-2,1}}
    }

TEST ///
  needs "FindEquivalencesV2.m2"
  assert(# signedPermutations{1,2,3} === 48)
  assert(# permutations{1,2,3} === 6)
  assert(cartesian{{1,2}} == {{1}, {2}})
  assert(cartesian{{1,2}, {3,4}} === {{1, 3}, {1, 4}, {2, 3}, {2, 4}})
  assert(cartesian{{1,2},{3,4},{5,6}} === {{1, 3, 5}, {1, 3, 6}, {1, 4, 5}, {1, 4, 6}, {2, 3, 5}, {2, 3, 6}, {2, 4, 5}, {2, 4, 6}})
  assert(cartesian{{1,2},{3},{5,6}} === {{1, 3, 5}, {1, 3, 6}, {2, 3, 5}, {2, 3, 6}})
///

TEST ///
  factor det hessian F1
  factorsByType det hessian F1
///

TEST ///
  restart
  needs "FindEquivalencesV2.m2"
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ)
  (A,phi) = genericLinearMap RQ

  -- Here is how we construct these polynomials (using h11=3 database).  
  --(L1, F1) = (c2Form Xs#(26,0), cubicForm Xs#(26,0))
  --(L2, F2) = (c2Form Xs#(37,0), cubicForm Xs#(37,0))
  
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
  restart
  needs "FindEquivalencesV2.m2"
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ)
  (A,phi) = genericLinearMap RQ

  -- Here is how we construct these polynomials (using h11=3 database).  
  -- (L1, F1) = (c2Form Xs#(190,0), cubicForm Xs#(190,0))
  -- (L2, F2) = (c2Form Xs#(193,0), cubicForm Xs#(193,0))
  
  (L1, F1) = (8*a-4*b+36*c,2*a^3-3*a^2*b-3*a*b^2+8*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*a*c^2)
  (L2, F2) = (-4*a+8*b+36*c,8*a^3-3*a^2*b-3*a*b^2+2*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*b*c^2)
  factorsByType det hessian F1
  
  MD = matchingData {
      L1 => L2,
      F1 => F2
      }
  (src, tars) = matches MD
  J = equivalenceIdeal(src, tars#0, RQ, (A,phi))
  A % J
  assert(first invertibleMatrixOverZZ(A, J) == CONSISTENT)
///


///
  restart
  needs "FindEquivalencesV2.m2"

  DBNAME = "./Databases/cys-ntfe-h11-3.dbm"
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ);
  elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);
  (A,phi) = genericLinearMap RQ

  labs = {(183, 0), (194, 1), (195, 0), (197, 0)}
  factorsByType det hessian cubicForm Xs#(labs#0)
  factorsByType det hessian cubicForm Xs#(labs#1)
///

 (L1, F1) = (c2Form Xs#(54,0), cubicForm Xs#(54,0))
 (L2, F2) = (c2Form Xs#(234,0), cubicForm Xs#(234,0))



  FQ1 = sub(F1, RQ)
  FQ2 = sub(F2, RQ)
  sing1 = trim saturate(ideal FQ1 + ideal jacobian FQ1)
  sing2 = trim saturate(ideal FQ2 + ideal jacobian FQ2)
  if sing1 == 1 then return {};
  comps1 := (decompose sing1)/trim;
  comps2 := (decompose sing2)/trim;
  idealsByBettiV2(comps1, comps2)
  ans := append(idealsByBettiV2(comps1, comps2), {sing1} => {sing2});


///

MD = matchingData {
    F1 => F2, 
    L1 => L2,
    {M1,M2,M3} => permutations {N1,N2,N3},
    matrix{{1,2,3}} => matrix{{1,-2,1}}
    }
(src,tar) = matches MD
(A, phi) = genericLinearMap RQ
J = equivalenceIdeal(src, tar#0, RQ, (A,phi))
tryEquivalencesV2(MD, RQ, (A,phi))

debugLevel = 1
MD1 = matchingData {
    a+2*b+3*c => a-b
    }
MD2 = matchingData {
    transpose matrix{{1,2,3}} => transpose matrix{{1,-1,0}}
    }

(src, tar) = matches MD1

J1 = equivalenceIdeal(src, tar#0, RQ, (A,phi))


(src, tar) = matches MD2
J2 = equivalenceIdeal(src, tar#0, RQ, (A,phi))
J1 == J2

MD3 = matchingData {
    a+2*b+3*c => a-b,
    a-b => c,
    c => a+b+c
    }
(src, tar) = matches MD3
J = equivalenceIdeal(src, tar#0, RQ, (A,phi))
J + ideal((det A) - 1)
A % oo



cartesian (MD/(x -> if instance(x#1, List) then x#1 else {x#1}))

MD = matchingData {
    F1 => F2, 
    L1 => L2,
    {M1,M2,M3} => signedPermutations {N1,N2,N3},
    matrix{{1,2,3}} => matrix{{1,-2,1}}
    }
cartesian (MD/(x -> if instance(x#1, List) then x#1 else {x#1}))

MD/(x -> x#0)//flatten
MD/(x -> x#1)
netList oo

