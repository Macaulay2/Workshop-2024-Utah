-------------------------------------------------------
-- OLD: replaced with IntegerEquivalences package -----
-- TO BE REMOVED, do not use --------------------------
-------------------------------------------------------
debug needsPackage "StringTorics"
-- routine to take saturation of singF, and its components, to their images

-- TODO: add in A as an argument, allow column matrices as well.
-- Then this can be used with findMaps from gv invariants too.
equivalenceIdeal = method()
equivalenceIdeal(List, List, RingMap) := Ideal => (List1, List2, phi) -> (
    TR := target phi; -- TODO: check: this is also source phi.
    if #List1 == 0 then return ideal(0_TR);
    RQ := ring List1_0;
    B := coefficientRing TR;
    toTR := map(TR, RQ, vars TR);
    toB := map(B, TR);
    List1 = List1/(I -> if ring I =!= RQ then toTR (sub(I, RQ)) else toTR I);
    List2 = List2/(I -> if ring I =!= RQ then toTR (sub(I, RQ)) else toTR I);
    -- List1 and List2 are lists with the same length, consisting of RingElement's and Ideal's.
    -- List1 and List2 should each have RingElement's and Ideal's in the same spot.
    ids := for i from 0 to #List1-1 list (
        if instance(List1#i, RingElement) then (
            if not instance(List2#i, RingElement) then 
                error "expected both lists to consist of Ideal's and RingElement's in the same spots";
            ideal toB (last coefficients(phi List1#i - List2#i))
        ) else if instance(List1#i, Ideal) then (
            if not instance(List2#i, Ideal) then 
                error "expected both lists to consist of Ideal's and RingElement's in the same spots";
            ideal toB (last coefficients((gens phi List1#i) % List2#i))
            )
        );
    sum ids
    )

equivalenceIdealMain = method()
equivalenceIdealMain(Thing, Thing, HashTable, Sequence) := Ideal => (lab1, lab2, Xs, Aphi) -> (
    X1 := Xs#lab1;
    X2 := Xs#lab2;
    (L1, F1) := (c2Form X1, cubicForm X1);
    (L2, F2) := (c2Form X2, cubicForm X2);
    (A,phi) := Aphi;
    equivalenceIdeal({L1, F1}, {L2, F2}, phi)
    )

factorsByType = method()
factorsByType RingElement := HashTable => F -> (
    facs := factors F;
    faclist := for fx in facs list (fx#0, sum first exponents fx#1, fx#1);
    H := partition(x -> {x#0, x#1}, faclist);
    hashTable for k in keys H list k => for x in H#k list x_2
    )   

-- TODO: make V2
idealsByBetti = method()
idealsByBetti(List, List) := List => (J1s, J2s) -> (
    H1 := partition(J -> betti gens J, J1s);
    H2 := partition(J -> betti gens J, J2s);
    if sort keys H1 =!= sort keys H2 then return {};
    list1 := flatten for k in sort keys H1 list H1#k;
    list2perms := cartesian for k in sort keys H2 list permutations H2#k;
    (list1, list2perms)
    )

idealsByBettiV2 = method()
idealsByBettiV2(List, List) := List => (J1s, J2s) -> (
    H1 := partition(J -> betti gens J, J1s);
    H2 := partition(J -> betti gens J, J2s);
    if sort keys H1 =!= sort keys H2 then return {};
    for k in sort keys H1 list (
        H1#k => permutations H2#k
    ))

-- TODO: this takes too much time, it seems.  If there are 4 factors, it creates 384 things to compute.
-- each is easy, but they add up.
hessianMatches = method()
hessianMatches(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
    (A,phi) := Aphi;
    X1 := Xs#lab1;
    X2 := Xs#lab2;
    (L1, F1) := (c2Form X1, cubicForm X1);
    (L2, F2) := (c2Form X2, cubicForm X2);
    fac1 := factorsByType(det hessian F1);
    fac2 := factorsByType(det hessian F2);
    keyset1 := sort for k in sort keys fac1 list if k === {1,0} then (cH1 = fac1#k; continue) else k;
    keyset2 := sort for k in sort keys fac2 list if k === {1,0} then (cH2 = fac2#k; continue) else k;
    if keyset1 === {{1,4}} then return "Do not use Hessian method on irreducible polynomials of degree > 4";    
    if keyset1 =!= keyset2 then error "expected same factors and their degrees";
    list1 := flatten for k in keyset1 list fac1#k;
    list2 := for k in keyset2 list fac2#k;
    list2perms := cartesian for fs in list2 list signedPermutations fs;
    (list1, list2perms)
    )

singularLocusMatches = method()
singularLocusMatches(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
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
    (list1, list2perms) := idealsByBetti(comps1, comps2);
    (join({L1, F1, sing1}, list1), cartesian{{{L2, F2, sing2}}, list2perms})
    )

checkFormatV2 = method()
-- TODO: format has changed.  Fix this function!
checkFormatV2(List, List) := (list1, list2perms) -> (
    -- TODO: list1 can be a single Ideal or RingElement or list of...
    -- TODO: 
    n := #list1;
    list1Type := for a in list1 list (
        if instance(a, RingElement) then RingElement
        else if instance(a, Ideal) then Ideal);
    for list2a in list2perms do (
        assert(#list2 == n);        
        list2Type := for a in list2a list (
            if instance(a, RingElement) then RingElement
            else if instance(a, Ideal) then Ideal);
        assert(list1Type === list2Type);
        );
    )
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
        fac1#k => signedPermutations fac2#k;
    --checkFormatV2 ans; -- not functional yet?  Check that...
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
    -- checkFormatV2 ans; -- not functional yet?  Check that...
    ans
    )

-- singularAndHessianMatches = method()
-- singularAndHessianMatches(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
--     (A,phi) := Aphi;
--     X1 := Xs#lab1;
--     X2 := Xs#lab2;
--     (L1, F1) := (c2Form X1, cubicForm X1);
--     (L2, F2) := (c2Form X2, cubicForm X2);
--     RQ := QQ (monoid ring L1);
--     FQ1 = sub(F1, RQ);
--     FQ2 = sub(F2, RQ);
--     fac1 := factorsByType(det hessian F1);
--     fac2 := factorsByType(det hessian F2);
--     keyset1 := sort for k in sort keys fac1 list if k === {1,0} then continue else k;
--     keyset2 := sort for k in sort keys fac2 list if k === {1,0} then continue else k;
--     --if keyset1 === {{1,4}} then return "Do not use Hessian method on irreducible polynomials of degree > 4";    
--     if keyset1 =!= keyset2 then error "expected same factors and their degrees";
--     list1 := flatten for k in keyset1 list fac1#k;
--     list2 := for k in keyset2 list fac2#k;
--     list2perms := cartesian for fs in list2 list signedPermutations fs;
--     hessPart := (list1, list2perms);
--     sing1 := trim saturate(ideal FQ1 + ideal jacobian FQ1);
--     sing2 := trim saturate(ideal FQ2 + ideal jacobian FQ2);
--     if sing1 == 1 then return {};
--     comps1 := (decompose sing1)/trim;
--     comps2 := (decompose sing2)/trim;
--     (list1, list2perms) = idealsByBetti(comps1, comps2);
--     singPart := (join({L1, F1, sing1}, list1), cartesian{{{L2, F2, sing2}}, list2perms});
--     secondPart := flatten for a in singPart#1 list for b in hessPart#1 list (a | b);
--     ((singPart#0 | hessPart#0), secondPart)
--     )

equivalenceByHessian = method()
-- equivalenceByHessian(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
--     (A,phi) := Aphi;
--     X1 := Xs#lab1;
--     X2 := Xs#lab2;
--     (L1, F1) := (c2Form X1, cubicForm X1);
--     (L2, F2) := (c2Form X2, cubicForm X2);
--     (list1, list2perms) := hessianMatches(lab1, lab2, Xs, Aphi);
--     badJs := {};
--     ans := for i from 0 to #list2perms-1 do (
--         << "doing " << i << endl;
--         J := trim equivalenceIdeal(join({L1,F1},list1), join({L2,F2},list2perms#i), phi);
--         if J == 1 then continue;
--         << "doing " << i << " not <1>" << endl;
--         A0 := A % J;
--         -- check A0 is integer
--         if support A0 =!= {} then badJs = append(badJs, J);
--         try (A0 = lift(A0, ZZ)) else continue;
--         if det A0 != 1 and det A0 != -1 then continue;
--         break A0
--         );
--     if ans === null and #badJs > 0 then return badJs;
--     ans
--     )

invertibleMatrixOverZZ = method()
invertibleMatrixOverZZ(Matrix, Ideal) := Sequence => (A, J) -> (
    -- returns (determinacy, A0), or (INCONSISTENT, null) or ...?
    if J == 1 then 
        (INCONSISTENT, null)
    else (
        A0 := A % J;
        detA0 := (det A0) % J;
        if liftable(detA0, ZZ) and all(flatten entries A0, f -> liftable(f, ZZ)) then
            (CONSISTENT, sub(A0, ZZ));
        (INDETERMINATE, J)
        )
    )

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

  findMatrix5 = method(Options => {Pair => {0,1}})
  findMatrix5(CalabiYauInToric, CalabiYauInToric, List, Sequence) := opts -> (X1, X2, pts2, Aphi) -> (
      -- pts2 is a list of possible entries {p,q} for the first 2 columns of a matrix A0
      (A, phi) := Aphi;
      F1 := cubicForm X1;
      L1 := c2Form X1;
      F2 := cubicForm X2;
      L2 := c2Form X2;
      F1 = F1 // polynomialContent F1;
      F2 = F2 // polynomialContent F2;
      L1 = L1 // polynomialContent L1;
      L2 = L2 // polynomialContent L2;
      L1' := sub(L1, target phi);
      L2' := sub(L2, target phi);
      F1' := sub(F1, target phi);
      F2' := sub(F2, target phi);
      result := null;
      count := 0;
      i1 := opts.Pair#0;
      j1 := opts.Pair#1;
      for pq in pts2 do (
          count = count+1;
          A1 := matrix {for i from 0 to 4 list 
            if i == i1 then transpose matrix{pq#0}
            else if i == j1 then transpose matrix{pq#1}
            else A_{i}
            };
          --A1 = (transpose matrix{pq#0, pq#1}) | A_{2,3,4};
          psi1 := map(target phi, target phi, transpose A1);
          Ja := ideal last coefficients (psi1 L1' - L2');
          Jb := ideal last coefficients (psi1 F1' - F2');
          J := sub(Ja + Jb, ring A1);
          (cons, A0) := invertibleMatrixOverZZ(A1, J);
          if cons === CONSISTENT then (result = (cons, A0); break);
          if cons === INCONSISTENT then (continue);
          if cons === INDETERMINATE then (
              error "debug me";
              << "found indeterminate case" << endl; continue
              );
          );
      if result === null then return null;
      A0 := result#1;
      if isEquivalent(X1, X2, transpose A0) then << "** Found one ** " << toString A0 << endl;
      << "took " << count << " tries" << endl;
      result
      )

tryEquivalences = method()
tryEquivalences(List, List, Sequence) := (list1, list2perms, Aphi) -> (
    badJs := {};
    inconsistentMatrix := null;
    for i from 0 to #list2perms-1 do (
        << "doing " << i << endl;
        J := trim equivalenceIdeal(list1, list2perms#i, phi);
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

-- invertibleMatrixOverZZ(A0, J):
--   check to see if A0 is an invertible matrix over ZZ (modulo J).
--   (det A0 % J) == 1, or -1.
--   if not return (INDETERMINATE, null)
--   if support A0 is {}, then lift to QQ, then to ZZ.
--     if lifts: return (CONSISTENT, A0).
--     
--  if A0 is a matrix over ZZ:
--    create it over ZZ
--    if det A0 != 1, -1 (perhaps mod J?) then return (INCONSISTENT, null)
--    make sure (A - A0) % J == 0.  If so, then return (CONSISTENT, A0)

-- I want a function which, given an equivalence ideal, returns either INCONSISTENT, CONSISTENT (solution over ZZ).
-- if J == 1 then INCONSISTENT
-- else:
--  try A0 = A % J.
--  (consistency, A0) = invertibleMapOverZZ(A0, J)
--  if consistency != INDETERMINATE then return (consistency, A0);
--  if A0 is not a matrix over ZZ:
--    cJ = decompose J -- we need one of these to succeed
--    for each j in cJ do:
--      A0 := A % j
--      if A0 is a matrix over ZZ
--  
equivalenceBySingularLocus = method()
-- equivalenceBySingularLocus(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
--     (A,phi) := Aphi;
--     X1 := Xs#lab1;
--     X2 := Xs#lab2;
--     (L1, F1) := (c2Form X1, cubicForm X1);
--     (L2, F2) := (c2Form X2, cubicForm X2);
--     (list1, list2perms) := singularLocusMatches(lab1, lab2, Xs, Aphi);
--     badJs := {};
--     ans := for i from 0 to #list2perms-1 do (
--         J := trim equivalenceIdeal(join({L1,F1},list1), join({L2,F2},list2perms#i), phi);
--         if J == 1 then continue;
--         A0 := A % J;
--         -- check A0 is integer
--         if support A0 =!= {} then badJs = append(badJs, J);
--         try (A0 = lift(A0, ZZ)) else continue;
--         if det A0 != 1 and det A0 != -1 then continue;
--         break A0
--         );
--     if ans === null and #badJs > 0 then return badJs;
--     ans
--     )

equivalenceBySingularLocus(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
    (A,phi) := Aphi;
    X1 := Xs#lab1;
    X2 := Xs#lab2;
    (L1, F1) := (c2Form X1, cubicForm X1);
    (L2, F2) := (c2Form X2, cubicForm X2);
    (list1, list2perms) := singularLocusMatches(lab1, lab2, Xs, Aphi);
    list1 = join(list1, {L1, F1});
    list2perms = for list2 in list2perms list join(list2, {L2, F2});
    tryEquivalences(list1, list2perms, Aphi)
    )

-- This is not working well...  i.e. no working examples yet!
-- equivalenceBySingularLocusAndHessian = method()
-- equivalenceBySingularLocusAndHessian(Thing, Thing, HashTable, Sequence) := List => (lab1, lab2, Xs, Aphi) -> (
--     << "Calling new code" << endl;
--     (A,phi) := Aphi;
--     X1 := Xs#lab1;
--     X2 := Xs#lab2;
--     (L1, F1) := (c2Form X1, cubicForm X1);
--     (L2, F2) := (c2Form X2, cubicForm X2);
--     (list1, list2perms) := singularAndHessianMatches(lab1, lab2, Xs, Aphi);
--     tryEquivalences(list1, list2perms, Aphi)
--     )

equivalenceByGVCone = method(Options => {DegreeLimit => 15})
equivalenceByGVCone(Thing, Thing, HashTable, Sequence) := Matrix => opts -> (lab1, lab2, Xs, Aphi) -> (
    -- return value of null means they might still be the same, we just don't know.
    X1 := Xs#lab1;
    X2 := Xs#lab2;
    gv1 := partitionGVConeByGV(X1, opts);
    gv2 := partitionGVConeByGV(X1, opts);
    if gv1 === null or gv2 === null then return null;
    findLinearMaps(gv1, gv2) -- what does this return?  Not what we want, certainly!
    -- let's put in the L1->L2, F1->F2 equations, maybe det+1, det-1 too, and stop when we find a match.
    )
-*
equivalenceIdealSingularSet = method()
equivalenceIdealSingularComponents = method()

  equivalenceIdealsSingularLocus = method()
  equivalenceIdealsSingularLocus(Sequence, Sequence, RingMap) := List => (LF1, LF2, phi) -> (
      -- returns a list of ideals in target phi (=== source phi).
      (L1, F1) := LF1;
      (L2, F2) := LF2; -- TODO: the rings of all 4 of these polynomials should be the same, and be RZ.
      RZ := ring L1;
      RQ := QQ (monoid RZ);
      toRQ := map(RQ, RZ);
      -- Get singular locus and components (over QQ) of F1.
      singF1 := saturate(ideal F1 + ideal jacobian F1);
      singF1Q := trim toRQ singF1;
      compsS1 := (decompose singF1Q)/trim;
      -- Get singular locus and components (over QQ) of F2.
      singF2 := saturate(ideal F2 + ideal jacobian F2);
      singF2Q := trim toRQ singF2;
      compsS2 := (decompose singF2Q)/trim;
      -- From invariants, we expect singF1Q, singF2Q to have the same degree generators (and their number).
      -- TODO: should we check this?
      -- Similarly, we expect compsS1, compsS2 to have the same number of components, and same degree gens each.
      if betti singF2Q != betti singF1Q then (
          << "degrees of generators of sing loci ideals differ: returning no ideals" << endl;
          return {};
          );
      error "debug me";
      )

  equivalenceIdealsSingularLocus(Thing, Thing, HashTable, RingMap) := List => (lab1, lab2, Xs, phi) -> (
      X1 := Xs#lab1;
      X2 := Xs#lab2;
      (L1, F1) := (c2Form X1, cubicForm X1);
      (L2, F2) := (c2Form X2, cubicForm X2);
      equivalenceIdealsSingularLocus((L1,F1), (L2,F2), phi)
      )
*-

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
    if #Ls == 1 then return for p in Ls#0 list p;
    Ls1 := cartesian drop(Ls,1);
    flatten for p in Ls#0 list for q in Ls1 list join(p, q)
    )

TEST ///
  R = ZZ[a,b,c,d]
  assert(allSigns{1,2,3} === {{1, 2, 3}, {-1, 2, 3}, {1, -2, 3}, {-1, -2, 3}, {1, 2, -3}, {-1, 2, -3}, {1, -2, -3}, {-1, -2, -3}})
  assert(allSigns{} === {{}})
  assert(allSigns{a+b,a+c} === {{a+b, a+c}, {-a-b, a+c}, {a+b, -a-c}, {-a-b, -a-c}})

  assert(signedPermutations {a,b} === {{a, b}, {-a, b}, {a, -b}, {-a, -b}, {b, a}, {-b, a}, {b, -a}, {-b, -a}})

  assert(
  cartesian{signedPermutations{a,b}, permutations{c,d}}
  ==
  {{a, b, c, d}, {a, b, d, c}, {-a, b, c, d}, {-a, b, d, c}, 
      {a, -b, c, d}, {a, -b, d, c}, {-a, -b, c, d}, {-a, -b, d, c}, 
      {b, a, c, d}, {b, a, d, c}, {-b, a, c, d}, {-b, a, d, c}, 
      {b, -a, c, d}, {b, -a, d, c}, {-b, -a, c, d}, {-b, -a, d, c}}
  )
  
///

-------------------------------
-- Given a list of labels, and the hash table Xs, and (A,phi) (do we really need this last one??)
-- Returns a list of 2 lists:
-- The firt is a list {L0, L1, ..., Ld}, where each Li is a list of labels (and possible (label, matrix).
--   such that for each i, the labels in Li are all equivalent (there is a matrix over ZZ of unit det which
--   maps one to the other,
--   for each i != j, any label in Li is NOT equivalent to any in Lj.
-- The second is a list of labels where we don't know if they are equiv or not equiv to any one of the above sets.
--
-- Maybe: allow an initial result to be created
separateAndCombineByFcn = method()
separateAndCombineByFcn(List, HashTable, Sequence, Function) := (Ls, Xs, Aphi, F) -> (
    -- F: (X,Y) --> List (of {CONSISTENT, ...}, {INCONSISTENT, ...}, {INDETERMINATE, ...}.
    resultList := new MutableList;
    resultUnknowns := {};
    for lab in Ls do (
        known := false;
        hasbad := false;
        i := 0;
        for i from 0 to #resultList-1 do (
            lab0 := first resultList#i;
            ans := F(lab0, lab, Xs, Aphi);
            if first ans === CONSISTENT then (
                -- then we know this one, and can stash it, and go on to the next one.
                << "-- note: " << lab << " and " << lab0 << " are equivalent via: " << ans#1 << endl;
                resultList#i = append(resultList#i, lab);
                known = true;
                break;
                )
            else if first ans =!= INCONSISTENT then (
                hasbad = true;
                );
            );
        if not known then (
            if not hasbad then (
                << "-- note: " << lab << " is not equivalent to ones before" << endl;                
                resultList#(#resultList) = {lab}
                )
            else (
                << "-- note: " << lab << " is inconclusive" << endl;
                resultUnknowns = append(resultUnknowns, lab);
                )
            );
        );
    {toList resultList, resultUnknowns}
    )

end--

-- Tests of this code, 12 Jan 2023.
TEST ///
-*
  restart
  needsPackage "StringTorics"
*-
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  DBNAME = "./Databases/cys-ntfe-h11-5.dbm"
  RZ = ZZ[a,b,c,d,e]
  RQ = QQ (monoid RZ);
  needs "../FindEquivalence.m2"
  (A,phi) = genericLinearMap RQ
  elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);
  REPS = value get "inequiv-reps-h11-5"
  SETS = new HashTable from {
      {{{4,1},{4,1}},{{5,1}}} => {50,55,65,66,70,72,73,80,85,91,92,95,96,100,103,116,117,186,193}, 
      {{{5,1}},{{5,1}}} => {7,9,10,14,22,23,24,27,31,36,37,49,51,56,71,93,99,101,102,108,115,130}, 
      {{{4,1},{4,1},{4,1}},{{5,1}}} => {52,77,155,171,180,217}, 
      {{{5,1}},{{1,1},{1,1},{3,1}}} => {18,25,32,35,41,59,60,75,111,112,142,143,153,154,178,188,208,215,221,226,227}, 
      {{{4,2}},{{1,1},{2,1},{2,1}}} => {79}, 
      {{{4,1}},{{1,1},{1,1},{1,1},{1,2}}} => {118,129,137,138,144,145,190,191,200}, 
      {{{5,1}},{{1,1},{1,1},{1,1},{1,1},{1,1}}} => {8,28,61,86,88,89,120,131,132,136,146,149,158,160,166,173,201,202,205,207,211,213,218,219,220,224,225}, 
      {{{4,1},{4,1},{4,1}},{{1,1},{1,1},{1,1},{1,1},{1,1}}} => {124,156,168}, 
      {{{4,1},{4,1}},{{1,1},{4,1}}} => {44,47,177}, 
      {{{3,2},{4,1}},{{1,1},{1,2},{2,1}}} => {135,164}, 
      {{{4,2}},{{5,1}}} => {6,15}, 
      {{{5,1}},{{1,1},{4,1}}} => {0,3,4,5,11,13,17,21,29,30,34,38,39,40,42,45,57,68,81,83,90,94,98,104,107,109,110,114,119,128,140,163,169,172,195,212,223}, 
      {{{4,2}},{{1,1},{1,1},{3,1}}} => {48,105}, 
      {{{4,1},{4,1},{4,1}},{{1,1},{4,1}}} => {62,63,127,151}, 
      {{{3,2},{4,1}},{{1,1},{4,1}}} => {53,54,76,122,134}, 
      {{{3,2},{4,1}},{{1,2},{3,1}}} => {78}, 
      {{{4,1}},{{5,1}}} => {1,2,19,20,26,64,67,97,121,123}, 
      {{{4,2}},{{1,1},{4,1}}} => {69,82,113,147,152,189,199,206,214}, 
      {{{4,1}},{{1,1},{1,2},{2,1}}} => {84,139,161,162,197}, 
      {{{4,1}},{{1,2},{3,1}}} => {159,179,181,182,183,184,192,194,198}, 
      {{{4,1}},{{1,1},{4,1}}} => {12,43,187,216}, 
      {{{3,2}},{{1,1},{1,2},{2,1}}} => {125,126,148,165,175,185,196,204}, 
      {{{4,1},{4,1},{4,1},{4,1}},{{5,1}}} => {222}, 
      {{{3,2}},{{1,1},{4,1}}} => {74,106,141,157,174,210}, 
      {{{3,2}},{{1,2},{3,1}}} => {46}, 
      {{{4,1},{4,1},{4,1},{4,1}},{{1,1},{4,1}}} => {16}, 
      {{{4,1},{4,1}},{{1,1},{1,1},{1,1},{2,1}}} => {133,170,209}, 
      {{{5,1}},{{1,1},{1,1},{1,1},{2,1}}} => {33,58}, 
      {{{2,2}},{{1,3},{2,1}}} => {87}, 
      {{{4,2}},{{1,1},{1,1},{1,1},{2,1}}} => {150,167,176,203}}  
  (sort keys SETS)/(k -> {k#0, k#1, #SETS#k, SETS#k})//netList
  -- column one: structure of sing locus
  -- column 2: factorization of Hessian.
  -- column 3: # of REP sets
  -- column 4: list of indices into REPS.
  
  -- We can use this collection to test FindEquivalence.m2 code:
  -- Goal: for each set of labels REPS#i, determine if these are the same topology or different.
  -- Note: some are very easy, some I can't yet do.
///

TEST ///
-*
  restart
  needsPackage "StringTorics"
*-
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  RZ = ZZ[a,b,c,d,e]
  RQ = QQ (monoid RZ);
  needs "../FindEquivalence.m2"
  load "../FindEquivalence.m2"
  (A,phi) = genericLinearMap RQ
  elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);
  
  -- Test of: idealsByBettiV2 -- XXX
  lab1 = (879,0)
  lab2 = (910,0)
  X1 = Xs#lab1;
  X2 = Xs#lab2;
  (L1, F1) = (c2Form X1, cubicForm X1);
  (L2, F2) = (c2Form X2, cubicForm X2);
  --RQ := QQ (monoid ring L1);
  FQ1 = sub(F1, RQ);
  FQ2 = sub(F2, RQ);
  sing1 = trim saturate(ideal FQ1 + ideal jacobian FQ1);
  sing2 = trim saturate(ideal FQ2 + ideal jacobian FQ2);
  comps1 = (decompose sing1)/trim;
  comps2 = (decompose sing2)/trim;
  (list1, list2perms) = idealsByBetti(comps1, comps2)

  comps1 = {ideal (d, c, b, a), ideal (e, d, b, a)}
  comps2 = {ideal (d, c, b, a), ideal (e, d, b, a)}
  val = idealsByBettiV2(comps1, comps2) 
  ans = {{comps1_0, comps1_1} => {comps2, reverse comps2}} -- this is because comps1, comps2 have 2 ideals of same type.
  assert(ans === val)

  singularLocusMatches(lab1, lab2, Xs, (A,phi))
  singularLocusMatchesV2(lab1, lab2, Xs, (A,phi))
  netList for x in oo list {x#0, netList x#1}

  idealsByBettiV2 -- 
///

TEST ///
-*
  restart
  needsPackage "StringTorics"
*-
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  RZ = ZZ[a,b,c,d,e]
  RQ = QQ (monoid RZ);
  needs "../FindEquivalence.m2"
  load "../FindEquivalence.m2"
  (A,phi) = genericLinearMap RQ
  elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);

-- Here is a potentially hard one:
-- Actually, we can find a map between them, this is done below.
  (lab1, lab2) = toSequence {(133, 0), (165, 0)} -- EQUIVALENT!
  X1 = Xs#lab1;
  X2 = Xs#lab2;
  (L1, F1) = (c2Form X1, cubicForm X1);
  (L2, F2) = (c2Form X2, cubicForm X2);
  FQ1 = sub(F1, RQ);
  FQ2 = sub(F2, RQ);
  sing1 = trim saturate(ideal FQ1 + ideal jacobian FQ1);
  sing2 = trim saturate(ideal FQ2 + ideal jacobian FQ2);
  comps1 = (decompose sing1)/trim;
  comps2 = (decompose sing2)/trim;
  classifyExtremalCurves X1, classifyExtremalCurves X2

  -- These give no information
  L1
  L2
  Rp = (ZZ/3) (monoid RQ)
  F1p = sub(F1, Rp)
  F2p = sub(F2, Rp)
  factor det hessian F1p
  factor det hessian F2p
  sing1 = trim saturate(ideal F1p + ideal jacobian F1p);
  sing2 = trim saturate(ideal F2p + ideal jacobian F2p);

  needsPackage "LLLBases"
  (cL1, M1) = gcdLLL ((listForm L1)/last)
  f1 = map(RZ, RZ, transpose M1)
  f1 L1
  L1 = f1 L1
  F1 = f1 F1

  (cL2, M2) = gcdLLL ((listForm L2)/last)
  f2 = map(RZ, RZ, transpose M2)
  f2 L2
  L2 = f2 L2
  F2 = f2 F2

  use RZ
  sing1 = trim saturate(ideal F1 + ideal jacobian F1);
  sing2 = trim saturate(ideal F2 + ideal jacobian F2);
  factor leadCoefficient sing1_0
  F1e = sub(F1, e => 0)
  F2e = sub(F2, e => 0)
  SZ = ZZ[a,b,c,d]
  SQ = QQ (monoid SZ)
  F1e = sub(F1e, SZ)
  F2e = sub(F2e, SZ)
  factor det hessian F1e
  factor det hessian F2e
  polynomialContent det hessian F1e
  polynomialContent det hessian F2e
  (A1, phi1) = genericLinearMap SQ
  J = ideal last coefficients(phi1 (sub(F1e, target phi1)) - sub(F2e, target phi1))
  --gbTrace=3
  --gb J;
  sing1 = saturate (ideal F1e + ideal jacobian F1e)
  sing2 = saturate (ideal F2e + ideal jacobian F2e)
  factor leadCoefficient sing1_0
  factor leadCoefficient sing2_0
  sing1 : (sing1 : 17)
  sing2 : (sing2 : 17)

  sing1 : (sing1 : 19)
  sing2 : (sing2 : 19)

  sing1 : (sing1 : 173)
  sing2 : (sing2 : 173)

  sing1 : (sing1 : 89)
  sing2 : (sing2 : 89)

  Sp = (ZZ/89) (monoid SZ)
  (Ap, phip) = genericLinearMap Sp
  sub(F1e, target phip)
  J = ideal last coefficients(phip (sub(F1e, target phip)) - sub(F2e, target phip))
  see J  
  J1 = trim sub(  sing1 : (sing1 : 89), Sp)
  J2 = trim sub(  sing2 : (sing2 : 89), Sp)
  J' = ideal last coefficients(phip (gens sub(J1, target phip)) % gens sub(J2, target phip))
  J = J + J'
  J = sub(J, coefficientRing ring J)
  elapsedTime gbJ = groebnerBasis J;
  see ideal gbJ
  comps = decompose ideal gbJ -- 8 components, but 4 don't have a point over ZZ/89.  The other 4 have a unique point.
  A1 = sub(Ap % comps_0, ZZ)
  fp = map(SZ, SZ, transpose A1)
  fp F1e - F2e

  -- let's try map that leaves e fixed:
  use RZ
  fR = map(RZ, RZ, {-b-c, c, -a, d, e})
  fR F1 - F2
  
  U = source phi
  T = coefficientRing U
  use T; use U
  fR = map(U, U, {-b-c + t_(1,1)*e, c + t_(1,2)*e, -a + t_(1,3)*e, d + t_(1,4)*e, e})
  trim ideal last coefficients(fR sub(F1, U) - sub(F2, U))
  fU = map(U, U, fR.matrix % oo)
  fU sub(F1, U) - sub(F2, U)
    
  A1 = sub(A1, coefficientRing source phi)
  A2 = (transpose A1 | transpose matrix{{t_(1,1), t_(1,2), t_(1,3), t_(1,4)}}) || matrix{{0,0,0,0,1}}
  psi = map(source phi, source phi, transpose A2)
  use target psi
  psi (e)    
  psi F1
  psi sub(F1, source psi) - sub(F2, source psi)
///


TEST ///
-*
  restart
  needsPackage "StringTorics"
*-
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  RZ = ZZ[a,b,c,d,e]
  RQ = QQ (monoid RZ);
  needs "../FindEquivalence.m2"
  load "../FindEquivalence.m2"
  (A,phi) = genericLinearMap RQ
  elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);

-- Here is a potentially hard one:
-- Actually, we can find a map between them, this is done below.
  (lab1, lab2) = toSequence {(159, 2), (219, 2)}
  X1 = Xs#lab1;
  X2 = Xs#lab2;
  (L1, F1) = (c2Form X1, cubicForm X1);
  (L2, F2) = (c2Form X2, cubicForm X2);
  FQ1 = sub(F1, RQ);
  FQ2 = sub(F2, RQ);
  sing1 = trim saturate(ideal FQ1 + ideal jacobian FQ1)
  sing2 = trim saturate(ideal FQ2 + ideal jacobian FQ2)
  comps1 = (decompose sing1)/trim
  comps2 = (decompose sing2)/trim
  classifyExtremalCurves X1, classifyExtremalCurves X2

  needsPackage "LLLBases"
  (cL1, M1) = gcdLLL ((listForm L1)/last)
  f1 = map(RZ, RZ, transpose M1)
  f1 L1
  L1 = f1 L1
  F1 = f1 F1

  (cL2, M2) = gcdLLL ((listForm L2)/last)
  f2 = map(RZ, RZ, transpose M2)
  f2 L2
  L2 = f2 L2
  F2 = f2 F2

  use RZ
  sing1 = trim saturate(ideal F1 + ideal jacobian F1);
  sing2 = trim saturate(ideal F2 + ideal jacobian F2);
  factor leadCoefficient sing1_0
  F1e = sub(F1, e => 0)
  F2e = sub(F2, e => 0)
  SZ = ZZ[a,b,c,d]
  SQ = QQ (monoid SZ)
  F1e = sub(F1e, SZ)
  F2e = sub(F2e, SZ)
  factor det hessian F1e
  factor det hessian F2e
  polynomialContent det hessian F1e
  polynomialContent det hessian F2e
  (A1, phi1) = genericLinearMap SQ
  J = ideal last coefficients(phi1 (sub(F1e, target phi1)) - sub(F2e, target phi1))
  --gbTrace=3
  --gb J;
  sing1 = saturate (ideal F1e + ideal jacobian F1e)
  sing2 = saturate (ideal F2e + ideal jacobian F2e)
  factor leadCoefficient sing1_0
  factor leadCoefficient sing2_0
  sing1 : (sing1 : 4243)
  sing2 : (sing2 : 4243)

  sing1 : (sing1 : 11)
  sing2 : (sing2 : 11)

  sing1 : (sing1 : 31)
  sing2 : (sing2 : 31)

  Sp = (ZZ/31) (monoid SZ)
  (Ap, phip) = genericLinearMap Sp
  sub(F1e, target phip)
  J = ideal last coefficients(phip (sub(F1e, target phip)) - sub(F2e, target phip))
  see J  
  J1 = trim sub(  sing1 : (sing1 : 31), Sp)
  J2 = trim sub(  sing2 : (sing2 : 31), Sp)
  J' = ideal last coefficients(phip (gens sub(J1, target phip)) % gens sub(J2, target phip))
  J = J + J'
  J = sub(J, coefficientRing ring J)
  elapsedTime gbJ = groebnerBasis J; -- 23 sec
  see ideal gbJ
  comps = decompose ideal gbJ -- 8 components, but 4 don't have a point over ZZ/89.  The other 4 have a unique point.
  A1 = sub(Ap % comps_4, ZZ)
  A1 = 1/6 * A1
  A1 = sub(A1, ZZ) 
  fp = map(Sp, Sp, transpose sub(A1, coefficientRing Sp) )
  fp sub(F1e, Sp) - sub(F2e, Sp)
  fZZ = map(SZ, SZ, transpose A1)
  fZZ F1e - F2e -- 0

  -- let's try map that leaves e fixed:
  use RZ
  fR = map(RZ, RZ, {-b-c, c, -a, d, e})
  fR F1 - F2
  
  U = source phi
  T = coefficientRing U
  use T; use U
  A = (A1 | transpose matrix{{t_(1,1), t_(1,2), t_(1,3), t_(1,4)}} ) || matrix{{0,0,0,0,1}}
  fU = map(U, U, transpose A)
  trim ideal last coefficients(fU sub(F1, U) - sub(F2, U))
  fU = map(U, U, fU.matrix % oo)
  fU sub(F1, U) - sub(F2, U)
    
  A1 = sub(A1, coefficientRing source phi)
  A2 = (transpose A1 | transpose matrix{{t_(1,1), t_(1,2), t_(1,3), t_(1,4)}}) || matrix{{0,0,0,0,1}}
  psi = map(source phi, source phi, transpose A2)
  use target psi
  psi (e)    
  psi F1
  psi sub(F1, source psi) - sub(F2, source psi)
///

-- In this file, we work on code for finding n x n integer invertible matrices A
-- such that phi L1 = L2, phi F1 = F2,
-- where: phi (x) = i-th row of (transpose A)*x (CHECK).
-- and L1, L2 are c2-forms for CY3's X1, X2
-- and F1, F2 are cubic forms for X1, X2.

findEquivalenceIdeals -- might come up with several ideals, any which could work.
findEquivalence -- uses findEquivalenceIdeals, checks to find integer solutions for each component.
  -- returns only one of these if found.
  -- if cannot solve, return list of ideals, ones that cannot be solved.
  -- if no solutions at all, return null.
  -- if find a solution, return the matrix.
  
partitionCY3sByEquivalence
  -- input: a list of CY3 labels, and a hash table Xs
  -- output: a list of lists
  -- each of the sublists has the form:
  -- {{label, {label, matrix}, ...}
  -- also returns the unknown ones?
  
  -- e.g.: if input is a list {lab1, lab2}
  -- it calls findEquivalence on these two CY's
  -- if null: returns {{lab1}, {lab2}} -- these are distinct topologies
  -- if a matrix A0: returns {{lab1, {lab2, A0}}
  -- if an ideal: returns {{{lab1, lab2}}} -- uughh... indicating that these may be same.

isEquivalent = method()
isEquivalent(Sequence, Sequence, Matrix) := Boolean => (LF1, LF2, A) -> (
    (L1,F1) := LF1;
    (L2,F2) := LF2;
    R := ring L1;
    if R =!= ring F1 or R =!= ring L2 or R =!= ring F2 then 
        error "excepted c2 and cubic forms to be in the same ring";
    phi := map(R, R, A);
    phi L1 == L2 and phi F1 == F2
    )
isEquivalent(CYToolsCY3, CYToolsCY3, Matrix) := 
isEquivalent(CalabiYauInToric, CalabiYauInToric, Matrix) := Boolean => (X1, X2, A) -> (
    -- X1, X2 are CalabiYauInToric's (of the same h11 = h11(X1) = h11(X2)).
    -- A is an h11 x h11 matrix over ZZ, with determinant 1 or -1.
    -- if the topological data of X1, X2 are equivalent via A, then true is returned.
    if hh^(1,2) X1 =!= hh^(1,2) X2 then return false;
    if hh^(1,1) X1 =!= hh^(1,1) X2 then return false;
    LF1 := (c2Form X1, cubicForm X1);
    LF2 := (c2Form X2, cubicForm X2);
    isEquivalent(LF1, LF2, A)
    )

-- not: mapIsIsomorphism is essentially identical (but works with CYToolsCY3, CalabiYauInToric

genericLinearMap = method(Options => {Variable => null})
genericLinearMap Ring := opts -> R -> (
    -- R should be a polynomial ring in n variables.
    n := numgens R;
    K := coefficientRing R;
    t := if opts.Variable === null then getSymbol "t" else opts.Variable;
    T := K[t_(1,1)..t_(n,n)];
    TR := T [gens R, Join => false];
    A := map(T^n,,transpose genericMatrix(T, T_0, n, n));
    phi := map(TR, TR, transpose A);
    (A, phi)
    )

-- linearEquationConstraints, linearEquationConstraintsIdeal: not quite what we want.

-- findMaps: keep this one?
-- findIsomorphism -- same, I think, or close.
-- determineIsomorphism

  getEquivalenceIdealHelper = (LF1, LF2, A, phi) -> (
      T := target phi;
      B := ring A;
      LF1' := LF1/(f -> sub(f, T));
      LF2' := LF2/(f -> sub(f, T));
      I0 := sub(ideal last coefficients(phi LF1'_0 - LF2'_0), B);
      A0 := A % I0;
      phi0 := map(T, T, transpose A0);
      trim(I0 + sub(ideal last coefficients (phi0 LF1'_1 - LF2'_1), B))
      )

  getEquivalenceIdeal = method()
  getEquivalenceIdeal(Thing, Thing, HashTable) := Sequence => (lab1, lab2, Xs) -> (
      X1 := Xs#lab1;
      X2 := Xs#lab2;
      LF1 := (c2Form X1, cubicForm X1);
      LF2 := (c2Form X2, cubicForm X2);
      RQ := QQ (monoid ring LF1_0);
      (A,phi) := genericLinearMap RQ;
      (getEquivalenceIdealHelper(LF1, LF2, A, phi), A)
      )

---------------------------------------------
-- examples for h11=4
end--
restart
load "../FindEquivalence.m2"

DB4 = "../Databases/cys-ntfe-h11-4.dbm"
RZ = ZZ[a,b,c,d]
RQ = QQ (monoid RZ);
(Qs, Xs) = readCYDatabase(DB4, Ring => RZ);
allXs = sort keys Xs;
(A, phi) = genericLinearMap RQ;

  h12 = 70
  thisXs = select(sort keys Xs, lab -> hh^(1,2) Xs#lab == h12)  
  elapsedTime thisInvSet = partition(lab -> invariantsAll(Xs#lab), thisXs); 
  reps = select(keys thisInvSet, k -> #thisInvSet#k > 1)
  netList select(values thisInvSet, k -> #k > 1)

  thisInvSet#(reps#1)

  equivalenceIdealsSingularLocus ((418,0), (434,0), Xs, phi)
  tally apply(sort keys Xs, lab -> hh^(1,2) Xs#lab)
  saturate (ideal F1 + ideal jacobian F1)
  sub(oo, RQ)
  decompose oo
  
 
  ideal last coefficients (phi gens sub(singF1Q, source phi) % sub(singF2Q, source phi))

(A, phi) = genericLinearMap RQ
equivalenceIdealMain((418,0),(434,0),Xs,(A,phi)) -- this is meant to be tacked on to others?

(L1,F1) = (c2Form Xs#(418,0), cubicForm Xs#(418,0))
(L2,F2) = (c2Form Xs#(434,0), cubicForm Xs#(434,0))

gens gb equivalenceIdeals({L1, F1, ideal(b,c,d)}, {L2, F2, ideal(d,c,a-b)}, phi)


S1 = trim sub(ideal F1 + ideal jacobian F1, RQ)
S2 = trim sub(ideal F2 + ideal jacobian F2, RQ)

gens gb equivalenceIdeals({L1, F1, ideal(b,c,d)}, {L2, F2, ideal(d,c,a-b)}, phi)

decompose ideal gens gb equivalenceIdeals({S1, L1, F1}, {S2, L2, F2}, phi)
for j in oo list A % j
oo/det
A % ideal oo

-- These are the same, which we can determine with GV's
(X1, X2) = ((422, 0), (431, 0))/(lab -> Xs#lab)
(L1,F1) = (c2Form X1, cubicForm X1)
(L2,F2) = (c2Form X2, cubicForm X2)
--elapsedTime ideal gens gb equivalenceIdeals({L1, F1}, {L2, F2}, phi);
gv1 = partitionGVConeByGV(X1, DegreeLimit => 15)
gv2 = partitionGVConeByGV(X2, DegreeLimit => 15)
findLinearMaps(gv1, gv2)
phi0 = map(RZ, RZ, sub(first oo, ZZ))
phi0 L1 - L2
phi0 F1 - F2

equivalenceIdealSingularSet(X1, X2, phi)

-- Hessians --
  these = thisInvSet#(reps#2) -- {(425, 6), (430, 11)}
(A,phi) = genericLinearMap RQ  
(X1, X2) = these/(lab -> Xs#lab)//toSequence
elapsedTime equivalenceIdealsHessian(these#0, these#1, Xs, (A, phi));
elapsedTime equivalenceByHessian(these#0, these#1, Xs, (A, phi))
  these = thisInvSet#(reps#3) -- {(425, 6), (430, 11)}
  reps#3
  ours = take(these, 2)
  (list1, list2p) = singularLocusMatches(ours#0, ours#1, Xs, (A, phi))
  for J in flatten for p in list2p list decompose trim equivalenceIdeal(list1, p, phi) list A % J
  equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))

reps#0 -- use sing locus

for ours in subsets(thisInvSet#(reps#0), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))

for ours in subsets(thisInvSet#(reps#1), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- neither works yet

for ours in subsets(thisInvSet#(reps#1), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- neither works yet

for ours in subsets(thisInvSet#(reps#2), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi))
for ours in subsets(thisInvSet#(reps#3), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- works, but should use a better algorithm than just trying all pairs?

reps#4 
#thisInvSet#(reps#4)
for ours in subsets(thisInvSet#(reps#4), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) 
for ours in subsets(thisInvSet#(reps#4), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi)) 

reps#5 -- neither should work here
#thisInvSet#(reps#5)
for ours in subsets(thisInvSet#(reps#5), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) 
for ours in subsets(thisInvSet#(reps#5), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi)) 

rep = reps#6  -- this one is subtle: the ideals are nontrivial.  BUG? Did I do decompose above?
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) 
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi)) 

rep = reps#7  
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) 
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))  -- shows they are different.

rep = reps#8
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- shows there are 2 classes of 2 each...  Takes longer than I would prefer.
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))  -- 

rep = reps#9 -- neither works
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- 
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))  -- 

rep = reps#10 -- all 4 are equiv, both methods work
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- 
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))  -- 

rep = reps#11 -- all 2 are equiv, both work
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- 
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))  -- 

rep = reps#12
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- gives ideal, using decompose would work here
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))  -- works

rep = reps#13
thisset = thisInvSet#(rep) -- 28 in this set
for i from 1 to 27 list (ans := equivalenceByHessian(thisset#0, thisset#i, Xs, (A, phi)); print ans; ans) -- fast, all equiv to first one.

rep = reps#14
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- gives ideal, using decompose would work here
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))  -- works

rep = reps#14 -- one that goes to ellipse.
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi))
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi)) 

rep = reps#15 -- hessian and sing locus do not work here
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi))
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi)) 

rep = reps#16 -- hessian and sing locus do not work here
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi))
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi)) 

rep = reps#17 -- hessian and sing locus do not work here
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi))
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi)) 

rep = reps#18 -- hessian and sing locus do not work here
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi))
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi)) 

rep = reps#19
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- nope
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))  -- works

rep = reps#20 -- this one is interesting: but it seems like adding in det + 1, det - 1 is a good plan.
#thisInvSet#(rep)
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByHessian(ours#0, ours#1, Xs, (A, phi)) -- 
for ours in subsets(thisInvSet#(rep), 2) list equivalenceBySingularLocus(ours#0, ours#1, Xs, (A, phi))  -- 

-- all the same
for ours in subsets(thisInvSet#(rep), 2) list equivalenceByGVCone(ours#0, ours#1, Xs, (A, phi)) -- 

