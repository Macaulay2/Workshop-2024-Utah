-- Code to compute the topology of a CY 3-fold.
-- Also contains code to help determine whether two CY3's are homeomorphic
-- (this appears to be equivalent to being diffeomorphic).

topologicalData = method()
-- topologicalData CalabiYauInToric := TopologicalDataOfCY3 => X -> (
--     -- TODO: this does not consider torision in H_2(X, ZZ) or H_3(X, ZZ)
--     << "calling topologicalData" << endl;
--     elapsedTime new TopologicalDataOfCY3 from {
--         "h11" => hh^(1,1) cyPolytope X,
--         "h21" => hh^(2,1) cyPolytope X,
--         "c2" => c2 X,
--         "intersection numbers" => intersectionNumbers X
--         }
--     )

topologicalData CalabiYauInToric := TopologicalDataOfCY3 => X -> (
    -- TODO: this does not consider torision in H_2(X, ZZ) or H_3(X, ZZ)
    new TopologicalDataOfCY3 from {
        c2Form X,
        cubicForm X,
        hh^(1,1) X,
        hh^(1,2) X
        }
    )
    -- << "calling topologicalData" << endl;
    -- elapsedTime new TopologicalDataOfCY3 from {
    --     "h11" => hh^(1,1) cyPolytope X,
    --     "h21" => hh^(2,1) cyPolytope X,
    --     "c2" => c2 X,
    --     "intersection numbers" => intersectionNumbers X
    --     }
    -- )


-- this is the older, alternate version of this function.
-- this is to be removed.
-- topologicalData(CalabiYauInToric, Ring) := TopologicalDataOfCY3 => (X, RZ) -> (
--     V := ambient X;
--     Q := X#"polytope data";
--     P := polytope Q;
--     data := elapsedTime topologyOfCY3(V, basisIndices X);
--     -- this data above computes intersection numbers for all toric divisors. 
--     -- So we consider only the ones whose indices are contained in basis indices:
--     new TopologicalDataOfCY3 from {
--         "h11" => elapsedTime hh^(1,1) Q,
--         "h21" => elapsedTime hh^(2,1) Q,
--         "c2" => sub(data_3, vars RZ),
--         "cubic intersection form" => sub(data_2, vars RZ)
--         }
--     )

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

------------------------------------
-- Separating a set of topologies --
------------------------------------
-- This takes a set of labels (and matrices defining equivalences)
-- "Sets" field: is a list of buckets.  Two topologies in different buckets are definitely not the same.
-- Each bucket is a list of lists.
--   Two topologies in one element of this list are definitely the same.
--   Two topologies in different elements of this list are possibly the same possibly different.
-- two types of routines:
--  separate: this will only create new number of buckets (by separating buckets which are there.
--  combine: this will only coalesce sets in any given bucket.
--     (adding in matrices that show this).
--     note: TODO: if two sets are combined, we need to multiply all the matrices of one set by the new change of basis matrix.
TopologySet = new Type of MutableHashTable

-- TODO XXX: T#"Sets"#i is a list of {lab, {lab2, map2}, {lab3,map3}, ...}
--           mapj is a matrix over the integers of size h11 x h11.
topologySet = method()
topologySet(List, HashTable) := TopologySet => (labels, Xs) -> (
    new TopologySet from {
        "Sets" => {for k in labels list {k}},
        "CYHash" => Xs
        }
    )

info TopologySet := T -> (
    << "Total number of objects considered:         " << T#"Sets"/(x -> (x/length//sum))//sum << endl;
    << "Number of known different topologies:       " << #T#"Sets" << endl;
    << "Maximum possible # of different topologies: " <<  T#"Sets"/(x -> length x)//sum << endl;
    << "Largest number in one set:                  " << T#"Sets"/(x -> (x/length//max))//max << endl;
    t := T#"Sets"/length//tally;
    for i in keys t do (
        if i == 1 then 
            << "  Number of buckets with 1 class           " << t#1 << endl
        else 
            << "  Number of buckets with "|i|" classes         " << t#i << endl;
            );
    -- bigones := for i in keys t list if i >= 10 then t#i else continue;
    -- if #bigones > 0 then
    --   << "  Number of buckets with >= 10 classes     " << sum bigones << endl;
    )

representatives = method(Options => {IgnoreSingles => true})
representatives TopologySet := opts -> T -> (
    -- NOTE: ignores those with only one class in a bucket!
    sort for t1 in T#"Sets" list if opts.IgnoreSingles and #t1 == 1 then continue else t1/first//sort
    )
equivalences = method(Options => {IgnoreSingles => true})
equivalences TopologySet := opts -> T -> (
    hashTable flatten for t1 in T#"Sets" list for t2 in t1 list (
        if #t2 == 1 and opts.IgnoreSingles then continue;
        (first t2) => drop(t2, 1)
        )
    )

separateIfDifferent = method()
separateIfDifferent(TopologySet, Function) := (T, fun) -> (
    -- fun takes a CalabiYauInToric, and returns some value.
    -- if fun X1 =!= fun X2, then X1 and X2 are distinct topologies.
    -- if the same, then nothing is asserted.
    Xs := T#"CYHash";
    -- the 'flatten' on the next line makes one list of all definitely distinct topologies
    newsets := flatten for L in T#"Sets" list (
        -- L is a list of labels, all are equivalent (so maybe only one, but always >= 1).
        if #L === 1 then {L}
        else (
            -- TODO: is this correct? XXX Just changed, not fixed....
            P := partition(lab -> fun Xs#(first lab), L);
            values P
            )
        );
    new TopologySet from {
        "Sets" => newsets,
        "CYHash" => T#"CYHash"
        }
    )

combineIfSame = method()
combineIfSame(TopologySet, Function) := (T, fun) -> (
    -- fun takes a CalabiYauInToric, and returns some value.
    -- if fun X1 === fun X2, then X1 and X2 are the same topology
    -- if different, then nothing is asserted.
    Xs := T#"CYHash";
    -- the 'flatten' on the next line makes one list of all definitely distinct topologies
    newsets := for L in T#"Sets" list (
        -- L is a list of labels, all are equivalent (so maybe only one, but always >= 1).
        P := partition(lab -> fun Xs#(first lab), L);
        (values P)/flatten
        -- TODO XXX: the previous line should add in identity maps
        );
    new TopologySet from {
        "Sets" => newsets,
        "CYHash" => T#"CYHash"
        }
    )

-- REMOVE THIS ONE: use combineByGV...
separateByGV = method(Options => {DegreeLimit => 15})
-- separateByGV TopologySet := opts -> T -> (
--     Xs := T#"CYHash";
--     newSets := for Ls in T#"Sets" list (
--         L1s := Ls/first; -- these are the ones we want to split up
--         Indices := hashTable for i from 0 to #Ls-1 list first Ls#i => i;
--         << "----------------" << endl;
--         print L1s;
--         print Indices;
--         << "----------------" << endl;
--         P := partitionByTopology(L1s, Xs, opts.DegreeLimit);
--         newlist := for k in keys P list {k}|(P#k);
--         newlist
--         -- what is the best way to get the ones that are the same into the same set?
--         );
--     new TopologySet from {
--         "CYHash" => Xs,
--         "Sets" => newSets
--         }
--     )

combineBucketByGV = method(Options => {DegreeLimit => 15})
combineBucketByGV(List, HashTable) := opts -> (Ls, Xs) -> (
    if #Ls === 1 then Ls
    else (
        L1s := Ls/first; -- these are the ones we want to split up
        L1rest := hashTable for L in Ls list (
            {first L, drop(L, 1)}
            );
        -- print L1s;
        -- << "----------------" << endl;
        P := partitionByTopology(L1s, Xs, opts.DegreeLimit);
        newlist := for k in keys P list (
            {k} | P#k | L1rest#k | flatten for x in P#k list L1rest#(first x)
            );
--        if #(keys P) > 1 then error "debug me";
        newlist
    ))

combineByGV = method(Options => {DegreeLimit => 15})
combineByGV TopologySet := opts -> T -> (
    Xs := T#"CYHash";
    newSets := for Ls in T#"Sets" list (
        combineBucketByGV(Ls, Xs, DegreeLimit => opts.DegreeLimit)
        );
    new TopologySet from {
        "CYHash" => Xs,
        "Sets" => newSets
        }
    )


findMapsFROMBELOW = (top1, top2, A, phi, RQ) -> (
    TR := target phi;
    n := numgens TR;
    T := coefficientRing TR;
    toTR := f -> sub(f, TR);
    toQQ := f -> sub(f, RQ);
    evalphi := (F,G) -> trim sub(ideal last coefficients((phi toTR F) - toTR G), T);
    (L1, F1, h11, h12) := toSequence top1;
    (L2, F2, l11, l12) := toSequence top2;    
    RZ := ring L1;
    if RZ =!= ring F1 or RZ =!= ring L2 or RZ =!= ring F2 then error "expected polynomials over the same ring";
    if h11 != l11 or h12 != l12 then return null;
    I := (evalphi(L1, L2) + evalphi(F1, F2));
    if I == 1 then return null;
    -- first see if there is a unique solution.
    -- if codim I === n*n and degree I === 1 then (
    --     A0 := A % I;
    --     if support A0 === {} then (
    --         A0 = lift(A0, QQ);
    --         phi0 := map(RQ, RQ, transpose A0);
    --         if phi0 toQQ L1 != toQQ L2 or phi0 toQQ F1 != toQQ F2 then error "map is not correct!";
    --         return (A0, phi0)
    --         );
    --     );
    -- now let's look through all of the components for a smooth point.
    compsI := decompose I;
    As := for c in compsI list A % c;
    As = for a in As list try lift(a, ZZ) else continue; -- grab the ones that lift.
    As = select(As, a -> (d := det a; d == 1 or d == -1));
    if #As > 0 then (
        A0 := As#0;
        phi0 := map(RZ, RZ, transpose A0);
        if phi0 L1 != L2 or phi0 F1 != F2 then error "map is not correct!";
        (A0, phi0)
        )
    else (
        if any(compsI, c -> codim c < n*n or degree c =!= 1) then (
            << "warning: there might be a map in this case!" << endl;
            << netList compsI << endl;
            << "----------------------------------" << endl;
            compsI 
            )
        else null
        )
    )

separateAndCombineViaAnsatz = method()
separateAndCombineViaAnsatz(TopologySet, HashTable, Ring) := TopologySet => (tops, Ts, RQ) -> (
    )
separateAndCombineViaAnsatz(List, HashTable, Ring) := List => (Ls, Ts, RQ) -> (
    -- Ls: is a list of labels.
    -- Ts: is a hash table containing for each lab in labels, the (c2, cubic, h11, h12) of Xs#lab.
    -- restL#lab is a list of either pairs: {lab, matrix}, or of simply: labels.
    -- RQ is QQ[h11 variables].
    -- Result: list of distinct topologies found, and in each list: each element is a list of equivalent topologies.
    (A, phi) := genericLinearMap RQ;
    if #Ls === 1 then return Ls;
    distinctTops := new MutableHashTable; -- label => list of {label, matrix}, those with the same topology
    for i in Ls do (
        Ti := Ts#i;
        prev := keys distinctTops;
        isFound := false;
        for j in prev do (
            Tj := Ts#j;
            ans := findMaps(Tj, Ti, A, phi, RQ); -- either a matrix of integers or null, or "unknown"
            if ans === null then (
                -- Ti is distinct from Tj
                )
            else if class first ans === Matrix then (
                -- we have a match!
                (A0, phi0) := ans;
                -- The following seems to be a redundant check!
                if all(flatten entries A0, a -> liftable(a, ZZ))
                then (
                    isFound = true;
                    distinctTops#j = append(distinctTops#j, {i, lift(A0, ZZ)});
                    break;
                    )
                )
            else (
                << (i,j) << " might be the same, might not *** " << endl;
                )
            );
        if not isFound then (
            distinctTops#i = {};
            << "found new top: " << i << endl;
            );
        );
    -- take distinctTops, so something with them....
    )

-- partitionH113sByTopology(List, HashTable, Ring) := HashTable => (Ls, Ts, RQ) -> (
--     -- Ls is a list of labels to separate.
--     -- Ts is a hash table of label => {c2, cubicform, h11, h12}
--     -- This function first separates these by the invariants: invariantsAll.
--     -- The for each pair in each set, it attempts to find a map between them.
--     -- output: a hashtable, keys are labels, values are lists of {label, matrix}
--     (A, phi) := genericLinearMap RQ;
--     if #Ls === 1 then return hashTable {Ls#0 => {}};
--     H := partition(lab -> invariantsAll toSequence Ts#lab, Ls);
--     distinctTops := new MutableHashTable; -- label => list of {label, matrix}, those with the same topology
--     for i in Ls do (
--         Ti := Ts#i;
--         -- now we attempt to match this with each key of distinctTops
--         << "trying " << i << endl;
--         prev := keys distinctTops;
--         isFound := false;
--         for j in prev do (
--             Tj := Ts#j;
--             ans := findMaps(Tj, Ti, A, phi, RQ);
--             --if j == (115,0) and i == (120,0) then error "debug me";
--             if ans === null then (
--                 -- Ti is distinct from Tj
--                 )
--             else if class first ans === Matrix then (
--                 -- we have a match!
--                 (A0, phi0) := ans;
--                 if all(flatten entries A0, a -> liftable(a, ZZ))
--                 then (
--                     isFound = true;
--                     distinctTops#j = append(distinctTops#j, {i, lift(A0, ZZ)});
--                     break;
--                     )
--                 )
--             else (
--                 << (i,j) << " might be the same, might not" << endl;
--                 )
--             );
--         if not isFound then (
--             distinctTops#i = {};
--             << "found new top: " << i << endl;
--             );
--         );
--     new HashTable from distinctTops
--     )


-- REMOVE?  This is the start of a union-find algorithm.  But we are not using it, I think.
-- combineSet = (Ls, binfun) -> (
--     -- binfun(X1,X2) should return a matrix if these are the same topology, null if we don't know.
--     -- note: the number of newsets is identical.  We might just be coalescing elements in one set.
--     -- NOTE: currently the matrix is lost, we need to be able to keep it!
--     upnode := new MutableList from 0..#Ls-1;
--     nodesize := new MutableList from (#Ls : 1);
--     find := x -> (
--         root := x;
--         while upnode#root != root do root = upnode#root;
--         while upnode#x != root do (
--             par := upnode#x;
--             upnode#x = root; 
--             x = par;
--             );
--         root
--         );
--     union := (x,y) -> (
--         x = find x;
--         y = find y;
--         if x === y then return;
--         if nodesize#x  < nodesize#y then (
--             (x,y) = (y,x);
--             );
--         upnode#y = x;
--         nodesize#x = nodesize#x + nodesize#y;
--         nodesize#y = 0;
--         );
--     for i from 0 to #Ls - 2 do 
--       for j from i+1 to #Ls - 1 do (
--           aij := binfun(Ls#i,Ls#j);
--           if aij =!= null then union(i,j);
--           );
--     P := partition(x -> find x, toList(0..#Ls-1));
--     for p in values P list (for p1 in p list Ls#p1)
--     )

-- doit = (T) -> (
--     Xs := T#"CYHash";
--     combineSet((T#"Sets"#0), (lab1,lab2) -> (
--             X1 := Xs#(first lab1);
--             X2 := Xs#(first lab2);
--             c2Form X1 == c2Form X2 and cubicForm X1 == cubicForm X2
--             ))
--     )


///
  -- Analyze h11=3 examples
restart
debug needsPackage "StringTorics"
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ);
  --(Qs, Xs) = readCYDatabase("../m2-examples/foo-cys-ntfe-h11-3.dbm", Ring => RZ);
  (Qs, Xs) = readCYDatabase("./Databases/cys-ntfe-h11-3.dbm", Ring => RZ);
  -- First, let's only consider those with torsion free class group.
   torsions = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k != ZZ^3;
      if istor then k else continue
      )

  allXs = sort select(keys Xs, x -> not member(x#0, torsions))
  allXs = sort keys Xs
  allT = topologySet(allXs, Xs);

  allT1 = combineIfSame(allT, X -> (c2Form X, cubicForm X))
  info allT1 

  elapsedTime allT = separateIfDifferent(allT, invariantsAll) -- 17 sec
  info allT

  elapsedTime allT3 = separateByGV allT2 -- 44 sec
  info allT3

  onestocheck = for x in allT3#"Sets" list if #x == 1 then continue else (
      x/first
      )
  
  flatten for x in onestocheck list (
      flatten for y in subsets(x, 2) list (
        print y;
        ans := getEquivalenceIdeal(y#0, y#1, Xs);
        print ans;
        ans
        )
      )
  
///  
  
///
restart
debug needsPackage "StringTorics"
  R = ZZ[a,b,c,d]
  RQ = QQ (monoid R);
  (Qs, Xs) = readCYDatabase("mike-ntfe-h11-4.dbm", Ring => R);

  allT = topologySet(sort keys Xs, Xs);
  #allT#"Sets" == 1
  #allT#"Sets"#0 == 1994
  info allT

  allT1 = combineIfSame(allT, X -> (c2Form X, cubicForm X))
  -- only 14 are the same as any other, and if two are the same, they also have same hh^(1,2).
  -- (checked this explicitly).
  info allT1 -- note very few are the same!
  netList allT1#"Sets"
  for x in allT1#"Sets" list (x/length)//tally

  elapsedTime allT2 = separateIfDifferent(allT1, invariantsAll) -- 270 sec
  elapsedTime allT3 = separateByGV allT2 -- 870 sec
  info allT3

  for x in allT3#"Sets" list if #x == 1 then continue else (
      x/first
      )
  onestocheck = {
      {(163, 2), (163, 3)}, 
      {(811, 0), (882, 3)}, 
      {(387, 1), (387, 3)}, 
      {(339, 1), (337, 1)}, 
      {(1001, 0), (982, 0)}, 
      {(436, 1), (433, 1)}, 
      {(364, 0), (344, 1)}, 
      {(1094, 2), (1090, 2)}, 
      {(246, 0), (249, 1)}, 
      {(1147, 0), (1146, 0)}, 
      {(436, 0), (433, 0)}, 
      {(831, 0), (806, 0)}, 
      {(1121, 0), (1122, 0)}, 
      {(1067, 0), (1076, 0)}, 
      {(930, 6), (927, 2)}, 
      {(981, 0), (997, 0)}, 
      {(364, 1), (377, 0)}, 
      {(455, 0), (451, 0)}, 
      {(403, 0), (408, 2)}, 
      {(1090, 0), (1094, 0)}, 
      {(714, 0), (707, 0), (695, 5), (709, 0)}, 
      {(991, 0), (998, 3)}, 
      {(1002, 1), (979, 0), (1005, 4)}, 
      {(716, 3), (705, 0)}, 
      {(851, 0), (814, 0), (810, 0)}, 
      {(1123, 0), (1124, 0)}, 
      {(478, 1), (473, 0)}, 
      {(1085, 0), (1082, 0), (1086, 0)}, 
      {(334, 1), (329, 1)}, 
      {(316, 0), (322, 0)}, 
      {(1004, 0), (994, 0)}, 
      {(265, 1), (250, 0), (254, 5), (228, 1)}, 
      {(647, 0), (645, 0)}, 
      {(529, 0), (510, 0)}, 
      {(935, 4), (915, 0)}, 
      {(628, 0), (653, 3)}, 
      {(943, 0), (958, 0)}, 
      {(329, 0), (334, 0)}, 
      {(319, 0), (321, 0), (302, 0)}, 
      {(875, 0), (854, 0)}, 
      {(1183, 2), (1182, 0)}, 
      {(449, 2), (419, 0)}, 
      {(383, 0), (356, 1)}, 
      {(250, 1), (213, 1)}, 
      {(938, 3), (933, 0)}, 
      {(1082, 5), (1077, 3)}, 
      {(668, 1), (649, 1), (626, 0)}, 
      {(80, 3), (72, 0), (80, 6)}, 
      {(616, 0), (606, 0)}, 
      {(540, 0), (545, 0)}, 
      {(331, 0), (337, 0), (339, 2)}, 
      {(397, 0), (350, 0)}, 
      {(900, 0), (901, 0)}, 
      {(930, 1), (927, 8)}, 
      {(559, 0), (577, 0)}, 
      {(938, 10), (933, 4)}, 
      {(931, 0), (932, 2)}, 
      {(884, 0), (820, 0)}, 
      {(650, 3), (656, 1), (630, 0)}, 
      {(976, 0), (989, 0)}, 
      {(552, 0), (551, 2)}, 
      {(924, 0), (937, 0)}}

  176 == # flatten for x in onestocheck list flatten for y in subsets(x,2) list y
  
  flatten for x in onestocheck list (
      flatten for y in subsets(x, 2) list (
        print y;
        ans := getEquivalenceIdeal(y#0, y#1, Xs);
        print ans;
        ans
        )
      )
  select(oo, i -> numgens i > 1)
  getEquivalenceIdeal((163, 2), (163, 3), Xs)
  getEquivalenceIdeal((1123,0),(1124,0), Xs)

  lab1 = (163,2)
  lab2 = (163,3)
  X1 = Xs#lab1;
  X2 = Xs#lab2;
  LF1 = (c2Form X1, cubicForm X1);
  LF2 = (c2Form X2, cubicForm X2);
  RQ = QQ (monoid ring LF1_0);
  (A,phi) = genericLinearMap RQ;
  T = target phi;
  B = ring A;
  LF1' = LF1/(f -> sub(f, T));
  LF2' = LF2/(f -> sub(f, T));
  I0 = sub(ideal last coefficients(phi LF1'_0 - LF2'_0), B);
  A0 = A % I0;
  phi0 = map(T, T, A0);
  I1 = sub(ideal last coefficients (phi0 LF1'_1 - LF2'_1), B)
  gens gb(I0 + I1)

  select(allT3#"Sets", x -> #x > 1)
  -- now we need to check that (hopefully) each set is not homeomorphic to any other set.
  -- It seems that now an ansatz might work?
  
  TbyH12s = separateIfDifferent(allT, X -> hh^(1,2) X);
  info TbyH12s
  #TbyH12s#"Sets"
  for x in TbyH12s#"Sets" list #x
  
  
  smallset = select(sort keys Xs, lab -> hh^(1,2) Xs#lab == 148) -- 
  T = topologySet(sort smallset, Xs)
  allTsame = combineIfSame(T, X -> (c2Form X, cubicForm X))
  info allTsame
  elapsedTime T1 = separateIfDifferent(allTsame, X -> invariantsAll X)
  info T1
  netList T1#"Sets"
  T2 = separateByGV T1
          
///

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

------------------------------------


factors = method()
factors RingElement := (F) -> (
     facs := factor F;
     facs//toList/toList/reverse
     )

invariants = method()
invariants List := (f) -> (
    RQ := QQ[gens ring first f];
    facs := select((factors f_1 )/toList/last, g -> support g != {});
    l1 := sub(f_0, RQ);
    f1 := sub(f_1, RQ);
    d := dim saturate ideal jacobian f1;
    nc := # decompose ideal(l1, f1);
    singZ := flatten entries gens gb saturate(ideal(f_1) + ideal jacobian f_1);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else first badp;
    {badp, (trim content f_0)_0, (trim content f_1)_0, #facs, d, nc, f_2, f_3}
    )



pointCount = method()
pointCount(RingElement, ZZ) := ZZ => (F, p) -> (
    -- F is a polynomial in 3 variables (FIXME: any number of variables)
    -- p is a prime number
    kk := ZZ/p;
    R := kk (monoid ring F);
    Fp := sub(F, R);
    allpts := allPoints(p, numgens ring F);
    allmaps := allpts/(pt -> map(kk, R, pt));
    ans1 := # for a in allpts list (phi := map(kk, R, a); if phi Fp == 0 then a else continue);
    --ans2 := # for x in (0,0,0)..(p-1,p-1,p-1) list if sub(Fp, matrix{{x}}) == 0 then x else continue;
    --if ans1 != ans2 then << "My previous code was incorrect" << endl;
    ans1
    )

allPointMaps = method()
allPointMaps(ZZ, Ring) := (p, R) -> (
    N := numgens R;
    K := ZZ/p;
    pts := allPoints(p,N);
    for a in pts list map(K, R, a)
    )



-- allPointMaps(ZZ, ZZ, Ring) := (p, n, R) -> (
--     N := numgens R;
--     K := GF(p,n);
--     pts := allPoints(p,N);
--     for a in pts list map(K, R, a)
--     )

-- pointCount = method()
-- pointCount(RingElement, RingElement, List) := ZZ => (L, F, pts) -> (
--     -- pts is a list of list of ring maps
--     --   each list is generally all of the points in k^N, for k = ZZ/p, p = 2,3,5,7,11,13, maybe k = GF 4, ...
    
--     -- F is a polynomial in 3 variables (FIXME: any number of variables)
--     -- p is a prime number
--     kk := ZZ/p;
--     R := kk (monoid ring F);
--     Fp := sub(F, R);
--     allpts := allPoints(p, numgens ring F);
--     allmaps := allpts/(pt -> map(kk, R, pt));
--     ans1 := # for a in allpts list (phi := map(kk, R, a); if phi Fp == 0 then a else continue);
--     --ans2 := # for x in (0,0,0)..(p-1,p-1,p-1) list if sub(Fp, matrix{{x}}) == 0 then x else continue;
--     --if ans1 != ans2 then << "My previous code was incorrect" << endl;
--     ans1
--     )

invariants CalabiYauInToric := List => X -> (
    L := c2Form X;
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    RZ := ring L;
    if ring F =!= RZ then 
        error "expected same rings for c2 and cubic form";
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    facs := select((factors F)/toList/last, g -> support g != {});
    LQ := sub(L, RQ);
    FQ := sub(F, RQ);
    d := dim saturate ideal jacobian FQ;
    nc := # decompose ideal(LQ, FQ);
    singZ := flatten entries gens gb saturate(ideal(F) + ideal jacobian F);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else sub(first badp, ZZ);
    ptcounts := for p in {2, 3, 5, 7, 11} list pointCount(F, p);
    {h11, h12, ptcounts, badp, (trim content L)_0, (trim content F)_0, #facs, d, nc}
    )

invariants0 = method()
invariants0 CalabiYauInToric := List => X -> (
    L := c2Form X;
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    contentL := (trim content L)_0;
    contentF := (trim content F)_0;
    {h11, h12} |  {contentL, contentF}
    )

invariants1 = method()
invariants1 CalabiYauInToric := List => X -> (
    L := c2Form X;
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    contentL := (trim content L)_0;
    contentF := (trim content F)_0;
    ptcounts := for p in {2, 3, 5, 7, 11} list pointCount(F, p);
    {h11, h12} |  ptcounts | {contentL, contentF}
    )

invariants2 = method()
invariants2 CalabiYauInToric := List => X -> (
    L := c2Form X;
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    RZ := ring L;
    if ring F =!= RZ then 
        error "expected same rings for c2 and cubic form";
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    facs := select((factors F)/toList/last, g -> support g != {});
    LQ := sub(L, RQ);
    FQ := sub(F, RQ);
    d := dim saturate ideal jacobian FQ;
    nc := # decompose ideal(LQ, FQ);
    singZ := flatten entries gens gb saturate(ideal(F) + ideal jacobian F);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else sub(first badp, ZZ);
    {h11, h12, badp, (trim content L)_0, (trim content F)_0, #facs, d, nc}
    )

invariants3 = method()
invariants3 CalabiYauInToric := List => X -> (
    -- these are some invariants only involving the cubic form, not the c2 form...
    -- (but currently also involving h11 and h12.
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    RZ := ring F;
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    facs := select((factors F)/toList/last, g -> support g != {});
    FQ := sub(F, RQ);
    d := dim saturate ideal jacobian FQ;
    singFZ := ideal gens gb saturate(ideal(F) + ideal jacobian F);
    sing := (F) -> ideal F + ideal jacobian F;
    linearcontent := (I) -> (
        if I == 0 then return 0;
        lins := select(I_*, f -> f != 0 and first degree f <= 1);
        if #lins == 0 then return 0;
        gcd for ell in lins list (trim content ell)_0
        );
    lincontent := linearcontent saturate sing F;
    singZ := flatten entries gens gb saturate(ideal(F) + ideal jacobian F);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else sub(first badp, ZZ);
    {h11, h12, badp, (trim content F)_0, #facs, d, lincontent}
    )

invariants4 = method()
invariants4 CalabiYauInToric := List => X -> (
    L := c2Form X;
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    RZ := ring L;
    if ring F =!= RZ then 
        error "expected same rings for c2 and cubic form";
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    facs := select((factors F)/toList/last, g -> support g != {});
    LQ := sub(L, RQ);
    FQ := sub(F, RQ);
    d := dim saturate ideal jacobian FQ;
    sing := (F) -> ideal F + ideal jacobian F;
    linearcontent := (I) -> (
        if I == 0 then return 0;
        lins := select(I_*, f -> f != 0 and first degree f <= 1);
        if #lins == 0 then return 0;
        gcd for ell in lins list (trim content ell)_0
        );
    lincontent := linearcontent saturate sing F;
    nc := # decompose ideal(LQ, FQ);
    singZ := flatten entries gens gb saturate(ideal(F) + ideal jacobian F);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else sub(first badp, ZZ);
    {h11, h12, badp, (trim content L)_0, (trim content F)_0, #facs, d, nc, lincontent}
    )

-- This one contains the best info we have to date (which isn't quite good enough).
polynomialContent = method()

-- Assumption: F is a polynomial over ZZ (not over a field).
-- Outout: the integer content of F.
polynomialContent RingElement := F -> (trim content F)_0

integerPart = method()
integerPart Ideal := (I) -> (
    -- expected: I is an ideal in a polynomial ring over ZZ.
    gs := select(flatten entries gens gb I, f -> support f === {});
    if #gs == 0 then 0 
    else if #gs == 1 then lift(gs_0, ZZ) 
    else error "internal error: somehow have two generators in ZZ in this GB"
    )

--factorShape = method()
-- This one is WRONG: lift(xxx, ZZ) could be positive or negative.  Those cannot be different.
-- factorShape RingElement := F -> (
--     facs := factors F;
--     sort for x in facs list if support x#1 == {} then 
--             {0, lift(x#1, ZZ)}
--         else
--             {first degree x#1, x#0}
--     )

factorShape = method()
factorShape RingElement := F -> (
    facs := factors F;
    sort for x in facs list if support x#1 == {} then 
            {0, abs lift(x#1, ZZ)}
        else
            {first degree x#1, x#0}
    )

-- factorShape RingElement := F -> (
--     facs := factors F;
--     con := trim content F;
--     con = con_0;
--     posfactors := sort for x in facs list if support x#1 == {} then 
--                      continue
--                    else
--                      {first degree x#1, x#0};
--     prepend({0, con},  posfactors)
--     )


invariantsAllX = method()
invariantsAllX(RingElement, RingElement, ZZ, ZZ) := (L, F, h11, h12) -> (
    RZ := ring L;
    if ring F =!= RZ then 
        error "expected same rings for c2 and cubic form";
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    toQQ := F -> sub(F, vars RQ);
    sing := (cod, I) -> trim(I + minors(cod, jacobian I));
    linearcontent := (I) -> (
        if I == 0 then return 0;
        lins := select(I_*, f -> f != 0 and first degree f <= 1);
        if #lins == 0 then return 0;
        gcd for ell in lins list (trim content ell)_0
        );
    FQ := toQQ F;
    LQ := toQQ L;
    inv0 := polynomialContent L;
    inv1 := polynomialContent F;
    -- dimension and degree of each component of the singular loci over QQ.
--    inv2 := sort for c in decompose sing_1 ideal FQ list {codim c, degree c};
--    inv3 := sort for c in decompose sing_2 ideal(LQ, FQ) list {codim c, degree c};
--    inv4 := betti res saturate sing_1 ideal FQ;
    -- inverse system of FQ
--    inv5 := betti res inverseSystem FQ; -- not clear this one is worthwhile
    -- integer parts of singular loci.
    conductF := integerPart saturate sing_1 ideal F;
    conductLF := integerPart saturate sing_2 ideal(L,F);
    inv6 := conductF;
    inv7 := conductLF;
    inv8 := factorShape det hessian F;
    inv9 := linearcontent saturate sing_1 F;
    hashTable {"h11" => h11,
     "h12" => h12,
     "c(L)" => inv0, 
     "c(F)" => inv1, 
--     "comps sing FQ" => inv2, 
--     "comps sing LFQ" => inv3, 
--     "bettti sing LFQ" => inv4,
--     "betti inv F" => inv5,
     "conduct(F)" => inv6,
     "conduct(L,F)}" => inv7,
     "hessian shape" => inv8,
     "lincontent sing F" => inv9
     }
    )

invariantsAll = method()
invariantsAll(RingElement, RingElement, ZZ, ZZ) := (L, F, h11, h12) -> (
    RZ := ring L;
    if ring F =!= RZ then 
        error "expected same rings for c2 and cubic form";
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    toQQ := F -> sub(F, vars RQ);
    sing := (cod, I) -> trim(I + minors(cod, jacobian I));
    linearcontent := (I) -> (
        if I == 0 then return 0;
        lins := select(I_*, f -> f != 0 and first degree f <= 1);
        if #lins == 0 then return 0;
        gcd for ell in lins list (trim content ell)_0
        );
    FQ := toQQ F;
    LQ := toQQ L;
    inv0 := polynomialContent L;
    inv1 := polynomialContent F;
    -- dimension and degree of each component of the singular loci over QQ.
    inv2 := sort for c in decompose sing_1 ideal FQ list {codim c, degree c};
    inv3 := sort for c in decompose sing_2 ideal(LQ, FQ) list {codim c, degree c};
    inv4 := betti res saturate sing_1 ideal FQ;
    -- inverse system of FQ
--    inv5 := betti res inverseSystem FQ; -- not clear this one is worthwhile
    -- integer parts of singular loci.
    conductF := integerPart saturate sing_1 ideal F;
    conductLF := integerPart saturate sing_2 ideal(L,F);
    inv6 := conductF;
    inv7 := conductLF;
    inv8 := factorShape det hessian F;
    inv9 := linearcontent saturate sing_1 F;
    hashTable {"h11" => h11,
     "h12" => h12,
     "c(L)" => inv0, 
     "c(F)" => inv1, 
     "comps sing FQ" => inv2, 
     "comps sing LFQ" => inv3, 
     "bettti sing LFQ" => inv4,
--     "betti inv F" => inv5,
     "conduct(F)" => inv6,
     "conduct(L,F)}" => inv7,
     "hessian shape" => inv8,
     "lincontent sing F" => inv9
     }
    )

invariantsAll CalabiYauInToric := X -> invariantsAll(c2Form X, cubicForm X, hh^(1,1) X, hh^(1,2) X)

mapIsIsomorphism = method()
mapIsIsomorphism(Matrix, CalabiYauInToric, CalabiYauInToric) :=
mapIsIsomorphism(Matrix, CYToolsCY3, CYToolsCY3) := Boolean => (M, X1, X2) -> (
    -- M is a matrix over the base field, a possible map giving
    -- an isomorphism of topologies.
    -- T1, T2 are two topologies.
    F1 := cubicForm X1;
    F2 := cubicForm X2;
    RZ := ring F1;
    if RZ =!= ring F2 then error "expected the same picard ring";
    phi := map(RZ, RZ, M);
    (phi c2Form X1 == c2Form X2) and (phi F1 == F2)
    )

partitionByTopology = method()
partitionByTopology List := LGVs -> (
    -- LGVs is a list of X => gvPartition.
    -- output: a hashtable, keys are labels, values are lists of {label, matrix}
    labels := for x in LGVs list label first x;
    hashXs := hashTable for y in LGVs list (label first y) => y;
    distinctTops := new MutableHashTable; -- label => list of labels.
    for i in labels do (
        Xi := first hashXs#i;
        GVi := last hashXs#i;
        << "trying " << i << endl;
        prev := keys distinctTops;
        isFound := false;
        for j in prev do (
            Xj := first hashXs#j;
            GVj := last hashXs#j;
            -- compare CY's i, j
            Ms := findLinearMaps(GVj, GVi);
            if Ms === {} then continue;
            Ms = Ms/(m -> lift(m, ZZ));
            Ms = select(Ms, m -> (d := det m; d === 1 or d === -1));
            isIsos := Ms/(m -> mapIsIsomorphism(m, Xj, Xi));
            if any(isIsos, x -> true) then (
                mi := position(isIsos, x -> true);
                << "found isomorphism" << endl;
                distinctTops#j = append(distinctTops#j, {i, Ms#mi});
                isFound = true;
                break;
                ));
        if not isFound then (
            distinctTops#i = {};
            << "found new top: " << topologicalData Xi << endl;
            );
        );
    new HashTable from distinctTops
    )

partitionByTopology(List, HashTable, ZZ) := (Ls, Xs, degreelimit) -> (
    -- output: a hashtable, keys are labels, values are lists of {label, matrix}
    labels := Ls;
    if #Ls === 1 then return hashTable {Ls#0 => {}};
    hashXs := Xs;
    GVs := hashTable for lab in labels list lab => partitionGVConeByGV(hashXs#lab, DegreeLimit => degreelimit);
    distinctTops := new MutableHashTable; -- invariants => list of labels.
    for i in labels do (
        Xi := hashXs#i;
        GVi := GVs#i;
        << "trying " << i << endl;
        prev := keys distinctTops;
        isFound := false;
        if GVi === null then prev = {}; -- we cannot use GV with non-favorables currently.
        for j in prev do (
            Xj := hashXs#j;
            GVj := GVs#j;
            if GVj === null then continue;  -- we cannot use GV with non-favorables currently.
            -- compare CY's i, j
            Ms := findLinearMaps(GVj, GVi);
            if Ms === {} then continue;
            --Ms = Ms/(m -> lift(m, ZZ));
            Ms = for m in Ms list try lift(m, ZZ) else continue;
            Ms = select(Ms, m -> (d := det m; d === 1 or d === -1));
            isIsos := Ms/(m -> mapIsIsomorphism(m, Xj, Xi));
            if any(isIsos, x -> x == true) then (
                mi := position(isIsos, x -> x == true);
                << "found isomorphism from " << j << " to " << i << ": " << Ms#mi << endl;
                distinctTops#j = append(distinctTops#j, {i, Ms#mi});
                isFound = true;
                break;
                ));
        if not isFound then (
            distinctTops#i = {};
            << "found new top: " << i << endl;
            );
        );
    new HashTable from distinctTops
    )

linearEquationConstraints = method()
linearEquationConstraints(Matrix, RingMap, List, List) := Sequence => (A, phi, Ls, pts) -> (
    -- each entry of Ls is a list/sequence of length 2: {F, G}
    -- where F, G are polynomials in a ring R, (A, phi) are obtained from
    -- genericLinearMap.  We return the ideal of constraints in T
    -- for which phi(F) = G, for all pairs {F,G} in Ls.
    -- We also return A0, phi0 corresponding to these constraints.
    T := ring A;
    A0 := A;
    TR := target phi;
    I := sum for L in Ls list ideal sub(last coefficients(phi L_0 - L_1), T);
    if I != 0 then A0 = sub(A, T) % (trim I);
    J := sum for pq in pts list minors(2, 
        (A0 * (transpose matrix {pq#1})) | transpose matrix {pq#0});
    J1 := if J != 0 then I+J else I;
    if J != 0 then A0 = A0 % (trim J1);
    (A0, map(TR, TR, transpose A0), trim ideal gens gb J1)
    )

linearEquationConstraintsIdeal = method()
linearEquationConstraintsIdeal(Matrix, RingMap, List, List) := Sequence => (A, phi, Ls, pts) -> (
    -- each entry of Ls is a list/sequence of length 2: {F, G}
    -- where F, G are polynomials in a ring R, (A, phi) are obtained from
    -- genericLinearMap.  We return the ideal of constraints in T
    -- for which phi(F) = G, for all pairs {F,G} in Ls.
    -- We also return A0, phi0 corresponding to these constraints.
    T := ring A;
    A0 := A;
    TR := target phi;
    I := sum for L in Ls list ideal sub(last coefficients(phi L_0 - L_1), T);
    I)

TEST ///
  debug StringTorics
  R = QQ[a..d]
  (A, phi) = genericLinearMap R
  TR = target phi
  assert(source phi === TR)
  assert(ring A === coefficientRing TR)
  for i from 0 to 3 do 
    assert(phi TR_i == (A^{i} * (transpose vars TR))_(0,0))
    
  linearEquationConstraints(A, phi, {}, {
          {{1,0,0,0}, {1,1,0,0}},
          {{0,1,0,0}, {1,1,3,7}},
          {{0,0,1,0}, {5,6,-2,8}},
          {{0,0,0,1}, {0,1,0,0}}}
      )

  F1 = -2*a^3-6*a^2*b+6*b^2*c-12*a^2*d+12*a*b*d+12*b^2*d+36*b*c*d+30*a*d^2+60*b*d^2+54*c*d^2+76*d^3
  F7 = -2*a^3+6*a^2*b-6*a*b^2+2*b^3-6*a*c^2-4*c^3+6*a^2*d-6*a*d^2+2*d^3

  (A0, phi0, I1) = linearEquationConstraints(A, phi, {
          {b+2*d, a-d},
          {b+3*d, a+c},
          {a+b+2*d, a-b},
          {F1, F7}
          }, {
          }
      )
  phi0 F1 == F7

  -- this one isn't correct yet.
  (A0, phi0, I1) = linearEquationConstraints(A, phi, {
          {b+2*d, a-d},
          {b+3*d, a+c},
          {a+b+2*d, a-b}
          }, {
          {{0,0,1,0}, {1,1,-1,1}}
          }
      )
--          

  A0 % sub(trim ideal last coefficients(phi0 F1 - F7), coefficientRing TR)
///

findMaps = (top1, top2, A, phi, RQ) -> (
    TR := target phi;
    n := numgens TR;
    T := coefficientRing TR;
    toTR := f -> sub(f, TR);
    toQQ := f -> sub(f, RQ);
    evalphi := (F,G) -> trim sub(ideal last coefficients((phi toTR F) - toTR G), T);
    (L1, F1, h11, h12) := toSequence top1;
    (L2, F2, l11, l12) := toSequence top2;    
    RZ := ring L1;
    if RZ =!= ring F1 or RZ =!= ring L2 or RZ =!= ring F2 then error "expected polynomials over the same ring";
    if h11 != l11 or h12 != l12 then return null;
    I := (evalphi(L1, L2) + evalphi(F1, F2));
    if I == 1 then return null;
    -- first see if there is a unique solution.
    -- if codim I === n*n and degree I === 1 then (
    --     A0 := A % I;
    --     if support A0 === {} then (
    --         A0 = lift(A0, QQ);
    --         phi0 := map(RQ, RQ, transpose A0);
    --         if phi0 toQQ L1 != toQQ L2 or phi0 toQQ F1 != toQQ F2 then error "map is not correct!";
    --         return (A0, phi0)
    --         );
    --     );
    -- now let's look through all of the components for a smooth point.
    compsI := decompose I;
    As := for c in compsI list A % c;
    As = for a in As list try lift(a, ZZ) else continue; -- grab the ones that lift.
    As = select(As, a -> (d := det a; d == 1 or d == -1));
    if #As > 0 then (
        A0 := As#0;
        phi0 := map(RZ, RZ, transpose A0);
        if phi0 L1 != L2 or phi0 F1 != F2 then error "map is not correct!";
        (A0, phi0)
        )
    else (
        if any(compsI, c -> codim c < n*n or degree c =!= 1) then (
            << "warning: there might be a map in this case!" << endl;
            << netList compsI << endl;
            << "----------------------------------" << endl;
            compsI 
            )
        else null
        )
    )


findIsomorphism = method()
findIsomorphism(CalabiYauInToric, CalabiYauInToric) := (X1, X2) -> (
    (L1, F1, h11, h12) := (c2Form X1, cubicForm X1, hh^(1,1) X1, hh^(1,2) X1);
    (L2, F2, l11, l12) := (c2Form X2, cubicForm X2, hh^(1,1) X2, hh^(1,2) X2);
    if h11 != l11 or h12 != l12 then return null;
    findIsomorphism((L1,F1), (L2,F2))
    )

findIsomorphism(Sequence, Sequence) := (LF1, LF2) -> (
    (L1, F1) := LF1;
    (L2, F2) := LF2;
    RZ := ring L1;
    if RZ =!= ring F1 or RZ =!= ring L2 or RZ =!= ring F2 then error "expected polynomials over the same ring";
    RQ := QQ (monoid RZ);
    (A, phi) := genericLinearMap RQ;
    TR := target phi;
    n := numgens TR;
    T := coefficientRing TR;
    toTR := f -> sub(f, TR);
    toQQ := f -> sub(f, RQ);
    -- Two ways we can use mappinginfo:
    -- the first is simple: we take phi(F) = G.
    -- the second is when we want one ideal to map to another:
    -- each generator of the first ideal must map to an element of the second ideal.
    evalPhi := (F,G) -> sub(ideal last coefficients((phi toTR F) - toTR G), T);
    evalPhiIdeal := (I1,I2) -> sub(ideal last coefficients((phi toTR I1) % toTR I2), T);
    I := (evalPhi(L1, L2) + evalPhi(F1, F2));
    if I == 1 then return null; -- this isn't so good: if we have trouble finding this,
      -- we must deal with that.
    -- first see if there is a unique solution.
    -- if codim I === n*n and degree I === 1 then (
    --     A0 := A % I;
    --     if support A0 === {} then (
    --         A0 = lift(A0, QQ);
    --         phi0 := map(RQ, RQ, transpose A0);
    --         if phi0 toQQ L1 != toQQ L2 or phi0 toQQ F1 != toQQ F2 then error "map is not correct!";
    --         return (A0, phi0)
    --         );
    --     );
    -- now let's look through all of the components for a smooth point.
    compsI := decompose I;
    As := for c in compsI list A % c;
    As = for a in As list try lift(a, ZZ) else continue; -- grab the ones that lift.
    As = select(As, a -> (d := det a; d == 1 or d == -1));
    if #As > 0 then (
        A0 := As#0;
        phi0 := map(RZ, RZ, transpose A0);
        if phi0 L1 != L2 or phi0 F1 != F2 then error "internal logic error: map is not correct!";
        (A0, phi0)
        )
    else (
        if any(compsI, c -> codim c < n*n or degree c =!= 1) then (
            << "warning: there might be a map in this case!" << endl;
            << netList compsI << endl;
            << "----------------------------------" << endl;
            compsI 
            )
        else null
        )
    )

determineIsomorphism = method()
determineIsomorphism(Sequence, Matrix, RingMap, Ring) := (Ts, A, phi, RQ) -> (
    -- Ts is a Sequence (Ti, Tj), where each Ti, Tj is a List:
    --  {h11, h12, c2 form, cubic form}
    -- The latter two are in RZ = ZZ[n vars], and
    -- RQ = QQ[same n vars]
    -- A is n x n generic matrix over a coeff ring `coefficientRing T`, over QQ.
    -- T = (this coeff ring)[same n vars].
    -- phi : T --> T, given by x |--> Ax, or perhaps (transpose A)*x.  TODO: GET THIS RIGHT!
    -- Returns either an n x n invertible matrix A0 over ZZ, or null.
    -- s.t. if phi0 : x |-> A0*x, then phi0 maps the c2 form L1 to L2, cubic form F1 to F2.
    -- If no such map exists, null is returned.  If there may be a map, but we can't 
    -- conclusively find one, then an ideal in the variables of A is returned.
    
    )

partitionH113sByTopology = method()
partitionH113sByTopology(List, HashTable, Ring) := HashTable => (Ls, Ts, RQ) -> (
    -- Ls is a list of labels to separate.
    -- Ts is a hash table of label => {c2, cubicform, h11, h12}
    -- This function first separates these by the invariants: invariantsAll.
    -- The for each pair in each set, it attempts to find a map between them.
    -- output: a hashtable, keys are labels, values are lists of {label, matrix}
    (A, phi) := genericLinearMap RQ;
    if #Ls === 1 then return hashTable {Ls#0 => {}};
    H := partition(lab -> invariantsAll toSequence Ts#lab, Ls);
    distinctTops := new MutableHashTable; -- label => list of {label, matrix}, those with the same topology
    for i in Ls do (
        Ti := Ts#i;
        -- now we attempt to match this with each key of distinctTops
        << "trying " << i << endl;
        prev := keys distinctTops;
        isFound := false;
        for j in prev do (
            Tj := Ts#j;
            ans := findMaps(Tj, Ti, A, phi, RQ);
            --if j == (115,0) and i == (120,0) then error "debug me";
            if ans === null then (
                -- Ti is distinct from Tj
                )
            else if class first ans === Matrix then (
                -- we have a match!
                (A0, phi0) := ans;
                if all(flatten entries A0, a -> liftable(a, ZZ))
                then (
                    isFound = true;
                    distinctTops#j = append(distinctTops#j, {i, lift(A0, ZZ)});
                    break;
                    )
                )
            else (
                << (i,j) << " might be the same, might not" << endl;
                )
            );
        if not isFound then (
            distinctTops#i = {};
            << "found new top: " << i << endl;
            );
        );
    new HashTable from distinctTops
    )


-------------------------------------------------------------------------
-- TODO: remove the following code (any reason to keep it?)
topologyOfCY3 = method(Options => {
        Variable => "x",
        Ring => null
        })

topologyOfCY3(NormalToricVariety, List) := opts -> (V, basisIndices) -> (
    -- input: 
    --   V: a simplicial resolution of a Fano toric 4-fold
    --      X is a (general) anti-canonical section of V.
    --   basisIndices: list of integer indicesas to which V_i will be in the 
    --      basis of Pic X that you choose.
    -- output: a hash table containing:
    --  a. triple intersection numbers (a hash table, H#{a,b,c}, with 0 <= a <= b <= c < h11(X))
    --  b. the h11 numbers: c2(X) . D_i, 0 <= i < h11
    --  c. the integers h11, h12
    --  d. the cubic form C(x,y,z) in a polynomial ring ZZ[x_0, ..., x_(h11-1)]
    --  e. a linear form L(x,y,z) in the same ring, representing c2(X).D_i
    --
    P := convexHull transpose matrix rays V;
    h11 := h21OfCY P; -- we want h11 of `polar P`.
    h21 := h11OfCY P;
    if #basisIndices != h11 then error("expected "|h11|" indices");
    H := CY3NonzeroMultiplicities V;
    -- basisInv := new MutableHashTable;
    -- for i from 0 to #basisIndices-1 do basisInv#(basisIndices#i) = i;
    -- H3 := hashTable for x in keys H list (
    --     if isSubset(x, basisIndices) then (
    --         x' := apply(x, i -> basisInv#i);
    --         x' => H#x 
    --         ) else continue
    --     );
    -- Now let's get the cubic form and the linear form directly from the intersection theory.
    -- For larger h11, this method will need to change.
    x := getSymbol opts.Variable;
    pt := base(x_0..x_(h11-1));
    A := intersectionRing pt; -- over QQ
    R := if opts#Ring =!= null then opts#Ring else ZZ (monoid A);
    if numgens R =!= h11 then error("expected a ring with "|toString h11|" variables");

    X := completeIntersection(V, {-toricDivisor V});
    Xa := abstractVariety(X, pt);
    IX := intersectionRing Xa;
    h := sum(h11, i -> A_i * IX_(basisIndices#i));
    C := sub(integral(h^3), vars R);
    L := integral((chern_2 tangentBundle Xa) * h);
    L = sub(L, vars R);
    (h11, h21, C, L)
    )

hh(Sequence, TopologicalDataOfCY3) := (pq, T) -> (
    (p,q) := pq;
    if p > q then (p, q) = (q, p);
    if p == 0 then (
        if q == 3 or q == 0 then 1 else 0
        )
    else if p == 1 then (
        if q == 1 then T#2
        else if q == 2 then T#3
        else 0
        )
    else if p == 2 then (
        if q == 1 then T#3
        else if q == 2 then T#2 
        else 0
        )
    else if p == 3 then (
        if q == 3 then 1
        else 0
        )
    )

c2Form TopologicalDataOfCY3 := T -> T#0
cubicForm TopologicalDataOfCY3 := T -> T#1
----- end of removing code TODO -----------------------------
