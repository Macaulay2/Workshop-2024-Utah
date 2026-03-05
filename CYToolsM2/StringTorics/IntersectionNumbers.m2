-------------------------------------------------
-- Intersection numbers for Calabi-Yau 3-folds --
-------------------------------------------------
-- currently is functional only for hypersurfaces
-- in a toric 4-fold.  Also, the polytope must be favorable.
-- TODO: remove favorable hypothesis
-- TODO: allow 4D CY's?
-- TODO: allow CI CY's in a Fano toric.

intersectionNumbers = method()
-- intersectionNumbers = method(Options => {Indices => null}) -- Indices: which elements to keep.
--   -- these will be reordered 0, 1, ...,, #opts.Indices-1.
--   -- default for CalabiYauInToric is `basisIndices X`

------------------------
-- New code, Nov/Dec 2022.
-- This will replace the older code below, which uses the intersection ring.
-- However, we want to keep that code to check the results.
-- Also, the interface is currently totally different.
--
-- TODO:
--   this function should stash its value in X?
--   we need to stash:
--    (a) intersectionNumbers (using basis, translating to 0, 1, ...): HashTable, lkeys are triples of ints
--        and values are integer intersection numbers.
--    (b) intersectionNumbers (full, using all toric divisors)
--    (c) c2 (a list of c2.D_(i1), ... where i1, ... is the basis of toric divisors.
--    (d) c2 (full)
-- Functions we need here:
--   the ring RZ should be given when X is created?  If not, M2 will create one?
--   cubicForm X
--   c2 X
--   to/from intersection number tables and cubic form (for testing easily with older code).
--   computeIntersectionNumbers(A, basis, T2)
--   computeIntersectionNumbers(X) -- stashes results into X,
--   intersectionNumbers X -- stashes value if not yet computed.
--   c2 X -- stashes value if not yet computed.
--   intersectionNumbers(X, Basis => Full)
--   c2(X, Basis => Full)
--   dump(X): should dump these values
--   cyData(String, ...): should read these values
--   topology X -- returns TopologicalDataOfCY3: h11, h12, c2, cubicForm.

-- dump: make sure we write out info about c2, intersectionNumbers as well.
--   and can read these back in.
--   also write out two face triangulation...
--
-- topologicalData: return a type with (h11, h12, intersectionnumbers, c2 list)
--   no ring involved.
------------------------

-- Old toBasisIntersectionNumbers moved to ScratchIntersectionNumbers.m2

toBasisIntersectionNumbers = (toricIntersectionNumbers, basIndices, nonfavsHash) -> (
    -- toricIntersectionNumbers: list of {i,j,k} => intersection number.
    -- basIndices: list of divisors and nonfavorable divisors making up the basis for Pic X.
    --   e.g. {0,1,2,(5,0),(5,1),(5,2)}
    -- nonfavsHash: hashtable, keys are nonfavorable toric divisors (so integers), and the value
    --   is the list {genus of 2-face, index of 2-face}.  This last is used to determine if two
    --   nonfavorable toric divisors lie on the same 2-face (if not, they definitely have intersection 0).
    -- Result: list of {i,j,k} => num, now they are all integers in the range 0..#basIndices-1.
    --   and num is a non-zero integer.  All other triple intersections are zero.
    allBasisDivisors := set sort unique (basIndices/(x -> if instance(x, ZZ) then x else x#0));
    basisIndex := hashTable for i from 0 to #basIndices-1 list basIndices#i => i;
    flatten for t in toricIntersectionNumbers list (
        if not isSubset(t#0, allBasisDivisors) then continue else (
            nonfavs := select(t#0, i -> nonfavsHash#?i);
            if #nonfavs === 0 then {t#0/(a -> basisIndex#a)//sort => t#1}
            else (
                if nonfavs/(i -> nonfavsHash#i#1)//unique//length > 1 then continue;
                g := nonfavsHash#(nonfavs#0)#0; -- they are all the same 2-face, so should all have the same genus.
                for ell from 0 to g list (
                    tnew := sort for t1 in t#0 list if nonfavsHash#?t1 then basisIndex#(t1,ell) else basisIndex#t1;
                    if t#1 % (1+g) != 0 then error "internal error: logic is messed up here!";
                    tnew => (t#1 // (1+g))
                    )
                )
            )
        )
    )

-- computeToricIntersectionNumbers: An internal function for computeIntersectionNumbers
-- this is the workhorse function for computing intersection numbers on X.
computeToricIntersectionNumbers = method()
computeToricIntersectionNumbers(Matrix, List) := (A, T2) -> (
    intnums1 := hashTable flatten for t in T2 list for s in t#2 list s => t#3+1;
    E := sort unique flatten for t in T2 list flatten for s in t#2 list subsets(s, 2);
    -- get intersection numbers {i,i,j} or {i,j,j}
    intnums2 := hashTable flatten for e in E list (
        (i,j) := toSequence e;
        -- find intersection numbers for {i,i,j}, {i,j,j}.
        A1 := A_e;
        b := sum for k from 0 to numcols A - 1 list (
            if k == i or k == j then continue else (
                t := sort{i,j,k};
                (if intnums1#?t then intnums1#t  else 0) * A_{k}
                )
            );
        x := flatten entries solve(A1,-b); -- this is done over ZZ.  I think that is fine here...
        -- TODO: check that x is correct: A1*x == -b?
        select({{i,i,j} => x#0, {i,j,j} => x#1}, x -> x#1 != 0)
        );
    -- Now we get the triple intersection numbers.
    intnums3 := hashTable for i from 0 to numcols A - 1 list (
        b := sum for k from 0 to numcols A - 1 list (
            if k == i then continue else (
                t := sort{i,i,k};
                (if intnums2#?t then (intnums2#t)  else 0) * A_{k}
                )
            );
        x := flatten entries solve(A_{i},-b); -- this is done over ZZ.  I think that is fine here...
        -- TODO: check that x is correct: A_{i}*x == -b?
        if x#0 == 0 then continue else  {i,i,i} => x#0
        );
    for x in sort join(pairs intnums1, pairs intnums2, pairs intnums3) list x#0 => x#1
    --FIXME: this should return a list of {i,j,k} => a, not ({i,j,k},a)
    -- (to match intersectionNumbers)
    )

-- computeC2: An internal function for computeIntersectionNumbers
computeC2 = method()
-- computeC2(List, List) := (toricIntersectionNumbers, basIndices) -> (
--     topval := toricIntersectionNumbers/first/max//max;
--     H := hashTable toricIntersectionNumbers;
--     for a in basIndices list (
--         -- add up all intersection numbers {a,i,j}, i<j
--         sum flatten for i from 0 to topval list for j from i+1 to topval list (
--             t := sort {a,i,j};
--             if H#?t then H#t else continue
--             )
--         )
--     )

computeC2(List, List, HashTable) := (toricIntersectionNumbers, basIndices, nonfavsHash) -> (
    topval := toricIntersectionNumbers/first/max//max;
    H := hashTable toricIntersectionNumbers;
    for a in basIndices list (
        -- add up all intersection numbers {a,i,j}, i<j
        -- BUT: if a is (alpha,ell) is nonfavorableon 2-face with genus g, want {alpha,i,j}//(1+g)
        if instance(a, ZZ) then (
            sum flatten for i from 0 to topval list for j from i+1 to topval list (
                t := sort {a,i,j};
                if H#?t then H#t else continue
                )
            )
        else (
            -- here, a is non-favorable.
            alpha := a#0;
            g := nonfavsHash#alpha#0;
            sum flatten for i from 0 to topval list for j from i+1 to topval list (
                t := sort {alpha,i,j};
                if H#?t then (
                    if H#t % (1+g) != 0 then error "internal error: logic is wrong for nonfavorables";
                    H#t // (1+g)
                    )
                else continue
                )
            )
        )
    )

-- computeIntersectionNumbers: An internal function for intersectionNumbers. toricIntersectionNumbers, and c2.
computeIntersectionNumbers = method()
computeIntersectionNumbers CalabiYauInToric := X -> (
    if not X.cache#?"toricIntersectionNumbers" then  (
        Q := cyPolytope X;
        basIndices := basisIndices Q;
        A := transpose matrix rays Q;
        nonfavs := hashTable findTwoFaceInteriorDivisors Q;
        T2 := restrictTriangulation X;
        result := computeToricIntersectionNumbers(A, T2);
        X.cache#"toricIntersectionNumbers" = result;
        X.cache#"intersectionNumbers" = toBasisIntersectionNumbers(result, basIndices, nonfavs);
        X.cache#"c2" = computeC2(result, basIndices, nonfavs);
        );
    )

--------------------------------------------
-- intersection number interface routines --
--------------------------------------------
intersectionNumbers CalabiYauInToric := X -> (
    computeIntersectionNumbers X;
    X.cache#"intersectionNumbers"
    )

toricIntersectionNumbers = method()
toricIntersectionNumbers CalabiYauInToric := X -> (
    computeIntersectionNumbers X;
    X.cache#"toricIntersectionNumbers"
    )

c2 = method();
c2 CalabiYauInToric := X -> (
    computeIntersectionNumbers X;
    X.cache#"c2"
    )

-- maybe: toricIntersectionNumbers, c2, c2Form, intersectionForm.
--
-- intersectionNumbers X, intersectionNumbers(X, Full => true), intersectionForm X
-- c2 X, c2Form X.

-- TODO: if X is not favorable, need to redo the basis, and intersection numbers (and also then the c2 form)


-----------------------------------------------
-- Utility functions --------------------------
-- Used to translate between data formats -----
-----------------------------------------------
exponentToProduct = exp -> (
    -- exp is a list of integers, e.g. {0,3,1}
    -- result is expanded to a product, e.g. {1,1,1,2}
    -- the result includes integers in 0..#exp - 1, in ascending order.
    -- e.g.
    --   exponentToProduct {0,3,1} == {1,1,1,2}
    flatten for i from 0 to #exp-1 list toList(exp#i : i)
    )

productToExponents = (prod, nvars) -> (
    -- prod is a list of ascending integers, in range 0..nvars-1
    -- e.g. {0,1,1,2,4}
    -- this is translated to an exponent vector,
    -- e.g. 
    --   productToExponents({0,1,1,2,4}, 6) == {1,2,1,0,1,0}
    -- this could be faster if needed
    T := tally prod;
    for i from 0 to nvars - 1 list if T#?i then T#i else 0
    )

multinomial = exp -> (
    -- exp is an exponent vector
    -- returns an integer
    -- e.g. 
    --   multinomial {1,0,2} == 3
    --   multinomial {1, 1, 1} == 6
    n := sum exp;
    den := product for i from 0 to #exp - 1 list (exp#i)!;
    n! // den
    )

toCOO = method()
toCOO RingElement := (F) -> (
    sort for f1 in listForm F list (
        e := first f1; -- exponents
        c := last f1; -- coeff
        d := multinomial e;
        e' := exponentToProduct e;
        e' => if c % d ==0 then c//d else c/d -- TODO: not the exponent vector!!
        )
    )

toRingElement = method()
toRingElement(List, Ring) := (f, RZ) -> (
    if #f == 0 then return 0_RZ;
    n := sum f#0#0;
    sum for f1 in f list (
        e' := f1#0; -- list of variable products, e.g. {0,1,1,2}
        c := f1#1; -- coeff (an integer, currently envisioned)
        e := productToExponents(e', numgens RZ);
        d := multinomial e;
        d * c * RZ_e
        )
    )

c2Form = method()
c2Form CalabiYauInToric := RingElement => X -> (
    RZ := picardRing X;
    ((vars RZ) * transpose matrix {c2 X})_(0,0)
    )

cubicForm = method()
cubicForm CalabiYauInToric := RingElement => X -> (
    RZ := picardRing X;
    toRingElement(intersectionNumbers X, RZ)
    )


------------------------------------
-- Intersection numbers via intersection ring in Schubert2
-- This is an alternate method, slower than intersectionNumbers, 
-- but that can be used to test intersectionNumbers.

-- Alternate to intersectionNumbers, slower.  But easier code.  Only
-- works for the case when the ambient toric variety has no torsion in
-- the Class group.

-- Simple subroutine for finding the list of indices for possible intersections.
-- e.g. if in the resulting list, {0,1,1} appears, then this will represent the
-- product H_0 . H_1 . H_1 (which is an integer)
monoms = (deg, lo, hi) -> (
    -- input: deg, lo, hi: all integers
    -- output: a list of lists of integers all of length 'deg',
    --   sorted in ascending order.
    if deg == 0 then {{}}
    else if lo === hi then {splice{deg:lo}}
    else
    flatten for i from lo to hi list (
        L1 := monoms(deg-1, i, hi);
        for t in L1 list prepend(i, t)
        )
    )

intersectionNumbersOfCY = method()

intersectionNumbersOfCY(Ring, List) := HashTable => (IX, basisIndices) -> (
    -- IX: should be a ring produced for Schubert2, having 'integral' function for top degree elements.
    -- basisIndices is a subList of {0, ..., numgens IX - 1}.
    -- WARNING: this is cubic in number of generators of IX.  This can be improved,
    -- using the toric structure of X as a hypersurface in a toric V.
    -- TODO WARNING: the 3 in here is for 3-folds...!
    mons := monoms(3, 0, #basisIndices-1);
    bas := for i in basisIndices list IX_i;
    for t in mons list (
        m := integral product(t, i -> bas_i);
        a := lift(integral product(t, i -> bas_i), ZZ);
        if a === 0 then continue else t => a
        )
    )

intersectionNumbersOfCY(NormalToricVariety, List) := (V, basisIndices) -> (
    X := completeIntersection(V, {-toricDivisor V});
    Xa := abstractVariety(X, base());
    IX := intersectionRing Xa;
    intersectionNumbersOfCY(IX, basisIndices)
    )
intersectionNumbersOfCY CalabiYauInToric := X -> (
    intersectionNumbersOfCY(ambient X, basisIndices X)
    )

-- end of intersectionNumbersOfCY, alternate slower method for
-- computing intersection numbers of a CY 3-fold.
------------------------------------------
    


------------------------------
-- tripleProductsCY --
------------------------------
-- Simpler code, used to debug the algorithm/implementation above.
tripleProductsCY = method()
tripleProductsCY NormalToricVariety := (V) -> (
    elapsedTime AV := abstractVariety(V, point);
    IV := intersectionRing AV;
    h := sum gens IV; -- Calabi-Yau hyperplane class in V.
    J := ideal select((ideal IV)_*, f -> size f == 1);
    forceGB gens J;
    gens gb J;
    A := (ring J)/J;
    monoms := ideal basis(3, A);
    elapsedTime JV := sub(monoms,IV);
    elapsedTime (JVh := h ** (gens JV));
    flatJVh := flatten entries JVh;
    hashTable for i from 0 to numgens monoms - 1 list (
        m := monoms_i;
        d := integral(flatJVh#i);
        if d > 0 then m => d else continue
        )
    )

------------------------------
-- possibleNonZeros --
------------------------------
possibleNonZeros = (V) -> (
    -- assumption currently: V has dim 4, is reflexive, and X is the anti-canonical CY3 divisor.
    -- returns a list of lists of 3 integers (0 <= i1 <= i2 <= i3 <= N-1)
    --  where N = #rays V.
    -- and all triples other than those on this list must have triple intersection
    -- on X being zero.
    P2 := convexHull transpose matrix rays V;
    F := annotatedFaces P2;
    faces2 := select(F, f -> f#0 == 2);
    faces2 = faces2/(x -> x#2); -- this is a list of all 2-faces in the polytope,
    -- with which rays are on each face.
    -- any triple not supported on a 2-face will have triple intersection zero.
    triangles := (max V)/(t -> subsets(t,3))//flatten//unique;
    edges := (max V)/(t -> subsets(t,2))//flatten//unique//sort;
    triples := sort flatten for f in faces2 list select(triangles, t -> isSubset(t,f));
    singles := for i from 0 to # rays V - 1 list {i,i,i};
    doubles := sort flatten for f in faces2 list select(edges, t -> isSubset(t,f));
    doubles = unique flatten for x in doubles list {{x#0,x#0,x#1},{x#0,x#1,x#1}};
    {singles,doubles,triples}
    )

--------------------------------------
-- CY3NonzeroMultiplicities --
--------------------------------------
  CY3NonzeroMultiplicities = method()
  CY3NonzeroMultiplicities NormalToricVariety := (V) -> (
      RAYS := transpose matrix rays V;
      P2 := convexHull RAYS;
      (singles,doubles,triples) := toSequence possibleNonZeros V;
      doubles = doubles/unique/sort//unique;
      singles = singles/unique/sort//unique;
      mult3 := new MutableHashTable from for x in triples list (
           x => 1 + genus(P2, minimalFace(P2, (rays V)_x))
           );
      multvec := ij -> (
          for ell from 0 to #rays V-1 list (
              if member(ell,ij) then 0
              else (
                  s := sort append(ij,ell);
                  if mult3#?s then mult3#s else 0
                  ))
          );
      for d in doubles do (
          RHS := - RAYS *  transpose (matrix{multvec d});
          vals := flatten entries solve(RAYS_d, RHS);
          d1 := prepend(d#0,d);
          d2 := append(d,d#1);
          if vals#0 != 0 then mult3#d1 = vals#0;
          if vals#1 != 0 then mult3#d2 = vals#1;
          );
      for d in singles do (
          s := {d#0,d#0};
          RHS := - RAYS *  transpose (matrix{multvec s});
          vals := flatten entries solve(RAYS_d, RHS);
          if vals#0 != 0 then mult3#{d#0,d#0,d#0} = vals#0;
          );
      new HashTable from mult3
      )

------------------------------
-- CY3Intersections --
------------------------------
CY3Intersections = method()
CY3Intersections(NormalToricVariety, List) := (V, indexOfDs) -> (
    H := CY3NonzeroMultiplicities V;
    loc := new HashTable from for i from 0 to #indexOfDs-1 list indexOfDs#i => i;
    tr := k -> sort for k1 in k list loc#k1;
    new Array from for kv in pairs H list (
      if not isSubset(kv#0, indexOfDs) then continue;
      new Array from append(tr kv#0, kv#1)
      )
  )

