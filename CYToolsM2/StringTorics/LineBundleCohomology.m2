-- several things here:
--  1. line bundle cohomology on V (NormalToricVarieties, cohomcalg)
--  2. basis as lists of (annotated) Laurent monomials)
--  3. matrix of multiplication by F.
--  4. compute ranks of matrices.  What should be stashed?  bases -- probably, matrices -- not sure, ranks -- not sure

----------------------------------------------------------------
-- Cohomology of line bundles: obtaining a basis of fractions --
----------------------------------------------------------------
toricCohomologySetup = method()
toricCohomologySetup NormalToricVariety := (X) -> (
    if not X.cache#?CohomologySetup then (
        -- part 1: compute free res of B^* (over ZZ^n)
        R := ring X;
        R1 := (coefficientRing R)(monoid[gens R, DegreeRank=>numgens R]);
        B1 := monomialIdeal sub(ideal X, vars R1);
        B1' := dual B1;
        C := resolution comodule B1';
        -- part 2: take the (fine) degrees in there, and place them depending on HH^i, in a hash table
        -- part 3: create multi-graded rings for each of these
        -- finally, stash this info into the cache table of X.
        -- outside of this: make a function that takes a degree (or an OO(D)?), and i, and computes the fractions for HH^i(V, OO(V)).
        rawDegs := sort flatten for j from 0 to length C list (degrees C_j)/(x -> (sum x - j,x));
        rawDegsWithRings := for x in rawDegs list (
            -- x is (i, I), where I has length n is an array of 0 and 1's.
            Ra := if all(x#1, i -> i == 0) then R else (
                ds := degrees R;
                degs := for i from 0 to #x#1-1 list if x#1#i == 0 then ds_i else -ds_i;
                (coefficientRing R)(monoid[generators R, Degrees=>degs])
                );
            (x#0, x#1, Ra)
            );
        H := partition(x -> x#0, rawDegsWithRings);
        X.cache#CohomologySetup = applyPairs(H, x -> (x#0, (x#1/(y -> drop(y,1)))));
        );
    X.cache#CohomologySetup
    )
cohomologyBasis = method()
cohomologyBasis(ZZ, NormalToricVariety, List) := (i,X,L) -> (
    H := toricCohomologySetup X;
    L = -L;
    KR := frac ring X;
    if not H#?i then return {};
    degs := H#i;
    flatten for x in degs list (
        -- x is (I, R), I is a vector of 0,1's of length n = numgens R, and R is a multigraded ring
        (I,Ra) := x;
        deg := degree(Ra_I);
        ---M := transpose matrix degrees Ra | transpose matrix{-L+deg};
        ---result := (normaliz(M, 5))#"gen";
        ---elems0 := if result === null then {} else entries result;
        ---elems1 := select(elems0, v -> v#-1 == 1);
        ---elems2 := elems1/(v -> drop(v,-1));
        ---exps := elems2;
        --<< "L=" << L << " deg=" << deg << " -L-deg=" << -L-deg << endl;
        elems := flatten entries basis(-L-deg, Ra);
        --print elems;
        exps := elems/exponents/first;
        -- now need to place these into frac R.
        for e in exps list (
            e1 := for i from 0 to #I-1 list if I#i == 0 then e#i else - e#i - 1;
            KR_e1
            )
        )
    )
genericCohomologyMatrix = method()

-- Compute the map HH^i(V, OO_V(D+K)) --> HH^i(V, OO_V(D))
-- as a mutable matrix, over a finite field.
-- given V, i, and D.
-- the map is given by multiplication by a random form F of HH^0(OO_V(K)).
-- Caveat: it is extemely unlikely, but could happen that the generic rank
--   does not occur: the actual generic rank could potentially be larger than 
--   the rank of this matrix.
-- 
genericCohomologyMatrix(ZZ, NormalToricVariety, List) := (i,V,D) -> (
    K := degree toricDivisor V;
    basis1 := cohomologyBasis(1, V, D+K);
    basis2 := cohomologyBasis(1, V, D);
    basis0 := first entries basis(-K, ring V);
    hash1 := hashTable for i from 0 to #basis1-1 list basis1#i => i;
    hash2 := hashTable for i from 0 to #basis2-1 list basis2#i => i;
    M := mutableMatrix(ZZ/32003, #basis2, #basis1);
    for p in basis0 do (
        c := random ring M;
        for k in keys hash1 do (
            m := p * k;
            if hash2#?m then M_(hash2#m, hash1#k) = c;
            )
        );
    M
    )

cohomologyMatrix = method()
cohomologyMatrix(ZZ, NormalToricVariety, List, RingElement) := (i,V,D,F) -> (
    -- Create the map HH^i(V, D - degree F) --> HH^i(V, D), induced by mult by F (in Cox ring).
    kk := coefficientRing ring F;
    deg := degree F; -- F should be an element of the Cox ring.
    basis1 := cohomologyBasis(i, V, D-deg);
    basis2 := cohomologyBasis(i, V, D);
    basis0 := terms F;
    hash1 := hashTable for i from 0 to #basis1-1 list basis1#i => i;
    hash2 := hashTable for i from 0 to #basis2-1 list basis2#i => i;
    M := mutableMatrix(kk, #basis2, #basis1);
    for tm in basis0 do (
        c := leadCoefficient tm;
        p := leadMonomial tm;
        for k in keys hash1 do (
            m := p * k;
            if hash2#?m then M_(hash2#m, hash1#k) = c;
            )
        );
    (M, basis2, basis1)
    )

cohomologyMatrixRank = method()
cohomologyMatrixRank(ZZ, NormalToricVariety, List, RingElement) := (i,V,deg,F) -> (
    z := cohomologyMatrix(i,V,deg,F);
    nrows := #z#1;
    ncols := #z#2;
    rk := rank z#0;
    (nrows, ncols, rk)
    )

-- The above code doesn't work on larger examples.  Let's see if we can do better.
toricOrthants = method()
toricOrthants NormalToricVariety := (X) -> (
        -- part 1: compute free res of B^* (over ZZ^n)
        R := ring X;
        R1 := (coefficientRing R)(monoid[gens R, DegreeRank=>numgens R]);
        B1 := monomialIdeal sub(ideal X, vars R1);
        B1' := dual B1;
        C := resolution comodule B1';
        -- part 2: take the (fine) degrees in there, and place them depending on HH^i, in a hash table
        -- part 3: create multi-graded rings for each of these
        -- finally, stash this info into the cache table of X.
        -- outside of this: make a function that takes a degree (or an OO(D)?), and i, and computes the fractions for HH^i(V, OO(V)).
        rawDegs := sort flatten for j from 0 to length C list (degrees C_j)/(x -> (sum x - j,x));
        rawDegs/(x -> (x#0, positions(x#1, i -> i != 0)))
    )

H1orthants = method()
H1orthants NormalToricVariety := (X) -> (
        R := ring X;
        R1 := (coefficientRing R)(monoid[gens R, DegreeRank=>numgens R]);
        B1 := monomialIdeal sub(ideal X, vars R1);
        B1' := dual B1;
        R2 := (coefficientRing R)(monoid[gens R]);
        B1'' := sub(B1', R2);
        I := ideal select(flatten entries gens B1'', m -> degree m === {2});
        C := resolution(I, DegreeLimit=>1);
        C1 := for i from 1 to length C list sub(C.dd_i, R1);
        C2 := new MutableList;
        C2#0 = map(target C1_0,, C1_0);
        for i from 1 to length C1-1 do 
          C2#i = map(source C2#(i-1),, C1#i);
        degs := flatten for phi in toList C2 list degrees source phi;
        return for d in degs list positions(d, a -> a != 0);
        -- part 2: take the (fine) degrees in there, and place them depending on HH^i, in a hash table
        -- part 3: create multi-graded rings for each of these
        -- finally, stash this info into the cache table of X.
        -- outside of this: make a function that takes a degree (or an OO(D)?), and i, and computes the fractions for HH^i(V, OO(V)).
        rawDegs := sort flatten for j from 0 to length C list (degrees C_j)/(x -> (sum x - j,x));
        rawDegs/(x -> (x#0, positions(x#1, i -> i != 0)))
    )

-- goal: given V, find a permutation P, a d x d invertible (over ZZ) matrix Q, s.t. QAP = [-I | C]
-- where A == transpose matrix degrees ring V (assuming we can compute that ring!  If not, we need to fix that later).
normalDegrees = method()
normalDegrees NormalToricVariety := (V) -> (
    if not V.cache.?normalDegrees then V.cache.normalDegrees = (
        S := ring V;
        A := transpose matrix degrees S;
        d := numRows A;
        n := numColumns A; -- also numgens S.
        good := select(1, max V, a -> (d := det A_a; d == 1 or d == -1));
        if good === {} then error "cannot find maximal cone with volume 1";
        cols := good#0; -- ordered list of d column indices of the (d x n) matrix A.
        others := sort toList (set (0..n - 1) - set cols);
        perm := join(cols, others); -- or is the inverse permutation!?
        permInv := hashTable for i from 0 to #perm-1 list perm#i => i;
        permInv = for i from 0 to #perm-1 list permInv#i;
        Q := A_cols;
        Anew := Q^-1 * A_perm;
        C := Anew_{d..n-1};
        -- newdeg: take a degree in ZZ^d (old basis, using A), and change it to 
        --  a degree in ZZ^d, wrt the new basis.
        newdeg := (deg1) -> flatten entries (Q^-1 * transpose matrix {deg1});
        -- newneg: take a subset of {0..n-1} (i.e. columns of A), and return the sorted
        --  list of the corresponding columns of Anew.
        newneg := (neg1) -> sort for x in neg1 list permInv#x;
        -- mon is a row vector  0..n-1 w.r.t the new set of columns in Anew.
        -- this function translates it back to a row vector w.r.t A,
        -- and then provides a monomial in S (all with respect to the original
        -- set of variables).
        tofrac := (mon) -> (
            v1 := flatten entries (mon^permInv); 
            product for i from 0 to #v1-1 list (
                if v1_i > 0 then 
                    S_i^(v1_i) 
                else if v1_i < 0 then 
                    (1/S_i^(-v1_i)) 
                else 1_S
            ));
        -- A is a map ZZ^n --> ZZ^d
        -- Anew is as well.
        -- These induce an isomorphism ZZ^d (old) --> ZZ^d (new).
        {Anew, perm, Q, newdeg, newneg, tofrac, permInv}
        );
    V.cache.normalDegrees#0
    )

normalDegree = method()
normalDegree(NormalToricVariety, List) := (V, deg1) -> (
    normalDegrees V; -- we don't use the value, just the cached values.
    V.cache.normalDegrees#3 deg1
    )

setOfColumns = method()
setOfColumns(NormalToricVariety, List) := (V, cols) -> (
    normalDegrees V; -- we don't use the value, just the cached values.
    V.cache.normalDegrees#4 cols
    )

getFraction = method()
getFraction(NormalToricVariety, List) := (V, mon) -> (
    normalDegrees V; -- we don't use the value, just the cached values.
    V.cache.normalDegrees#5 mon
    )

-- This one is now correct, I think.
findTope = method()
findTope(Matrix, List, List) := (normalA, normalNeg, normalDeg) -> (
    n := numColumns normalA;
    d := numRows normalA;
    C := normalA_{d..n-1};
    I := id_(ZZ^(n-d));
    alpha := if #normalNeg == 0 then normalDeg else
      normalDeg + flatten entries sum(normalNeg, p -> normalA_p);
    -- create hyperplanes matrix M, RHS b, polytope will be Mx <= b
    hypers := mutableMatrix(C || -I);
    RHS := mutableMatrix transpose matrix{join(alpha, toList(n-d: 0))};
    --return (hypers, RHS);
    for p in normalNeg do (
        rowMult(hypers, p, -1);
        rowMult(RHS, p, -1);
        );
    (polyhedronFromHData(matrix hypers, matrix RHS), hypers, RHS)
    )

cohomologyFractions = method()
cohomologyFractions(NormalToricVariety, List, List) := (V, negativeSet, deg) -> (
    -- returns a list of fractions in Ring S.
    Anew := normalDegrees V;
    d := numRows Anew;
    n := numColumns Anew;
    normalNeg := setOfColumns_V negativeSet;
    beta := normalDegree(V, deg);
    normalNegDegree := (sum for i in normalNeg list flatten entries Anew_i);
    gamma := beta + normalNegDegree;
    P := first findTope(Anew, normalNeg, beta);
    if not isCompact P then error "Your negative set doesn't correspond to a cohomology cone";
    if isEmpty P then return {};
    LP := latticePoints P;
    C := Anew_{d..n-1};
    for lp in LP list (
        mon := (transpose matrix{gamma} - C * lp) || lp;
        (V.cache.normalDegrees#5 mon)/(product(negativeSet, i -> (ring V)_i))
        )
    )

--------------------------------------------------------------
-- Cohomology for complete intersections in toric varieties --
--------------------------------------------------------------

cohomologyVector = method()
cohomologyVector(NormalToricVariety, ToricDivisor) := (Y,D) -> cohomologyVector(Y, degree D)
cohomologyVector(NormalToricVariety, List) := (Y, D) -> (cohomCalg(Y, {D}); first Y.cache.CohomCalg#D)

cohomology(ZZ, NormalToricVariety, List, RingElement) := opts -> (i,X,deg,F) -> (
    (nrows, ncols, rk1) := cohomologyMatrixRank(i, X, deg, F);
    (nrows2, ncols2, rk2) := cohomologyMatrixRank(i+1, X, deg, F);
    nrows-rk1 + ncols2 - rk2
    )

cohomology(ZZ, CompleteIntersectionInToric, List, RingElement) := opts -> (i,X,deg,F) -> (
    if dim X != dim ambient X - 1 then error "cohomology for line bundles of CI's of codimension >= 2 is not yet handled";
    V := ambient X;
    (nrows, ncols, rk1) := cohomologyMatrixRank(i, V, deg, F);
    (nrows2, ncols2, rk2) := cohomologyMatrixRank(i+1, V, deg, F);
    nrows-rk1 + ncols2 - rk2
    )

cohomologyVector(CompleteIntersectionInToric, List, RingElement) := (X, deg, F) -> (
    -- TODO: allow full complete intersection here... Not just a hypersurface.
    -- require currently that X has codim 1 in a toric variety
    if dim X != dim ambient X - 1 then error "cohomology for line bundles of CI's of codimension >= 2 is not yet handled";
    V := ambient X;
    rks := for i from 0 to dim ambient X list cohomologyMatrixRank(i, V, deg, F);
    -- rks is a list of (nrows, ncols, rk).
    for j from 0 to dim X list (
        (nrows1, ncols1, rk1) := rks#j;
        (nrows2, ncols2, rk2) := rks#(j+1);
        nrows1 - rk1 + ncols2 - rk2
        )
    )

cohomology(ZZ, CalabiYauInToric, List, RingElement) := opts -> (i,X,deg,F) -> (
    V := ambient X;
    (nrows, ncols, rk1) := cohomologyMatrixRank(i, V, deg, F);
    (nrows2, ncols2, rk2) := cohomologyMatrixRank(i+1, V, deg, F);
    nrows-rk1 + ncols2 - rk2
    )

cohomologyVector(CalabiYauInToric, List, RingElement) := (X, deg, F) -> (
    V := ambient X;
    rks := for i from 0 to dim ambient X list cohomologyMatrixRank(i, V, deg, F);
    -- rks is a list of (nrows, ncols, rk).
    for j from 0 to 3 list (
        (nrows1, ncols1, rk1) := rks#j;
        (nrows2, ncols2, rk2) := rks#(j+1);
        nrows1 - rk1 + ncols2 - rk2
        )
    )

cohomologyVector LineBundle := List => L -> (
    X := variety L;
    Fs := equations X;
    if #Fs =!= 1 then error "alas, complete intersection line bundle cohomology not yet implemented";
    cohomologyVector(X, degree L, Fs#0)
    )

hh(ZZ, LineBundle) := ZZ => (i,L) -> (
    X := variety L;
    Fs := equations X;
    cohomology(i, X, degree L, Fs#0)
    )

ScriptedFunctor ^* := (scriptedfun) -> if scriptedfun === hh then cohomologyVector else null

-------------------------------------------
-- Cohomology from short exact sequences --
-------------------------------------------
cohomologyFromLES = method()
cohomologyFromLES Matrix := (M) -> (
    les := flatten entries M;
    eqns := trim splitLES les;
    M % eqns
    )

splitLES = L -> (
    alternatingSum := (P) -> (sign := 1; sum for i from 0 to #P-1 list (result := sign * P#i; sign = -sign; result));
    pos := positions(L, x -> x == 0);
    if #pos == 0 then return ideal alternatingSum L;
    pos = {-1}|pos|{#L};
    ideal for i from 1 to #pos-1 list (
        alternatingSum L_{pos#(i-1)+1..pos#i-1}
        )
    )
removeUnusedVariables = method()
removeUnusedVariables List := (L) -> (
    -- L is a list of matrices
    R := ring L#0;
    keepthese := L/support/set//sum//toList//sort;
    A := (coefficientRing R)[keepthese];
    phi := map(A,R);
    for m in L list phi m
    )

collectLineBundles = method()
collectLineBundles(NormalToricVariety, List, List) := (Y, CI, Ds) -> (
    -- Y is a projective normal toric variety
    -- CI is a list of multidegrees, i.e. elements in Cl(Y).
    -- Ds is another such list
    -- Result: a list of (r, {...}) of cohomologies that need to be computed, where r is the number of
    --  CI divisors being used (0 <= r <= #CI).
    N := #CI;
    degRank := # rays Y - dim Y;
    koszuls := for v in subsets CI list (
        if #v > 0 then degree(-sum v) else toList(degRank:0)
        );
    rsort unique flatten for D in Ds list (
        d := if instance(D, ToricDivisor) then degree D else D;
        for h in koszuls list (h+d)
        )
    )

--------------------------------------------
-- Complete Intersection in Toric Variety --
--------------------------------------------
initializeCohomologies = method(Options => {Symbol => getSymbol "a"})
initializeCohomologies(CompleteIntersectionInToric, ZZ, ZZ) := opts -> (X, nlinebundles, nbundles) -> (
    -- remember that for each line bundle, we need (#CI-1) * (dim Y + 1) + (dim X + 1) variables
    -- and for each new bundle, we need (dim X + 1) variables
    nvars := nlinebundles * ((#X.CI - 1) * (dim X.Ambient + 1) + dim X + 1) + nbundles * (dim X + 1);
    A := QQ(monoid [VariableBaseName=>opts#Symbol, Variables => 3*nvars+1200]);
    X.cache.ring = A;
    X.cache.cohom = new MutableHashTable;
    X.cache.nextVar = 0; -- before a call to one of the routines here, nextVar == nextFinalVar, but this
      -- is used to 'reserve' variables while determining LES relations.
      -- After a set of cohomology vectors is completed, then they are 'squashed down': the variables appearing
      -- (that are not finalized) are placed starting at nextFinalVar.
    X.cache.nextFinalVar = 0; -- variables before this are finalized, and appear in some cohomology
    )
basicCohomologies = method(Options=>{Limit=>null})
basicCohomologies(CompleteIntersectionInToric) := opts -> (X) -> (
    -- these include all the D_i's, and 00_X too
    -- all that is needed on Y to get Omega_X^1...
    N := #X.CI;
    Y := ambient X;
    dimY := dim Y;
    dimX := dim X;
    if not X.cache.?ring then 
      initializeCohomologies(X, #(rays Y) + 1, 2);
    ndegrees := #(rays Y) - dimY;
    Ds := prepend(Y_0-Y_0, for i from 0 to #(rays Y)-1 list (-Y_i));
    Ds = join(Ds, for h in X.CI list -h);
    Ds = Ds/degree;
    linebundles := collectLineBundles(Y, X.CI, Ds);
    cohomCalg(Y, linebundles); -- actually computes the cohomologies needed over Y, cohoms is a MutableHashTable.
    for d in Ds do cohomologyVector(X, d);
    X.cache.cohom
    )
nextLineBundle = method()
nextLineBundle CompleteIntersectionInToric := (X) -> (
    -- returns a list of (#CI-1) cohom vectors of length: dim Y + 1
    --  and the last one has the same length, but is zero above dim X.
    Y := X.Ambient;
    CI := X.CI;
    nvecs := #CI-1;
    A := X.cache.ring;
    firstvecs := for i from 0 to #CI-2 list (
        result := flatten entries genericMatrix(A, A_(X.cache.nextVar), 1, dim Y + 1);
        X.cache.nextVar = X.cache.nextVar + dim Y + 1;
        result
        );
    lastvec := flatten entries genericMatrix(A, A_(X.cache.nextVar), 1, dim X + 1);
    lastvec = join(lastvec, toList(#CI : 0));
    X.cache.nextVar = X.cache.nextVar + dim X + 1;
    append(firstvecs, lastvec)
    )
nextBundle = method()
nextBundle CompleteIntersectionInToric := (X) -> (
    -- returns a single list, of length: dim X + 1, representing the
    -- cohomoogies of a sheaf or vector bundle on X.
    A := X.cache.ring;
    result := flatten entries genericMatrix(A, A_(X.cache.nextVar), 1, dim X + 1);
    X.cache.nextVar = X.cache.nextVar + dim X + 1;
    result
    )
finalizeVariables = (X, L) -> (
    -- X is a CompleteIntersectionInToric
    -- L: List of Matrix
    --  each one should be a matrix over X.cache.ring.
    -- Returned value: L': a list of matrices over the same ring, 
    --   where the variables have been moved
    --   up to not waste ring indeterminants.
    -- This updates X.cache.nextFinalVar, X.cache.nextVar as well.
    -- YYY
    R := ring L#0;
    firstvar := X.cache.nextFinalVar;
    -- We make a list of the variables that actually occur here
    keepthese := L/support/set//sum;
    keepthese = keepthese - set for i from 0 to firstvar-1 list X.cache.ring_i;
    keepthese = sort toList keepthese;
    -- now we create a ring map that maps these existing variables to their compacted counterpart
    phi := map(R,R,for v from 0 to #keepthese-1 list (keepthese#v => R_(firstvar + v)));
    X.cache.nextVar = X.cache.nextFinalVar = firstvar + #keepthese;
    for m in L list phi m
    )
-- Interface for cohomology in toric complete intersections

cohomologyVector(CompleteIntersectionInToric, List) := (X, degreeD) -> (
    if not X.cache.?cohom then basicCohomologies X;
    if not X.cache.cohom#?degreeD then X.cache.cohom#degreeD = (
        CI := X.CI;
        N := #CI;
        Y := ambient X;
        dimY := dim Y;
        linebundles := collectLineBundles(Y, CI, {degreeD});
        cohoms := cohomCalg(Y, linebundles); -- actually computes the cohomologies needed over Y, cohoms is a MutableHashTable.
        -- Next step is to create the short exact sequences we need.
        H0 := partition(s -> #s, subsets CI);
        H := applyPairs(H0, (r,v) -> (r, if r > 0 then apply(v, v0 -> degreeD - degree sum v0) else apply(v, v0->degreeD)));
        -- H is a hash table: keys are 0..#CI, values: lists of degrees of line bundles at that step in Koszul complex.
        -- Compute all of these cohomologies:
        H1 := applyPairs(H, (r,v) -> (r, sum for v1 in v list cohomologyVector(Y,v1)));
        -- Now create a ring with #CI*(dim Y + 1) number of variables
        -- L contains the ansatz cohomology vectors for the syzygy vector bundles
        L := reverse nextLineBundle X;
        ses := prepend(matrix transpose{H1#N, H1#(N-1), L#(N-1)},
            reverse for i from 0 to N-2 list matrix transpose{L#(i+1), H1#i, L#i});
        J := trim (sum for s in ses list splitLES flatten entries s);
        ses1 := for s in ses list (s%J);
        if debugLevel > 0 then << "exact sequences: " << ses1 << endl;
        M := (last ses1)_{2};
        M = first finalizeVariables(X,{M});
        result := take(flatten entries M, dim X + 1);
        if all(result, x -> liftable(x, ZZ))
        then result = result/(x -> lift(x,ZZ));
        result
        );
    X.cache.cohom#degreeD
    )
cohomologyVector(CompleteIntersectionInToric, ToricDivisor) := (X, D) -> cohomologyVector(X, degree D)
cohomologyVector CompleteIntersectionInToric := (X) -> cohomologyVector(X, 0 * (ambient X)_0)

cohomologyOmega1 = method()
cohomologyOmega1(CompleteIntersectionInToric) := (X) -> (
    -- Need to make 2 ses's, and two cohom vectors for dim X
    -- 
    if not X.cache.?cohom then basicCohomologies X;
    if not X.cache.cohom#?"omega1" then X.cache.cohom#"omega1" = (
        dimX := dim X;
        Y := ambient X;
        CI := X.CI;
        ndegrees := #(rays Y) - dim Y;
        cohoms := basicCohomologies X;
        cohomOOX := cohomologyVector(X, degree (0*Y_0));
        spot1 := sum for i from 0 to #(rays Y)-1 list cohomologyVector(X, - degree Y_i);
        spot2 := sum for h in CI list cohomologyVector(X, - degree h);
        -- note: we computed all the cohomologies above, now we consider vec1, vec2
        -- these are created after the above 4 lines, so that when we call finalizeVariables 
        --in cohomologyVector, they don't conflict with our choice of vec1 and and vec2 variables.
        vec1 := nextBundle X;
        vec2 := nextBundle X;
        ses1 := {
            vec1,
            spot1,
            ndegrees * cohomOOX};
        ses2 := {
            spot2,
            vec1,
            vec2
            };
        sess := {transpose matrix ses1, transpose matrix ses2};
        J1 := ideal(vec2_0 - cohomOOX_1, vec2_(dimX) - cohomOOX_(dimX-1));
        J := trim (J1 + (sum for s in sess list splitLES flatten entries s));
        sess = for s in sess list (s%J);
        if debugLevel > 0 then << "exact sequences: " << sess << endl;
        M := (last sess)_{2};
        M = first finalizeVariables(X,{M});
        vec := flatten entries M;
        if all(vec, x -> liftable(x,ZZ))
        then vec = vec/(x -> lift(x,ZZ));
        vec
        );
    X.cache.cohom#"omega1"
    )

hodgeDiamond = method()
hodgeDiamond CompleteIntersectionInToric := (X) -> (
    -- Assumptions: Y is smooth?
    --  Certainly want: X is smooth.
    -- if dimX <= 3, then we only need Omega1_X
    -- if dimX == 4, then we can either use Omega2_X, or the topological Euler characteristic
    dimX := dim X;
    if dimX >= 4 then <<  "warning: not yet implemented for dimension >= 4, -1's mean not computed" << endl;
    Y := ambient X;
    vec1 := cohomologyOmega1 X;
    vec0 := cohomologyVector(X, degree(0*Y_0));
    matrix for p from 0 to dimX list for q from 0 to dimX list (
        if p == 0 then vec0_q 
        else if q == 0 then vec0_p
        else if p == 1 then vec1_q
        else if q == 1 then vec1_p
        else if p == dimX then vec0_(dimX-q)
        else if q == dimX then vec0_(dimX-p)
        else if p == dimX-1 then vec1_(dimX-q)
        else if q == dimX-1 then vec1_(dimX-p)
        else -1
        )
    )

-- Ds are toric divisors in a toric variety
-- Computes the cohomology vector of the intersection
-- of these divisors.
hodgeVector = method()
hodgeVector List := (Ds) -> (
    if #Ds == 0 then error "expected at least one divisor";
    V := variety Ds#0;
    if not all(Ds, D -> instance(D,ToricDivisor) and variety D === V)
      then error "expected divisors all on the same toric variety";
    D := completeIntersection(V, Ds);
    cohomologyVector(D, degree (V_0-V_0))
    )

-- The following uses cohomology matrix and bases to compute
-- the cohomology dimensions.  This gets exact values, based on specific
-- polynomials in the Cox ring, however it needs specific polynomials.

makeCohomMatrix = method()
makeCohomMatrix(ZZ, NormalToricVariety, RingElement, RingElement) := (i,V,F,G) -> (
   (m1, base1, base2) := cohomologyMatrix(i,V,-degree F, G);
   (m2, base3, base4) := cohomologyMatrix(i,V,-degree G, F);
   if base2 != base4 then error "hmm, expected bases to be identical for columns";
   (matrix m1 || matrix m2, join(base1, base3), base2)
   )
cohomVectorOfCodim2 = method()
cohomVectorOfCodim2(NormalToricVariety, RingElement, RingElement) := (V,F,G) -> (
    maps := for i from 0 to dim V list makeCohomMatrix(i,V,F,G);
    rks := for i from 0 to dim V list (#maps_i_1, rank maps_i_0, #maps_i_2);
    if rks_0 != (0,0,0) then << "warning: expected H^0 part to be all zero" << endl;
    if rks_1_2 != rks_1_1 then << "warning: expected H^1 map to be an inclusion" << endl;
    if rks_(dim V)_1 != rks_(dim V)_0 then << "warning: expected H^" << dim V << " map to be a surjection" << endl;
    for i from 0 to dim V - 2 list (
        rks_(i+1)_0 + rks_(i+2)_2 - rks_(i+1)_1 - rks_(i+2)_1 + if i == 0 then 1 else 0
        )
    )
hodgeVector(List, List) := (Ds, Fs) -> (
    -- Ds: list of (effective) divisor classes on a toric variety
    -- Fs: list of polynomials, one for each class.
    -- currently: these are limited to #Ds == #Fs == 2
    -- invariant: degree Ds_i == degree Fs_i, all i.
    if #Ds =!= 2 or #Fs =!= 2 then error "expected codimension 2 complete intersection in a toric variety";
    if degree Ds_0 =!= degree Fs_0 or degree Ds_1 =!= degree Fs_1 then
      error "expected polynomials to be in the given divisor classes";
    V := variety Ds_0;
    cohomVectorOfCodim2(V, Fs_0, Fs_1)
    )
