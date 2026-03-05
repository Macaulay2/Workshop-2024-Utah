---------------------------------------
-- gvInvariants
-- Gopakumar-Vafa invariants (similar to Gromov-Witten invariants)
-- Contains code to call the external C++ program `computeGV` from CYTools
-- This requires computing some information first (intersection numbers, mori cone cap, etc).
---------------------------------------

toricMoriCone(NormalToricVariety, List) := Cone => (V, basisIndices) -> (
    IV := intersectionRing (abstractVariety V);
    Cs := matrix for x in orbits(V, 1) list (
        c := product(x, i -> IV_i);
        for j in basisIndices list integral(c * IV_j)
        );
    posHull transpose lift(Cs, QQ) -- TODO: lift to ZZ?
    )

toricMoriCone CalabiYauInToric := Cone => X -> (
    -- TODO: handle toric mori cones of non-favorables
    if isFavorable X then toricMoriCone(ambient X, basisIndices X)
    )

hilbertBasisGenerators = method()
hilbertBasisGenerators Cone := List => C -> (
    for x in hilbertBasis C list flatten entries x
    )

-- This function returns a very large heft vector.  Not so good!
-- TODO: this appears to be computing the toricMoriCone, which it is not using!?
heft CalabiYauInToric := List => X -> (
    C := toricMoriCone X;
    heft1 := sum entries transpose rays dualCone C;
    Ccap := posHull transpose matrix toricMoriConeCap X;
    heft2 := sum entries transpose rays dualCone Ccap;
    heft2
    )

-- TODO: use findProgram/runProgram methods in M2 to handle access to computeGV.
gvInvariants = method(Options => {
    Mori => null, -- null means: compute rays of the Mori cone of V (in ZZ^(h11))
    Heft => null, -- null means: compute it
    DegreeLimit => infinity,
    Precision => 150,
    FilePrefix => "foo",
    Executable => (options StringTorics).Configuration#"computeGV",
    KeepFiles => true
    })

-- The function to write the data needed by the computeGV program
gvInput = (moriGenerators, heftval, GLSM, intersectionnums, degreelimit, prec) -> (
    -- moriGenerators: list of lists. Hilbert basis of the cone of 
    --   irreducible curves induced from the toric variety.
    -- heftval: list of ints
    -- GLSM: list of list of ints
    -- intersectionnums: list of triples of ints
    -- degreelimit: infinity or positive integer
    -- prec: positive integer
    str1 := toString moriGenerators;
    str3 := toString heftval;
    str4 := toString GLSM;
    str5 := toString intersectionnums;
    str6 := toString ({
            if degreelimit === infinity then -1 else degreelimit, 
            prec,
            0,
            300000
            });
    concatenate between("\n", {str1, toString {}, str3, str4, toString {}, str5, str6})
    )

filenameCounter := 0; -- TODO: not used? or change to use it?

-- TODO: use findProgram/runProgram to get this...
-- TODO: Also: there should be one function which calls computeGV.
-- Here there are two...
gvInvariants(NormalToricVariety, List) := HashTable => opts -> (V, basisIndices) -> (
    -- Compute intersection numbers for X in V (using this basis)
    -- Compute mori cone (if needed) (?? requires basis too...)
    -- Compute a vector which dots positively with all these generators.
    -- Then write the file
    -- Execute the command
    -- Read the results, and return them
    intersectionnums := for t in intersectionNumbersOfCY(V, basisIndices) list append(t#0, t#1);
    -- X := completeIntersection(V, {-toricDivisor V});
    -- Xa := abstractVariety(X, base());
    -- IX := intersectionRing Xa;
    -- intersectionnums := for x in pairs intersectionNumbers(IX, basisIndices) list append(x#0, x#1);
    -- H := hashTable for i from 0 to #basisIndices-1 list basisIndices#i => i;
    -- intersectionnums := for x in pairs CY3NonzeroMultiplicities V list (
    --     if isSubset(x#0, basisIndices) then
    --         append(sort for a in x#0 list H#a, x#1)
    --     else 
    --         continue
    --     );
    mori := if opts.Mori =!= null then 
                opts.Mori 
            else 
                hilbertBasisGenerators toricMoriCone(V, basisIndices);
    heft := if opts.Heft =!= null then opts.Heft else (
      sum entries transpose rays dualCone posHull transpose matrix mori
      );
    -- OK, now we have computed everything we need.  Write it to a file
    infile := opts.FilePrefix | "-input";
    outfile := opts.FilePrefix | "-output";
    infile << gvInput(mori, heft, transpose degrees ring V, intersectionnums,
        opts.DegreeLimit, opts.Precision) << close;
    inputLine := opts.Executable | " <" | infile | " >" | outfile;
    print inputLine;
    run inputLine;
    -- Get the output, package as a hash table
    (lines get outfile)/value//hashTable
    )

gvInvariants CalabiYauInToric := HashTable => opts -> X -> (
    if not isFavorable X then return null;
    intersectionnums := for t in intersectionNumbers X list append(t#0, t#1);
    -- mori := if opts.Mori =!= null then 
    --             opts.Mori 
    --         else 
    --             hilbertBasisGenerators toricMoriCone(ambient X, basisIndices X);
    mori := if opts.Mori =!= null then 
                opts.Mori 
            else 
                hilbertBasisGenerators posHull transpose matrix toricMoriConeCap X;
    heft := if opts.Heft =!= null then opts.Heft else (
      sum entries transpose rays dualCone posHull transpose matrix mori
      );
    -- OK, now we have computed everything we need.  Write it to a file
    infile := temporaryFileName(); -- opts.FilePrefix | "-input" 
    outfile := temporaryFileName(); -- opts.FilePrefix | "-output" | filenameCounter;
    infile << gvInput(mori, heft, transpose degrees X, intersectionnums,
        opts.DegreeLimit, opts.Precision) << close;
    inputLine := opts.Executable | " <" | infile | " >" | outfile;
    print inputLine;
    run inputLine; -- TODO: run this as a program and if it crashes, return something reasonable.
    -- Get the output, package as a hash table
    contents := get outfile;
    if #contents == 0 then return null;
    (lines contents)/value//hashTable
    )

-- Not used anymore??  See `extremalRayGVs`
gvRay = method(Options => options gvInvariants)
gvRay(CalabiYauInToric, List) := HashTable => opts -> (X, curveClass) -> (
    -- This doesn't seem to be correct
    if not isFavorable X then return null;
    degvec := heft X;
    grad := dotProduct(degvec, curveClass);
    << "using DegreeLimit: " << 4*grad << endl;
    return gvInvariants(X,Mori => {curveClass}, DegreeLimit => 4 * grad, Heft => {0,1,0,0})
    )

gvCone = method(Options => options gvInvariants)
gvCone CalabiYauInToric := Cone => opts -> X -> (
    if not isFavorable X then return null;
    gv := gvInvariants(X, opts);
    if gv === null then return null;
    posHull transpose matrix ((keys gv)/toList)
    )

gvInvariantsAndCone = method(Options => options gvInvariants)
gvInvariantsAndCone(CalabiYauInToric, ZZ) := Sequence => opts -> (X, D) -> (
    -- D is the degree bound to start with.  We could start with 5, or DegreeLimit/2 or DegreeLimit/4, or ...
    if not isFavorable X then return null;
    degvec := heft X;
    gv := gvInvariants(X, opts);
    if gv === null then return null;
    keysgv := keys gv;
    H := hashTable for k in keysgv list k => dotProduct(k, degvec);
    firstSet := select(keys H, k -> H#k <= D);
    if debugLevel > 0 then << "The number of curves in the first set: " << #firstSet << endl;
    C := posHull transpose matrix (firstSet);
    Cdual := dualCone C;
    HC := transpose rays Cdual;
    curves := for k in keys H list if H#k > D then transpose matrix {k} else continue;
    set2 := select(curves, c  -> any(flatten entries (HC * c), a -> a < 0));
    if debugLevel > 0 then << "The number of curves not in the first cone: " << #set2 << endl;
    C2 := if #set2 == 0 then C else posHull (rays C | matrix{set2});
    if debugLevel > 0 and #set2 == 0 then (
        << "CY " << label X << " C = " << rays C  << endl
        )
    else
        << "*differs* CY " << label X << " C1 = " << rays C << " and C2 = " << rays C2 << endl;
    (gv, C2)
    )

partitionGVConeByGV = method(Options => options gvInvariants)
partitionGVConeByGV CalabiYauInToric := HashTable => opts -> X -> (
    -- return null if we cannot computr GV invariants (i.e. if non-favorable).
    if not isFavorable X then return null;
    gv := gvInvariants(X, opts); -- TODO: stash this?
    if gv === null then return null;
    C := posHull transpose matrix ((keys gv)/toList);
    gvX := entries transpose rays C;
    partition(f -> if gv#?(toSequence f) then gv#(toSequence f) else 0, gvX)
    )

partitionGVConeByGV(CYToolsCY3, ZZ) := HashTable => opts -> (X, D) -> (
    -- return null if we cannot computr GV invariants (i.e. if non-favorable).
    (gv, C) := gvInvariantsAndCone(X, D, opts);
    gvX := entries transpose rays C;
    partition(f -> if gv#?f then gv#f else 0, gvX)
    )

-- TODO: move to Topology.m2? file?
findLinearMaps = method()
findLinearMaps(HashTable, HashTable) := List => (gv1, gv2) -> (
    -- gv1, gv2: result of partitionGVConeByGV
    if sort keys gv1 =!= sort keys gv2 then return {};
    for k in keys gv1 do if #gv1#k =!= #gv2#k then return {};
    for k in keys gv1 do if #gv1#k >= 7 then return {}; -- do not waste time (1) trying to separate these?
    n := # (first values gv1)_0; -- we should check if all the values are lists of integers of this size.
    t := symbol t;
    T := QQ[t_(1,1)..t_(n,n)];
    M := genericMatrix(T, n, n);
    -- now we make the ideals for each key, and each permutation.
    ids := for k in keys gv1 list (
        perms := permutations(#gv1#k);
        mat1 := transpose matrix gv1#k;
        mat2 := transpose matrix gv2#k;
        for p in perms list (
            I := trim ideal (M * mat1 - mat2_p); 
            if I == 1 then continue else I
            )
        );
    topval := ids/(x -> #x - 1);
    zeroval := ids/(x -> 0);
    fullIdeals := for a in zeroval .. topval list (
        J := trim sum for i from 0 to #ids-1 list ids#i#(a#i);
        if J == 1 then continue else J
        );
    Ms := for i in fullIdeals list M % i;
    --newMs := select(Ms, m -> (d := det m; d == 1 or d == -1));
    --if any(newMs, m -> support m =!= {}) then << "some M is not reduced to a constant" << endl;
    Ms
    )

gvRay(HashTable, List, ZZ, List) := opts -> (GVHash, C, deglimit, degvector) -> (
    contentC := gcd C;
    if contentC =!= 1 then C = C // contentC;
    d := dotProduct(degvector, C);
    rayC := for i from 1 to floor(deglimit/d) list (
        Cseq := toSequence(i*C);
        if GVHash#?Cseq then GVHash#Cseq else 0
        );
    rayC
    )

gvRay(HashTable, List, ZZ, List) := opts -> (GVHash, C, deglimit, degvector) -> (
    contentC := gcd C;
    if contentC =!= 1 then C = C // contentC;
    d := dotProduct(degvector, C);
    rayC := for i from 1 to floor(deglimit/d) list (
        iC := i*C;
        if GVHash#?iC then GVHash#iC else 0
        );
    rayC
    )

count = 0; -- used to give a unique index to each ZERO ray.

classifyExtremalCurve = method()

classifyExtremalCurve List := rayC -> (
    -- rayC: a list of the gv invariants along the ray of a toric mori cone extremal curve.
    -- these are either extremal on the CY, or not effective on the CY.
    if #rayC <= 2 then return {"OTHER", rayC};
    if all(2..#rayC-1, i -> rayC#i == 0) then (
        -- only first two, possibly, are non-zero.
        if rayC#0 == 0 and rayC#1 == 0 then (count=count+1; return {"ZERO", count});
        if rayC#0 == -2 or rayC#1 == -2 then return {"TYPEIII0", rayC};
        if rayC#0 >= 0 and rayC#1 >= 0 then return {"FLOP", rayC}; -- these could be type IIIg as well?
        if rayC#0 < 0 or rayC#1 < 0 then return {"TYPEIIIg", rayC}; -- 
        )
    else return {"TYPEII", rayC}
    )

classifyExtremalCurve(HashTable, List, ZZ, List) := (GVHash, C, deglimit, degvector) -> (
    rayC := gvRay(GVHash, C, deglimit, degvector);
    if #rayC <= 2 then return {"OTHER", rayC};
    if all(2..#rayC-1, i -> rayC#i == 0) then (
        -- only first two, possibly, are non-zero.
        if rayC#0 == 0 and rayC#1 == 0 then (count=count+1; return {"ZERO", count});
        if rayC#0 == -2 or rayC#1 == -2 then return {"TYPEIII0", {rayC#0, rayC#1}};
        if rayC#0 >= 0 and rayC#1 >= 0 then return {"FLOP", {rayC#0, rayC#1}};
        if rayC#0 < 0 or rayC#1 < 0 then return {"TYPEIIIg", {rayC#0, rayC#1}};
        )
    else return {"TYPEII", {rayC#0, rayC#1, rayC#2, "..."}}
    )

gvTopMoriConeCapDegree = method()
gvTopMoriConeCapDegree CalabiYauInToric := X -> (
    if not isFavorable X then error "expected a favorable polytope";
    degvec := heft X;
    max for c in toricMoriConeCap X list dotProduct(degvec, c)
    )

classifyExtremalCurves = method(Options => {
        Verbose => 0,
        DegreeLimit => null,
        MoriHilbertGens => null
        })
classifyExtremalCurves(HashTable, List, ZZ, List) := (GVHash, Cs, deglimit, degvector) -> (
    partition(c -> classifyExtremalCurve(GVHash, c, deglimit, degvector), Cs)
    )

classifyExtremalCurves CalabiYauInToric := opts -> X -> (
    if not isFavorable X then error "expected favorable CY3-fold";
    mori := if opts.MoriHilbertGens === null then toricMoriConeCap X else opts.MoriHilbertGens;
    deglimit := 3 * gvTopMoriConeCapDegree X;
    if opts.Verbose > 1 then << "*** mori cone cap degree limit is " << deglimit << " ***" << endl;
    degvec := if opts.DegreeLimit === null then heft X else opts.DegreeLimit;
    gvX := gvInvariants(X, DegreeLimit => deglimit);
    partition(c -> classifyExtremalCurve(gvX, c, deglimit, degvec), mori)
    )

-- Good one here, I think.
extremalRayGVs = method(Options => {Limit => 4, Heft => null})
extremalRayGVs(CalabiYauInToric, List) := opts -> (X, curveClass) -> (
    if not isFavorable X then return null; -- later, maybe we can modify this...
    degvec := if opts.Heft =!= null then opts.Heft else heft X;
    deglimit := opts.Limit * dotProduct(degvec, curveClass);
    gvHash := gvInvariants(X,Mori => {curveClass}, DegreeLimit => deglimit, Heft => degvec);
    for i from 1 to opts.Limit list (
        c := toSequence(i * curveClass);
        if gvHash#?c then gvHash#c else 0
        )
    )

classifyExtremalCurves CalabiYauInToric := opts -> X -> (
    if not isFavorable X then error "expected favorable CY3-fold";
    if not X.cache#?"toric mori cone gvs" then (
        mori := toricMoriConeCap X;
        val := hashTable for c in mori list (
            gvs := extremalRayGVs(X, c, Limit => 4);
            c => classifyExtremalCurve gvs
            );
        X.cache#"toric mori cone gvs" = partition(c -> val#c, mori)
        );
    X.cache#"toric mori cone gvs"
    )

extremalCurveInvariant = method()
extremalCurveInvariant CalabiYauInToric := X -> (
    gv := classifyExtremalCurves X;
    sort for a in pairs gv list {a#0, #a#1}
    )

-----------------------------------------
-- GV by rays ---------------------------
-- this is a better data format for our intended usage.
-- it contains about the same info, but it also retains information
-- about degree limit, and degree vector used...
-- Possibly: want to allow to merge two tables, using different degree vectors
-- (originally in the Flops.m2 file)
-----------------------------------------

--------------------
-- GVInvariants ----
--------------------
-- Should this move to GVInvariants.m2 in StringTorics?
-- Probably...
gvByRay = method()
-- gvByRay: Take a "standard" hash table of GV invariants for various curve classes,
--  and data used to compute that, and create a sometimes more useful version of the data:
--  a hash table whose keys are primitive curve classes, and whose values are list of GV invariants
--  for a number of multiples of the primitive curve, up to some degree bound.
gvByRay(HashTable, ZZ, List) := HashTable => (GVHash, deglimit, degvector) -> (
      primcurve := curve -> (
          i := position(curve, a -> a != 0);
          c := curve // gcd curve;
          mult := curve#i // c#i;
          {c, mult}
          );
      curves := for k in keys GVHash list join(primcurve k, {GVHash#k});
      H := partition(val -> val#0, curves);
      H2 := hashTable for kv in pairs H list (
          val := kv#0 => kv#1/(x -> {x#1, x#2});
          val);
      -- now we change the values to be a list og gv's.
      toRayList := (primcurve, listOfPairs) -> (
          Hmult2gv := hashTable listOfPairs;
          d := max (listOfPairs/first);
          deg := dotProduct(degvector, primcurve);
          maxmult := max(d, deglimit//deg);
          for i from 1 to maxmult list if Hmult2gv#?i then Hmult2gv#i else 0
          );
      result := hashTable for kv in pairs H2 list (
          val := toRayList(kv#0, kv#1);
          kv#0 => val
          --if #val >= 2 then kv#0 => val else continue -- this misses curves that we might need for mori cone!
          );
      result
      )

GVTable = new Type of HashTable
expression GVTable := GVT -> (
    expression GVT.GVRays
    )

gvTable = method(Options => {DegreeLimit => 20, Heft => {}})

-- GVInvariants is a hash table constructed with gvInvariants function.
gvTable HashTable := GVTable => opts -> GVinvariants -> (
    if opts.DegreeLimit === null then error "need to provide DegreeLimit value, e.g. DegreeLimit => 20";
    if opts.Heft === null then error "need to provide Heft list, e.g. Heft => {1,3,6,2}";
    -- now we take the input hashtable of curveclass => GVinvariant
    -- and change it to primitiveCurve C => {list of GV invariants of C, 2*C, 3*C, ...}
    -- where the values go up max i*C whose heft vector value is <= degreelimit.
    gvs := gvByRay(GVinvariants, opts.DegreeLimit, opts.Heft);
    result := new GVTable from {
        symbol cache => new CacheTable,
        symbol DegreeLimit => opts.DegreeLimit,
        symbol Heft => opts.Heft,
        symbol GVRays => gvs
        };
    result
    )
gvTable CalabiYauInToric := GVTable => opts -> X -> (
    if opts.DegreeLimit === null then error "need to provide DegreeLimit value, e.g. DegreeLimit => 20";
    degvector := heft X;
    gvH := gvInvariants(X, DegreeLimit => opts.DegreeLimit);
    gvH2 := hashTable for kv in pairs gvH list {toList kv#0, kv#1};
    gvTable(gvH2, DegreeLimit => opts.DegreeLimit, Heft => heft X)
    )

gvRays = method()
gvRays GVTable := gvTable -> gvTable.GVRays

moriCone = method()
moriCone(GVTable, List) := Cone => (gvTable, negatedCurves) -> (
    -- We take the primitive curves in the table, negate the ones that need negating,
    -- and make the cone of all these
    primcurves := keys gvRays gvTable;
    negatedCurves = set negatedCurves;
    curves := matrix transpose for c in primcurves list if negatedCurves#?c then -c else c;
    posHull curves
    )


-- TODO: use the GV code to do these rays directly?  Is that possible?
gvRay(GVTable, List) := opts -> (gvTable, curve) -> (
    if (gvRays gvTable)#?curve then (gvRays gvTable)#curve else {}
    )

isNilpotent = method()
isNilpotent(GVTable, List) := (gvTable, curve) -> (
    gvrays := gvRays gvTable;
    if not gvrays#?curve then return false;
    thisray := gvrays#curve;
    (#thisray >= 4 and thisray#-1 == 0 and thisray#-2 == 0 )
      or (#thisray <= 3 and thisray#-1 == 0)
    )

-- Database-dependent GV tests and scratch code moved to ScratchGVInvariants.m2
