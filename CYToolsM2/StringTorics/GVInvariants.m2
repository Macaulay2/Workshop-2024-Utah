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




-------------------------------------------------------------------------
-- some tests -----------------------------------------------------------


///
  restart
  debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
  DB3 = "../Databases/cys-ntfe-h11-3.dbm"
  DB3 = "./StringTorics/Databases/cys-ntfe-h11-3.dbm"
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ);
  (Qs, Xs) = readCYDatabase(DB3, Ring => RZ);

  moris = for k in sort keys Xs list (
      X = Xs#k;
      if not isFavorable X then continue; -- only handles favorable polytopes and CY3's.
      classifyExtremalCurves(X, Verbose => 2)
      );

  moris/keys/sort//unique
      degvec = heft X;
      deglimit = max for c in toricMoriConeCap X list 3 * dotProduct(degvec, c);
      << "degree limit for " << k << " is " << deglimit << endl;
      gvX = gvInvariants(X, DegreeLimit => deglimit);
      moriClass = sort for c in toricMoriConeCap X list
          classifyExtremalCurve(gvX, c, deglimit, degvec);
      print netList (ans := prepend(k, moriClass));
      moriClass
      )     

  restart
  debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
  DB4 = "../Databases/cys-ntfe-h11-4.dbm"
  RZ = ZZ[a,b,c,d]
  RQ = QQ (monoid RZ);
  (Qs, Xs) = readCYDatabase(DB4, Ring => RZ);
  
  X = Xs#(137,0)
  for k in sort keys Xs list (
  deglimit = 14
  degvec = heft X
  gvX = gvInvariants(X, DegreeLimit => deglimit);

  netList for c in toricMoriConeCap X list
    sort classifyExtremalCurve(gvX, c, deglimit, degvec)

  for k in sort keys Xs list (
      X = Xs#k;
      if not isFavorable X then continue; -- only handles favorable polytopes and CY3's.
      deglimit = 14;
      degvec = heft X;
      gvX = gvInvariants(X, DegreeLimit => deglimit);
      moriClass = sort for c in toricMoriConeCap X list
          classifyExtremalCurve(gvX, c, deglimit, degvec);
      print moriClass;
      moriClass
      )     
  moricap = toricMoriConeCap X
  for c in toricMoriConeCap X list
    c => gvRay(gvX, c, 22, heft X)
  netList oo
  

///

///
  restart
  debug needsPackage "StringTorics"
  DB4 = "../Databases/cys-ntfe-h11-4.dbm"
  RZ = ZZ[a,b,c,d]
  RQ = QQ (monoid RZ);
  (Qs, Xs) = readCYDatabase(DB4, Ring => RZ);

  moris = for k in sort keys Xs list elapsedTime (
      X = Xs#k;
      if not isFavorable X then continue; -- only handles favorable polytopes and CY3's.
      mori1 = classifyExtremalCurves(X, Verbose => 2);
      print netList {mori1};
      mori1
      );


  X = Xs#(34,0)
  assert isFavorable X
  classifyExtremalCurves(Xs#(34,0), Verbose => 2)
  classifyExtremalCurves(Xs#(35,0), Verbose => 2)
  mori = toricMoriConeCap X;
  deglimit = 3 * gvTopMoriConeCapDegree X;
  deglimit
  degvec = heft X;
  gvX = gvInvariants(X, DegreeLimit => deglimit)
  partition(c -> classifyExtremalCurve(gvX, c, deglimit, degvec), mori)
///

TEST ///
-- Good test, TODO: place this back in once GVinvariants are working again.
-*
  restart
  needsPackage "StringTorics"
*-
  -- hh^(1,1) == 3
  -- label (5,0)
  rys = {{-1, -1, 1, 0}, {-1, -1, 1, 1}, {-1, -1, 2, 1}, {-1, 3, -2, -1}, {1, -1, 0, 0}, {2, -1, 0, 0}, {-1, 1, 0, 0}}
  cones4 = {{0, 1, 2, 4}, {0, 1, 2, 6}, {0, 1, 3, 4}, {0, 1, 3, 6}, {0, 2, 4, 5}, {0, 2, 5, 6}, {0, 3, 4, 5}, {0, 3, 5, 6}, {1, 2, 4, 5}, {1, 2, 5, 6}, {1, 3, 4, 5}, {1, 3, 5, 6}}
  Q = cyPolytope(rys, ID => 5)
  label Q
  assert(rays Q == rys)
  assert(hh^(1,1) Q == 3)
  assert(hh^(1,2) Q == 57)

  RZ = ZZ[a,b,c]
  Xs = findAllCYs(Q, Ring => RZ)
  X = Xs#0
  label X == (5,0)
  assert(max X === cones4)
  assert(rays X === rys)

  -- OK, now we are ready to test functions in this file
  C = toricMoriCone X;
  heft1 = sum entries transpose rays dualCone C

  rays toricMoriCone X
  morirays = toricMoriConeCap X
  assert(sort morirays === sort {{0, 0, 1}, {0, 1, 0}, {1, -1, 0}, {1, 0, -1}})
  assert(morirays === {{0, 0, 1}, {0, 1, 0}, {1, -1, 0}, {1, 0, -1}}) -- not required.
  heft X-- why not {2,1,1}??

  gv = gvInvariants(X, DegreeLimit => 15)
  (keys oo)/toList//matrix//transpose//posHull//rays
  gv#(0,1,1)
  
  assert(extremalRayGVs(X, {0,0,1}, Limit => 5) == {60, 0, 0, 0, 0})
  assert(extremalRayGVs(X, {0,1,0}, Limit => 5) == {6, 0, 0, 0, 0})
  assert(extremalRayGVs(X, {1,-1,0}, Limit => 5) == {252, -9252, 848628, -114265008, 18958064400})
  assert(extremalRayGVs(X, {1,0,-1}, Limit => 5) == {56, -272, 3240, -58432, 1303840})

  assert(extremalRayGVs(X, {1,-1,-3}, Heft => {5,1,1}, Limit => 5) == {0,0,0,0,0}) -- all zeros...

  -- Note: if the curve class is in the interior, the answer from extremalRayGVs is NOT correct.
  -- Also: a zero vector means either that the curve class is not effective on X,
  --  or that there is an elliptic ruled surface of some sort on X.  If one jiggles this X,
  --  such surfaces go away, ie. in general moduli, an extrmeal ray has zero GV's iff
  --  every class along that curve ray is NOT effective.

  --------------
  -- gvCone ----
  -- this is the cone of all non-zero GV curves.  This either matches the Mori cone, or is possibly
  -- a subset (as some rays can contain effective curves, all of whose GV invariants (in the ray)
  -- are zero.
  --------------
  assert(set entries transpose rays gvCone X == set toricMoriConeCap X) -- for this example,
  -- the intersection of the Mori cones for the various simplicial
  -- toric varieties (which is always contains the Mori cone of X),
  -- all has non-zero GV invariants, so if one believes the GV
  -- computation, this is exactly the Mori cone of X, so we know the
  -- nef/Kahler cone for X too

--  gvInvariantsAndCone(X, 15)-- doesn't work at the momemnt... BUG -- remove?

  partitionGVConeByGV X -- not completely what we want?

  classifyExtremalCurve extremalRayGVs(X, {0,0,1}, Limit => 5)
  classifyExtremalCurve extremalRayGVs(X, {0,1,0}, Limit => 5)
  classifyExtremalCurve extremalRayGVs(X, {1,-1,0}, Limit => 5)
  classifyExtremalCurve extremalRayGVs(X, {1,0,-1}, Limit => 5)

  classifyExtremalCurves X
  heft X

  debug needsPackage "StringTorics" -- for gvTopMoriConeCapDegree
  assert(gvTopMoriConeCapDegree X == 2) -- not exported
///


///
-- Tests of this code, 15 Jan 2023. Removed from tests, since it used created databases...
-*
  restart
  needsPackage "StringTorics"
*-
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  DBNAME = "./Databases/cys-ntfe-h11-5.dbm"
  DBNAME = "./StringTorics/Databases/cys-ntfe-h11-5.dbm"
  RZ = ZZ[a,b,c,d,e]
  RQ = QQ (monoid RZ);
--  needs "../FindEquivalence.m2"
  (A,phi) = genericLinearMap RQ
  elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);
  REPS = value get "./Analysis/inequiv-reps-h11-5"
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
  
  -- We can use this collection to test IntegerEquivalences package:
  -- Goal: for each set of labels REPS#i, determine if these are the same topology or different.
  -- Note: some are very easy, some I can't yet do.

  -- example: REPS#87 {(1835, 0), (1864, 0), (1876, 0)}
  -- example: REPS#222
  (X1, X2, X3) = REPS#87/(lab -> Xs#lab)//toSequence
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  (L3, F3) = (c2Form X3, cubicForm X3)

  (X1, X2) = REPS#222/(lab -> Xs#lab)//toSequence
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  
  hessianMatches(F1, F2)   | singularPointMatches(sub(F1, RQ), sub(F2, RQ)) | matchingData{L1 => L2, F1 => F2}
  selectLinear oo | matchingData{L1 => L2, F1 => F2}
  tryEquivalences(oo, RQ, (A, phi))
  
  GV1 = classifyExtremalCurves X1
  GV2 = classifyExtremalCurves X2
  netList {GV1, GV2}

  (X1, X2) = REPS#44/(lab -> Xs#lab)//toSequence
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  
  hessianMatches(F1, F2)   | singularPointMatches(sub(F1, RQ), sub(F2, RQ)) | matchingData{L1 => L2, F1 => F2}
  selectLinear oo | matchingData{L1 => L2, F1 => F2}
  tryEquivalences(oo, RQ, (A, phi))
  
  GV1 = classifyExtremalCurves X1
  GV2 = classifyExtremalCurves X2
  netList {GV1, GV2}


  (X1, X2) = REPS#45/(lab -> Xs#lab)//toSequence
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  
  hessianMatches(F1, F2)   | singularPointMatches(sub(F1, RQ), sub(F2, RQ)) | matchingData{L1 => L2, F1 => F2}
  selectLinear oo | matchingData{L1 => L2, F1 => F2}
  tryEquivalences(oo, RQ, (A, phi))
  
  GV1 = classifyExtremalCurves X1
  GV2 = classifyExtremalCurves X2
  netList {GV1, GV2}

  (X1, X2) = REPS#46/(lab -> Xs#lab)//toSequence
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  
  hessianMatches(F1, F2)   | singularPointMatches(sub(F1, RQ), sub(F2, RQ)) | matchingData{L1 => L2, F1 => F2}
  selectLinear oo | matchingData{L1 => L2, F1 => F2}
  tryEquivalences(oo, RQ, (A, phi))
  
  GV1 = classifyExtremalCurves X1
  GV2 = classifyExtremalCurves X2
  netList {GV1, GV2}

  matchingData{transpose matrix GV1#{"FLOP", {1,0,0,0}} => transpose matrix GV2#{"FLOP", {1,0,0,0}}}
  matchingData{{Permutations, transpose matrix GV1#{"FLOP", {3,0,0,0}}, transpose matrix GV2#{"FLOP", {1,0,0,0}}}

  -- hessians which should be straightforward
  set1 = SETS#{{{5,1}}, {{1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}}}  
  elapsedTime for a in take(set1,6) list (
      if #REPS#a != 2 then continue;
      << "DOING " << a << endl;
      (X1, X2) = REPS#a/(lab -> Xs#lab)//toSequence;
      (L1, F1) = (c2Form X1, cubicForm X1);
      (L2, F2) = (c2Form X2, cubicForm X2);
      md := hessianMatches(F1, F2) | {L1 => L2} | {F1 => F2};
      md1 := selectLinear md;
      ans := tryEquivalences(md1 | {F1 => F2}, RQ, (A, phi));
      --if ans#0 === INCONSISTENT then {{first a}, {last a}} else if ans#0 === CONSISTENT then {first a, last a, ans#1} else {a, ans}
      print (a => ans);
      ans
      )
  -- if a pair are the same, make a list {first, second, matrix}
  -- if a pair is not the same, make {{first}, {second}}      
  
  set1 = SETS#{{{5,1}}, {{1, 1}, {1, 1}, {1, 1}, {2, 1}}} -- {33, 58} both distinct.
  set1 = SETS#{{{4,1},{4,1}}, {{1, 1}, {1, 1}, {1, 1}, {2, 1}}} -- {133, 170, 209} distinct pairs
  set1 = SETS#{{{5,1}}, {{1, 1}, {1, 1}, {3, 1}}} -- 
  elapsedTime for a in take(set1,6) list (
      if #REPS#a != 2 then (
          << "DEFERRING " << a << " " << REPS#a << endl;
          continue;
          );
      << "DOING " << a << " " << REPS#a << endl;
      (X1, X2) = REPS#a/(lab -> Xs#lab)//toSequence;
      (L1, F1) = (c2Form X1, cubicForm X1);
      (L2, F2) = (c2Form X2, cubicForm X2);
      md := hessianMatches(F1, F2) | {L1 => L2} | {F1 => F2};
      md1 := selectLinear md;
      ans := tryEquivalences(md1 | {F1 => F2}, RQ, (A, phi));
      --if ans#0 === INCONSISTENT then {{first a}, {last a}} else if ans#0 === CONSISTENT then {first a, last a, ans#1} else {a, ans}
      print (a => ans);
      ans
      )
  
///
