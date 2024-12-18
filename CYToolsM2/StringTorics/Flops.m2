-- Notes to self:
--  Naomi's code is   /Users/mike/src/stringtorics/naomi-flop-code/Flop Code:
--  Jakob's code is in 

-- Code for constructing potential flops of a CY Hypersurface, or something 
-- constructed from that bvia a sequence of flops.

-- given a curve, gv invariant, topology.  Return the new topology.

-- Determine the types of arguments for:

-- find_nilpotent
-- find_nilpotent_outside_inf
-- is_symmetric_flop
-- find_all_flops

-- what about:
--   toric_curves.compute

--   two_face_triags.all_two_face_triangulations(p)

-- determine if a curve is gv nilpotent
-- determine if a nilpotent ray is "outside the infinity cone"
-- perform a flop


---- given:
----  Input: CY3, as given  by h11, h12, c2, cubic, and also non-zero GV classes up to some cutoff.
----  Input: a curve class, representing a flop.
----  Output: a new CY3
----    negate the curve class, take all effective curves other than that.  What if a curve class has 0 gv's?
----    use the gv of the flopped curve.

debug needsPackage "StringTorics"

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
          Hmult2gv = hashTable listOfPairs;
          d := max (listOfPairs/first);
          deg := dotProduct(degvector, primcurve);
          maxmult := max(d, deglimit//deg);
          for i from 1 to maxmult list if Hmult2gv#?i then Hmult2gv#i else 0
          );
      hashTable for kv in pairs H2 list (
          val := toRayList(kv#0, kv#1);
          if #val >= 2 then kv#0 => val else continue
          )
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
gvRay(GVTable, List) := opts -> (gvTable, curve) -> (gvRays gvTable)#curve

isNilpotent = method()
isNilpotent(GVTable, List) := (gvTable, curve) -> (
    gvrays := gvRays gvTable;
    if not gvrays#?curve then return false;
    thisray := gvrays#curve;
    (#thisray >= 4 and thisray#-1 == 0 and thisray#-2 == 0 )
      or (#thisray <= 3 and thisray#-1 == 0)
    )
--------------------------------------------
--------------------------------------------

-- this is the version of the GV hash table we will place into a "CY3" object:



-- CY3 -- a type which encapsulates the topology (h11, h12, cubic form, c2Form), as well
-- as some GV invariant information, as well as some information about flops.

-- remove this commented out code block.
-- CY3 = new Type of HashTable
-- makeCY3 = method(Options => {GV => null, -- of the original, toric, hypersurface
--         NegatedCurves => null,
--         Heft => null, 
--         DegreeLimit => null, 
--         Label => null, 
--         MoriCone => null
--         })
-- makeCY3(ZZ, ZZ, RingElement, RingElement) := opts -> (h11val, h12val, L, F) -> (
--     X := new CY3 from {
--         cache => new CacheTable,
--         h11 => h11val,
--         h12 => h12val,
--         c2Form => L, -- maybe store list of ints?
--         cubicForm => F -- maybe store intnums? (for basis).
--         };
--     if opts#NegatedCurves =!= null then X.cache.NegatedCurves = opts#NegatedCurves;
--     if opts#GV =!= null then X.cache.GV = opts#GV; -- should be a hash table: curve classes => integers.
--     if opts#MoriCone =!= null then X.cache.MoriCone = opts#MoriCone; -- should be a hash table: curve classes => integers.
--     if opts#Heft =!= null then X.cache.Heft = opts#Heft;
--     if opts#DegreeLimit =!= null then X.cache.DegreeLimit = opts#DegreeLimit;
--     if opts#Label =!= null then X.cache.Label = opts#Label;
--     X
--     )

CY3 = new Type of HashTable
makeCY3 = method(Options => {
        GVTable => null, -- of the original, toric, hypersurface
        DegreeLimit => null, -- if null, no GVs will be computed fresh
        NegatedCurves => null,
        Label => null
        })
makeCY3(ZZ, ZZ, RingElement, RingElement) := opts -> (h11val, h12val, L, F) -> (
    X := new CY3 from {
        cache => new CacheTable,
        h11 => h11val,
        h12 => h12val,
        c2Form => L, -- maybe store list of ints?
        cubicForm => F -- maybe store intnums? (for basis).
        };
    if opts.NegatedCurves =!= null then X.cache.NegatedCurves = opts.NegatedCurves;
    if opts#GVTable =!= null then X.cache#GVTable = opts#GVTable; -- should be a hash table: curve classes => integers.
    if opts.Label =!= null then X.cache.Label = opts.Label;
    X
    )
makeCY3 CalabiYauInToric := opts -> X -> (
    makeCY3(hh^(1,1) X, hh^(1,2) X, c2Form X, cubicForm X,
        GVTable => gvTable(X, DegreeLimit => opts.DegreeLimit),
        NegatedCurves => {}
        ))

cubicForm CY3 := X -> X#cubicForm
c2Form CY3 := X -> X#c2Form
label CY3 := X -> X#Label
negatedCurves = method()
negatedCurves CY3 := X -> X#negatedCurves

-- TODO: store info in X?
-- TODO: warn that these are the "generic" complex structure mori cone and nef cone.
gvTable CY3 := opts -> X -> X.cache#GVTable
gvRays CY3 := X -> gvRays gvTable X
gvRay(CY3, List) := opts -> (X, curve) -> (gvRays X)#curve
moriCone CY3 := X -> moriCone(gvTable X, X.cache.NegatedCurves)
nefCone CY3 := X -> dualCone moriCone X
isNilpotent(CY3, List) := (X, curve) -> isNilpotent(gvTable X, curve)
    
performFlop = method()
performFlop(CY3, List) := CY3 => (X, C) -> (
    if gcd C != 1 then error "expected a primitive curve class";
    -- perform a flop
    L := c2Form X;
    F := cubicForm X;
    R := ring L;
    linform := sum for i from 0 to numgens R - 1 list C_i * R_i;
    n := first gvRay(X, C); -- is this correct?  Maybe not. -- TODO: look at the entire ray.
    makeCY3(X#h11, X#h12, L + 2*n*linform, F - n * linform^3,
        Label => splice{X, {"flop via ", C}},
        GVTable => gvTable X,
        NegatedCurves => X.cache.NegatedCurves | {C}
        )
    )

-- performFlop(CY3, List) := CY3 => (X, C) -> (
--     if gcd C != 1 then error "expected a primitive curve class";
--     -- perform a flop
--     L := X#c2Form;
--     F := X#cubicForm;
--     R := ring L;
--     linform := sum for i from 0 to numgens R - 1 list C_i * R_i;
--     n := X.cache#GV#C; -- is this correct?  Maybe not. -- TODO: look at the entire ray.

--     primcurves := unique for curve in keys X.cache#GV list (
--         curve // gcd curve
--         );
--     negatedCurves := join(X.cache#NegatedCurves, {C});
--     --negatedCurvesSet := set negatedCurves;
--     curveclasses := for c in primcurves list if c == C then -c else c;
--     mori := entries transpose rays posHull transpose matrix curveclasses; -- can be made faster.
    
--     -- make the new GV table.
--     GV1 := hashTable for kv in pairs X.cache#GV list (
--         (curve, gvnum) := kv;
--         primcurve := curve // gcd curve;
--         if primcurve == C then (-curve,gvnum) else (curve,gvnum)
--         );
--     -- WARNING: We assume that all extremal rays have non-zero gv invariant.
--     makeCY3(X#h11, X#h12, L + 2*n*linform, F - n * linform^3,
--         Label => splice{X, {"flop via ", C}},
--         MoriCone => mori,
--         NegatedCurves => join(X.cache.NegatedCurves, C),
--         -- The following are take directly from X.
--         GV => GV1,
--         Heft => X.cache#Heft,
--         DegreeLimit => X.cache#DegreeLimit
--         )
--     )

  -- gvByRay CY3 := X -> (
  --     primcurve := curve -> curve // gcd curve; -- returns {mu
  --     mult := curve -> (
  --         i := position(curve, a -> a != 0);
  --         c := primcurve curve;
  --         curve#i // c#i
  --         );
  --     --primcurves := unique for curve in keys X.cache#GV list primcurve curve;
  --     partition((c,gv) -> primcurve c, pairs X.cache#GV)
  --     )

  -- -- Use this one, and clean it up
  -- gvByRay CY3 := X -> (
  --     GVHash := X.cache#GV;
  --     degvector := X.cache#Heft;
  --     deglimit := X.cache#DegreeLimit;
  --     primcurve := curve -> (
  --         i := position(curve, a -> a != 0);
  --         c := curve // gcd curve;
  --         mult := curve#i // c#i;
  --         {c, mult}
  --         );
  --     curves := for k in keys GVHash list join(primcurve k, {GVHash#k});
  --     H := partition(val -> val#0, curves);
  --     H2 := hashTable for kv in pairs H list (
  --         val := kv#0 => kv#1/(x -> {x#1, x#2});
  --         val);
  --     -- now we change the values to be a list og gv's.
  --     toRayList := (primcurve, listOfPairs) -> (
  --         Hmult2gv = hashTable listOfPairs;
  --         d := max (listOfPairs/first);
  --         deg := dotProduct(degvector, primcurve);
  --         maxmult := max(d, deglimit//deg);
  --         for i from 1 to maxmult list if Hmult2gv#?i then Hmult2gv#i else 0
  --         );
  --     hashTable for kv in pairs H2 list (
  --         val := toRayList(kv#0, kv#1);
  --         if #val >= 2 then kv#0 => val else continue
  --         )
  --     )
  
  -- -- TODO: change to use NegatedCurves
  -- gvRay(CY3, List) := List => opts -> (X, C) -> (
  --   contentC := gcd C;
  --   if contentC =!= 1 then C = C // contentC;
  --   degvector := X.cache#Heft;
  --   GVHash := X.cache#GV;
  --   deglimit := X.cache#DegreeLimit;
  --   d := dotProduct(degvector, C);
  --   for i from 1 to max(floor(deglimit/d), 1) list (
  --       iC:= i*C;
  --       if GVHash#?iC then GVHash#iC else 0
  --       )
  --   )
  -- gvByRay CY3 := X -> (
  --     curves := keys X.cache#GV;
  --     primcurves := unique for c in curves list c // gcd c;
  --     hashTable for c in primcurves list c => gvRay(X, c)
  --     );

  heftFunction = method()
  heftFunction CalabiYauInToric := X -> (
      mori := hilbertBasisGenerators toricMoriCone(ambient X, basisIndices X);
      sum entries transpose rays dualCone posHull transpose matrix mori
      )

  -- dot = method()
  -- dot(List, List) := (v,w) -> (
  --     if #v =!= #w then error "expected vectors of the same size";
  --     sum for i from 0 to #v-1 list v#i * w#i
  --     )


end--
-- Right now, we will do it on an example with h11=3.

-- load this in dir m2-examples.
restart
--load "../Flops.m2"
load "~/utah/CYToolsM2/StringTorics/Flops.m2"
  RZ = ZZ[a,b,c];
  (Qs, Xs) = readCYDatabase("~/StringDatabases/cys-ntfe-h11-3.dbm", Ring => RZ);
  #Qs
  #Xs

-- Step 1. Find all gv invariants up to a degree bound.
  X = Xs#(53,0)
  X1 = makeCY3(X, DegreeLimit => 20)
  rays moriCone X1
  (gvRays X1)
  gvRay(X1, {1,0,0})
  -- the following needs to be improved: if every other one is zero, it is still potent
  hashTable for c in keys gvRays X1 list if (gvRay(X1, c))#-1 == 0 then c => (gvRays X1)#c else continue
  gvRay(X1, {1,-1,0}) -- type III0
  gvRay(X1, {0,1,0}) -- type I, or III4
  gvRay(X1, {-2,1,1}) -- type II?
  hashTable  for c in keys gvRays X1 list if isNilpotent(X1, c) then c => gvRay(X1, c) else continue
  select(keys gvRays X1, isNilpotent_X1)
  
  X2 = performFlop(X1, {0,1,0})
  rays moriCone X2
  
  X3 = performFlop(X2, {1,0,0})
  rays moriCone X3

-- Design: What should a GVInvariants class look like?
--  1. Has hash table as it does now.
--  2. Knows its degree limit, and grading vector.
--  3. Can compute "infinity cone": actually, should be done for 2 or 3 different degrees,
--       then compare them?
--  4. Compute ray of GV values out some distance.
--    This should use special features of the code?  Does it work on non-extremal rays?
--  5. Determine what kind of extremal ray a ray is:
--    1. nilpotent (type I)
--    2. nilpotent (type II0, type IIg)
--    3. potent ray.
--    4. is a ray in the closure of the infinity cone?  Or can we not consider this possibility?
--  6. Find non-zero-gv cone (the Mori cone in the case when the CY3 is general in moduli.
--    Handle negated curve rays.
-- For non-general CY3's it is possible for a curve to be effective, but have gv ray all 0's.
