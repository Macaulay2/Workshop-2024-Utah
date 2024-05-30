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

-- Let's write a flop function
CY3 = new Type of HashTable
makeCY3 = method(Options => {GV => null, 
        NegatedCurves => null,
        Heft => null, 
        DegreeLimit => null, 
        Label => null, 
        MoriCone => null
        })
makeCY3(ZZ, ZZ, RingElement, RingElement) := opts -> (h11val, h12val, L, F) -> (
    X := new CY3 from {
        cache => new CacheTable,
        h11 => h11val,
        h12 => h12val,
        c2Form => L, -- maybe store list of ints?
        cubicForm => F -- maybe store intnums? (for basis).
        };
    if opts#NegatedCurves =!= null then X.cache.NegatedCurves = opts#NegatedCurves;
    if opts#GV =!= null then X.cache.GV = opts#GV; -- should be a hash table: curve classes => integers.
    if opts#MoriCone =!= null then X.cache.MoriCone = opts#MoriCone; -- should be a hash table: curve classes => integers.
    if opts#Heft =!= null then X.cache.Heft = opts#Heft;
    if opts#DegreeLimit =!= null then X.cache.DegreeLimit = opts#DegreeLimit;
    if opts#Label =!= null then X.cache.Label = opts#Label;
    X
    )

-- TODO: have a function which checks that C defines a flop.

performFlop = method()
performFlop(CY3, List) := CY3 => (X, C) -> (
    if gcd C != 1 then error "expected a primitive curve class";
    -- perform a flop
    L := X#c2Form;
    F := X#cubicForm;
    R := ring L;
    linform := sum for i from 0 to numgens R - 1 list C_i * R_i;
    n := X.cache#GV#C; -- is this correct?  Maybe not. -- TODO: look at the entire ray. 
    -- make the new GV table.
    GV1 := hashTable for kv in pairs X.cache#GV list (
        (curve, gvnum) := kv;
        primcurve := curve // gcd curve;
        if primcurve == C then (-curve,gvnum) else (curve,gvnum)
        );
    mori := entries transpose rays posHull transpose matrix keys GV1; -- can be made faster.
    -- WARNING: We assume that all extremal rays have non-zero gv invariant.
    makeCY3(X#h11, X#h12, L + 2*n*linform, F - n * linform^3,
        Label => splice{X, {"flop via ", C}},
        MoriCone => mori,
        NegatedCurves => join(X.cache.NegatedCurves, C),
        -- The following are take directly from X.
        GV => X.cache.GV,
        Heft => X.cache#Heft,
        DegreeLimit => X.cache#DegreeLimit
        )
    )

  -- TODO: change to use NegatedCurves
  gvRay(CY3, List) := List => opts -> (X, C) -> (
    contentC := gcd C;
    if contentC =!= 1 then C = C // contentC;
    degvector := X.cache#Heft;
    GVHash := X.cache#GV;
    deglimit := X.cache#DegreeLimit;
    d := dotProduct(degvector, C);
    for i from 1 to floor(deglimit/d) list (
        iC:= i*C;
        if GVHash#?iC then GVHash#iC else 0
        )
    )

  heftFunction = method()
  heftFunction CalabiYauInToric := X -> (
      mori := hilbertBasisGenerators toricMoriCone(ambient X, basisIndices X);
      sum entries transpose rays dualCone posHull transpose matrix mori
      )

  dot = method()
  dot(List, List) := (v,w) -> (
      if #v =!= #w then error "expected vectors of the same size";
      sum for i from 0 to #v-1 list v#i * w#i
      )

end--
-- Right now, we will do it on an example with h11=3.

-- load this in dir m2-examples.
restart
load "../Flops.m2"
  -- debug for e.g. hilbertBasisGenerators.
  RZ = ZZ[a,b,c];
  (Qs, Xs) = readCYDatabase("../Databases/cys-ntfe-h11-3.dbm", Ring => RZ);
  #Qs
  #Xs

-- Step 1. Find all gv invariants up to a degree bound.
  X1 = Xs#(53,0)
  gvX1 = hashTable for kv in pairs gvInvariants(X1, DegreeLimit => 20) list (toList kv#0, kv#1)
  X1 = makeCY3(hh^(1,1) X1, hh^(1,2) X1, c2Form X1, cubicForm X1,
      GV => gvX1,
      MoriCone => toricMoriConeCap X1,
      DegreeLimit => 20,
      Heft => heft X1)

  X1.cache.GV#{0,1,0}
  X1#h11
  X1#h12
  X1#c2Form
  gvRay(X1, {0,1,0})  
  
  X2 = performFlop(X1, {0,1,0})
  for c in X2.cache.MoriCone list c => gvRay(X2, c)

  
  X3 = performFlop(X2, {1,0,0})
  X3.cache.MoriCone
  GV = gvInvariants(X, DegreeLimit => 20)
  GVcone = (keys GV)/toList//matrix//transpose//posHull
  rays GVcone
  -- Question: what is the degree limit method?
  

  deglimit = 20
  hf = heftFunction X
  assert(hf == {1,1,3})
  (keys GV)/(v -> dot(toList v, hf))
  allGV = (keys GV)/toList//set
  select(keys allGV, v -> dot(v, hf) <= deglimit // 2)
  allrays = unique for v in keys allGV list (
      c := gcd v;
      for v1 in v list v1//c
      )
  goodray = x -> (
      d := dot(x, hf);
      topval := floor(deglimit/d);
      all(splice{1..topval}, i -> member(i * x, allGV))
      )
  H = partition(goodray, allrays);
  rays posHull transpose matrix H#true
  matrix{hf} * oo
  matrix{hf} * transpose matrix H#false
  for x in oo list (
      d := dot(x, hf);
      topval := floor(deglimit/d);
      if not all(splice{1..topval}, i -> member(i * x, allGV))
        then continue
        else x
      )  
  rays posHull transpose matrix oo

  for k in keys allGV list (
      -- we will select all of the ones
      )

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
