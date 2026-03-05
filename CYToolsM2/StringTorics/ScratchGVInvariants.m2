-- Scratch code from GVInvariants.m2
-- Tests using databases, too slow or requiring database files for regular tests

///
  restart
  debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
  DB3 = "../Databases/cys-ntfe-h11-3.dbm"
  DB3 = "./StringTorics/Databases/cys-ntfe-h11-3.dbm"
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ);
  (Qs, Xs) = readCYDatabase(DB3, Ring => RZ);
  sort keys Xs
  X = Xs#(1,0)
  deglimit = 10
  gvX = gvInvariants(X, DegreeLimit => deglimit)
  gvTable X

  mori = toricMoriConeCap X
  extremalRayGVs X
  classifyExtremalCurves X

  elapsedTime gvConeX = gvCone X
  partitionGVConeByGV(X, 5)
  netList oo


///

///
  restart
  debug needsPackage "StringTorics"
  DB4 = "../Databases/cys-ntfe-h11-4.dbm"
  RZ = ZZ[a,b,c,d]
  (Qs, Xs) = readCYDatabase(DB4, Ring => RZ);
  sort keys Xs

  X = Xs#(10, 0)
  ambient X
  deglimit = 10
  mori = toricMoriConeCap X
  gvX = gvInvariants(X, DegreeLimit => deglimit)
  degvec = heft X;
  gvX = gvInvariants(X, DegreeLimit => deglimit)
  partition(c -> classifyExtremalCurve(gvX, c, deglimit, degvec), mori)
///

"TEST"
///
-- Good test, TODO: place this back in once GVinvariants are working again.
-*
  restart
  debug needsPackage "StringTorics"
*-
  debug StringTorics
  DB3 = "../Databases/cys-ntfe-h11-3.dbm"
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ);
  (Qs, Xs) = readCYDatabase(DB3, Ring => RZ);

  gvOK = for lab in sort keys Xs list (
    << lab << endl;
    try gvInvariants(Xs#lab, DegreeLimit => 10) then lab else continue
    )

  X = Xs#(0,0)
  gvTable X
  gvRays X
  extremalRayGVs X
  gvRays X
  classifyExtremalCurves X
  gvCone X

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
  DB5 = "../Databases/cys-ntfe-h11-5.dbm"
  RZ = ZZ[a,b,c,d,e]
  (Qs, Xs) = readCYDatabase(DB5, Ring => RZ);

  allXs = sort keys Xs;
  gvOK = for lab in allXs list (
    << lab << endl;
    try gvInvariants(Xs#lab, DegreeLimit => 10) then lab else continue
    )
  gvBAD = toList(set allXs - set gvOK)
  for lab in gvBAD list (hh^(1,2) Xs#lab, intersectionNumbers Xs#lab)
  for lab in gvBAD list (hh^(1,2) Xs#lab, #intersectionNumbers Xs#lab)
  #gvOK
  for lab in gvOK list (
    << "doing " << lab << endl;
    elapsedTime extremalRayGVs Xs#lab
    )
  netList oo

  for lab in gvOK list (
    << "doing " << lab << endl;
    elapsedTime partitionGVConeByGV(Xs#lab, 5)
    )

  X1 = Xs#(24,0)
  X2 = Xs#(24,2)
  X3 = Xs#(24,3)
  gv1 = gvInvariants(X1, DegreeLimit => 10)
  gv2 = gvInvariants(X2, DegreeLimit => 10)
  gv3 = gvInvariants(X3, DegreeLimit => 10)
  partitionGVConeByGV(X1, 5)
  partitionGVConeByGV(X2, 5)
  partitionGVConeByGV(X3, 5)
  cubicForm X1
  cubicForm X2
  cubicForm X3
  c2Form X1
  c2Form X2
  c2Form X3

  -- let's check the ones that should be equivalent.
  H1 = hessian cubicForm X1
  H2 = hessian cubicForm X2
  H3 = hessian cubicForm X3
  -- the first step: X1, X3 have the same singular stuff.  X2 has different.
  singularPointMatches(cubicForm X1, cubicForm X2) -- none exist, so not equivalent
  singularPointMatches(cubicForm X1, cubicForm X3)
  --L2 = hessianMatches(cubicForm X1, cubicForm X2) -- slow, 120 sec, gives empty list
  L3 = hessianMatches(cubicForm X1, cubicForm X3)

  RQ = QQ[a,b,c,d,e]
  md = matchingData({singularPointMatches, hessianMatches}, cubicForm X1, cubicForm X3, RQ)
  md1 = matchingData({singularPointMatches, hessianMatches}, cubicForm X1, cubicForm X1, RQ)
  tryEquivalences(md, RQ, (c2Form X1, cubicForm X1), (c2Form X3, cubicForm X3))
  -- for a in L3 list (
  --     isEquivalent((c2Form X1, cubicForm X1), (c2Form X3, cubicForm X3), a)
  --     )
  -- oo/first
  -- result is: there exists a matrix relating them!  What is it?
  -- position(oo, x -> x)
  -- L3#20

  allXs = sort keys Xs;
  for lab in allXs list (
      << "doing " << lab << endl;
      --if ans === null then (gvInvariants(Xs#lab, DegreeLimit => 10); "OK") else continue
      a := try gvInvariants(Xs#lab, DegreeLimit => 10) then "OK" else "BAD";
      --if a == "BAD" then (print (hh^(1,2) Xs#lab, #intersectionNumbers Xs#lab));
      lab => a
      );
  tally (oo/last)

  for lab in gvOK list (
      X = Xs#lab;
      gvX = gvInvariants(X, DegreeLimit => 10);
      deglimit = 10;
      degvec = heft X;
      mori = toricMoriConeCap X;
      a = partition(c -> classifyExtremalCurve(gvX, c, deglimit, degvec), mori);
      --if ans#0 === INCONSISTENT then {{first a}, {last a}} else if ans#0 === CONSISTENT then {first a, last a, ans#1} else {a, ans}
      print (a => ans);
      ans
      )

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
      print (a => ans);
      ans
      )

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
      print (a => ans);
      ans
      )

///
