-- Scratch code from CYTools.m2
-- Database exploration, GV invariants, and topological classification

///
-- Let's try this
-*
  restart
  debug needsPackage "StringTorics"
  needs "CYTools.m2"
*-

  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ);

  -- readCYDatabase needs updating
  DB3 = databaseLOC | "/cy3-h11-3.dbm"
  (Qs, Xs) = readCYDatabase(DB3, Ring => RZ);

  sort keys Xs

  for lab in sort keys Xs list (
    << lab << endl;
    try gvInvariants(Xs#lab, DegreeLimit => 10) then lab else continue
    )

  for lab in sort keys Xs list (
    hh^(1,2) Xs#lab => intersectionNumbers Xs#lab
    )

  gvOK = for lab in sort keys Xs list (
    << lab << endl;
    try gvInvariants(Xs#lab, DegreeLimit => 10) then lab else continue
    )

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

  RQ = QQ[a,b,c]
  md = matchingData({singularPointMatches, hessianMatches}, cubicForm X1, cubicForm X3, RQ)
  md1 = matchingData({singularPointMatches, hessianMatches}, cubicForm X1, cubicForm X1, RQ)
  tryEquivalences(md, RQ, (c2Form X1, cubicForm X1), (c2Form X3, cubicForm X3))
  -- oo is a list.  False means: not equivalent under this set of matching data.
  -- result is: there exists a matrix relating them!  What is it?
  -- for a in L3 list isEquivalent((c2Form X1, cubicForm X1), (c2Form X3, cubicForm X3), a)
  -- oo/first
  -- position(oo, x -> x)
  -- L3#20

allXs#0 -- return hashtable: gv value => list of rays.

allXs = sort keys Xs;
for lab in gvOK list
  partitionGVConeByGV(Xs#lab, 5)

///
