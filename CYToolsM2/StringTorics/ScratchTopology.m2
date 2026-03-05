-- Scratch code from Topology.m2
-- Analysis of h11=3 and h11=4 examples using topology classification

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

  smallset = select(sort keys Xs, lab -> hh^(1,2) Xs#lab == 148) --
  T = topologySet(sort smallset, Xs)
  allTsame = combineIfSame(T, X -> (c2Form X, cubicForm X))
  info allTsame
  elapsedTime T1 = separateIfDifferent(allTsame, X -> invariantsAll X)
  info T1
  netList T1#"Sets"
  T2 = separateByGV T1

///
