TEST ///
  -- XX TODO: being worked on now 29 June 2023
  -- Checking the methods for CYPolytope's
  -- We eventually want to test this for: favorable, non-favorable, torsion V.

-*
  restart
  needsPackage "StringTorics"
*-
  debug needsPackage "StringTorics"
  
  tope = KSEntry "4 13  M:34 13 N:12 10 H:7,29 [-44] id:40
   1   0   0   0   0  -2  -2   1   2   1   2   2  -2
   0   1   0   0   0   2   1  -1  -2  -2   0  -2   0
   0   0   1  -1   0  -1   0  -1   1   0   1  -1   1
   0   0   0   0   1   1   2   1  -2  -1  -2   0   0
   "
  Q = cyPolytope(tope, ID => 40)  
  assert isFavorable Q
  basisIndices Q
  Q.cache#"toric basis indices"
  Q.cache#"basis indices"
  transpose matrix degrees Q
  X = makeCY Q
  toricIntersectionNumbers X
  intersectionNumbers X
  assert(c2 X === {-4, 10, 12, 10, 18, 0, 10})
  cubicForm X
  c2Form X
  assert(ring cubicForm X === ring c2Form X)

  nonfavTope = KSEntry "4 7  M:28 7 N:11 7 H:7,27 [-40] id:9
   1   0   3   0  -1  -3  -6
   0   1   2   0   0  -2  -5
   0   0   4   1   0  -2  -8
   0   0   0   2   2   2  -4
   "

  Q2 = cyPolytope nonfavTope   
  assert not isFavorable Q2
  assert(hh^(1,1) Q2 == 7)
  assert(hh^(1,2) Q2 == 27)
  -- the following will need to change...
  assert(basisIndices Q2 === {0, 1, 2, 3, 4, (9, 0), (9, 1)}) -- this could change if the algorithm changes.
  Q2.cache#"toric basis indices" === {0, 1, 2, 3, 4, 9}
  findTwoFaceInteriorDivisors Q2
  netList annotatedFaces Q2
  assert(
      rays Q2 === {{-1, -1, 1, -1}, 
      {-1, -1, 1, 1}, 
      {-1, -1, 2, -1}, 
      {-1, 0, 1, -1}, 
      {-1, 3, -1, 0}, 
      {1, 0, -1, 0}, 
      {3, -1, -2, 1}, 
      {-1, -1, 1, 0}, 
      {1, -1, 0, 0}, 
      {-1, 1, 0, 0}}
      )
  Xs = findAllCYs Q2; -- what if I only want the NTFE ones?  FIX.
  #Xs
  (netList toricIntersectionNumbers Xs#0, netList intersectionNumbers Xs#0)
  c2 Xs#0 -- fix me
  intersectionNumbers Xs#0 -- fix me
  ring cubicForm Xs#0 === ring c2Form Xs#0
  c2Form Xs#0
  
  vertices polytope(Q2, "M")
  vertices polytope(Q2, "N") -- notice these are NOT in the order of the rays of Q2!
  transpose matrix degrees Q2
  
  transpose matrix rays Q2
  
  cubicForm Xs#0
///

///
  -- Checking on the interface of the package.
-*
  restart
  needsPackage "StringTorics"
*-
  debug needsPackage "StringTorics"
  topes = kreuzerSkarke(5, Limit => 10000);
  #topes == 4990
  A = matrix topes_40 -- this will be vertices of a polytope in the M lattice 
  -- We need to get to a triangulation of the dual polytope...
  X0 = cyPolytope topes_40
  X = makeCY(X0, Ring => (RZ = ZZ[s_1..s_5]))

  -- X = calabiYau(A, Lattice => "M") -- A must define a reflexive polytope.

  V = ambient X
  aX = abstractVariety(X, base(a,b,c,d,e))
  intersectionRing aX -- defines integral.
  intersectionRing V -- defines integral.
  topX = topologicalData X
  cubicForm topX
  isFavorable cyPolytope X
  
  intersectionNumbers X
  toricIntersectionNumbers X
  c2 X
  cubicForm X
  c2Form X
///

TEST ///
-*
  restart
  needsPackage "StringTorics"
*-
  str4 =         "         1   0   0   0   1   1   0  -1  -1  -2  -4
                  0   1   1   0  -2   2   3  -1  -4   1  -1
                  0   0   2   0  -2   4   4  -1  -4  -2  -4
                  0   0   0   1   0  -2  -2   2   2   0   2  
                  "
  A = matrixFromString str4
  P = convexHull A
  V = reflexiveToSimplicialToricVariety P
  X = completeIntersection(V, {-toricDivisor V})
  pt = base(a,b,c)
  Xa = abstractVariety(X, pt)
  I = intersectionRing Xa
  a*t_1
  hodgeDiamond X
  assert(h11OfCY P == 6)
  assert(h21OfCY P == 50)
  a*t_1 -- should still work...
///


TEST ///
-*
  restart
*-
  debug needsPackage "StringTorics"
  -- augmentWithOrigin
    mat = "   1   0   0   1  -3   3   3  -3   5
            0   1   0   0   2  -4  -6   4  -8  
            0   0   1   0   2  -2  -2   0  -4 
            0   0   0   2   0  -4  -6   6  -6  "
  A = matrixFromString mat;
  assert(
      augmentWithOrigin A 
      == 
      matrix {
          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
          {1, 0, 0, 1, -3, 3, 3, -3, 5, 0}, 
          {0, 1, 0, 0, 2, -4, -6, 4, -8, 0}, 
          {0, 0, 1, 0, 2, -2, -2, 0, -4, 0},
          {0, 0, 0, 2, 0, -4, -6, 6, -6, 0}
          }
      )
///


TEST ///
  -- readSageTriangulation
  Ts = readSageTriangulations sageTri
  
  -- let's also switch to a different choice of rays
  -- Bstr is the lattice points (a.k.a. rays) used by Cody's code in Sage.
    Bstr = "[ 1 -1 -1  1 -1 -1 -1  1  1  0]
[ 0  1  0 -1  1  0  0 -1 -1  0]
[-1  0  0  0  1  1  0  0  0  0]
[ 2  0  0  1 -1 -1 -1 -1  0  0]"
  -- rays coming from M2
  Amat = matrix {{-1, -1, -1, -1, -1, 1, 1, 1, 1}, {0, 0, 0, 1, 1, -1, -1, -1, 0}, {0, 0, 1, 0, 1, 0, 0, 0, -1}, {-1, 0, -1, 0, -1, -1, 0, 1, 2}}
  Bmat = matrixFromString Bstr
  
  (fromM2, toM2) = matchNonZero(Amat, Bmat)
  applyPermutation(toM2, Ts)
  
  assert(
      applyPermutation(toM2, {0,1,2,3}) 
      == 
      {1,3,7,8}
      )
///

-- commenting this out.  It should really be a test in ReflexivePolytopesDB
///
  L1 = kreuzerSkarke(5,57, Access => "wget");
  assert(#L1 == 197)

  L2 = kreuzerSkarke(5,57, Access => "curl");
  assert(#L2 == 197)
  
  assert(L1 === L2)
///

TEST ///
  -- XXX  TODO: This test should be in Triangulations?
  -- Test functionality of triangulations, part 1. Basic tests
  -- 
-*  
restart
needsPackage "StringTorics"
*-
  -- WARNING: currently, we use the polar dual of a convex polytope to determine
  -- minimal faces, etc.  However, for this to work, the convex polytope
  -- MUST contain the origin in the interior.
  -- Here are the functions that won't work if this is not the case:
  --polar P -- doesn't work as expected, since the origin is not an interior point of the polytope.
  --  minimalFace(P, pts)
  --  dualFace(P,f)
  --  latticePointList(P,f)
  --  genus(P,f)
  --  annotatedFaces P
  --  annotatedFaces(i,P)
  -- TODO: it would be possible to get all of these to work, except genus, and dualFace,
  --  by implementing them a bit differently.  Is it worth it? Most polytopes of interest will contain
  -- the origin in the interior?

  A = transpose matrix{{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}}
  A1 = matrix{{8:1}} || A
  P = convexHull A
  assert(dim P == 3)

  vertices P
  latticePoints P
  faces P
  faces(1,P)

  LP = latticePointList P
  assert(set LP === set {{1, 1, 1}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 0}})
  assert(set vertexList P === set {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}})
  assert(entries transpose vertexMatrix P == vertexList P)

  assert(# faceList P == 27)
  faceDimensionHash P
  for f in faceList P do assert (2^(dim(P,f)) === #f)
  assert(faceList(0,P) == for i from 0 to 7 list {i})

  Amat = transpose matrix LP
  tri = regularFineTriangulation Amat
  assert naiveIsTriangulation tri
  assert topcomIsTriangulation(Amat, max tri)
  
  -- check what happens if Amat is homoogenized:
  AmatH = Amat || matrix{{8:1}}
  naiveIsTriangulation(AmatH, max tri) -- this should be false...? TODO: this is a bug!!
  assert not topcomIsTriangulation(AmatH, max tri) -- good! it complains that the index sets are not full dimensional (I think that is good?)
  
  assert(affineCircuits tri == affineCircuits(Amat, max tri))
  for x in affineCircuits tri list bistellarFlip(tri, x)
  neighbors tri
  
  generateTriangulations tri
  allTriangulations(Amat, RegularOnly => false, ConnectedToRegular => false, Fine => false)
  
  generateTriangulations(tri, RegularOnly => true)
  assert(# allTriangulations Amat == 74)
  assert(# generateTriangulations(tri, RegularOnly => true) == 74)

  -- let's check that this is a triangulation.
  -- part of what we are checking: calls relative to homogenization are correct, and types make sense.
  -- part 1: for each oriented circuit
  circs = orientedCircuits transpose matrix LP
  assert(circs == {
          {{0, 4}, {1, 2}}, {{0, 4, 5}, {1, 6}}, {{0, 4, 6}, {2, 5}}, {{0, 5}, {1, 3}}, 
          {{0, 5, 6}, {3, 4}}, {{0, 6}, {2, 3}}, {{0, 7}, {1, 2, 3}}, {{0, 7}, {1, 6}}, 
          {{0, 7}, {2, 5}}, {{0, 7}, {3, 4}}, {{0, 7}, {4, 5, 6}}, {{1, 2, 7}, {3, 4}}, 
          {{1, 3, 7}, {2, 5}}, {{1, 6}, {2, 3, 7}}, {{1, 6}, {2, 5}}, {{1, 6}, {3, 4}}, 
          {{1, 7}, {4, 5}}, {{2, 5}, {3, 4}}, {{2, 7}, {4, 6}}, {{3, 7}, {5, 6}}
          }
      )
  Ts = generateTriangulations(tri, RegularOnly => true)
  assert(Ts/isWellDefined//unique == {true})

  -- it is possible that another triangulation would be output.
  -- the following is more like possible code to decide if a subset is a triangulation.
  -- It works currently, so I'll keep it...
  assert(max tri == {{0, 1, 2, 3}, {1, 2, 3, 4}, {1, 3, 4, 5}, {2, 3, 4, 6}, {3, 4, 5, 6}, {4, 5, 6, 7}})
  for c in circs list (
      n1 := # select(max tri, t -> isSubset(c#0, t));
      n2 := # select(max tri, t -> isSubset(c#1, t));
      ok := n1 == 0 or n2 == 0 or n1 == #(max tri) or n2 == #(max tri);
      (n1,n2,ok)
      )

  walls = tri//max/(x -> subsets(x, #x-1))//flatten
  nfacets = tally walls
  facs = (faces(1,P))/first
  walls = partition(k -> nfacets#k, keys nfacets)
  facs = for f in facs list latticePointList(P, f)
  for w in walls#1 list (
      # select(facs, f -> isSubset(w, f))
      )
  for w in walls#2 list (
      # select(facs, f -> isSubset(w, f))
      )
  
  walls#2 -- 6 walls here.  Compute the vector for each.
  matrix {for w in walls#2 list (
      --w = {2,3,4} -- a wall
      circ := select(max tri, t -> isSubset(w, t));
      others := circ/(c -> toList(set c - set w));
      elems := (flatten others) | w;
      print elems;
      id_(ZZ^8)_elems * syz AmatH_elems
      )}
      --id_(ZZ^8)_{1,6,2,3,4} * syz AmatH_{1,6,2,3,4} 
///

TEST /// 
  -- XXX  
  -- simple test of the polyhedral functions here: on the cube with the origin as its only 
  -- interior point
  needsPackage "StringTorics"
  P = hypercube 3
  assert isCompact P -- These functions fail if P is not a polytope.
  assert(
    vertices P  == matrix(QQ, {
        {-1, 1, -1, 1, -1, 1, -1, 1}, 
        {-1, -1, 1, 1, -1, -1, 1, 1}, 
        {-1, -1, -1, -1, 1, 1, 1, 1}
        })
    )
  vertices2 = vertexList P
  vertices3 = entries transpose vertexMatrix P
  assert(vertices2 == vertices3)
  
  LP = latticePointList P
  LP2 = transpose entries lift(matrix {latticePoints P}, ZZ)
  assert(set LP === set LP2)
  assert(27 == # LP)
  assert(27 == # LP2)

  assert(# faceList P == 27)
  faceDimensionHash P
  for f in faceList P do assert (2^(dim(P,f)) === #f)
  assert(faceList(0,P) == for i from 0 to 7 list {i})

  hashTable for pt in LP list pt => minimalFace(P, pt)      
  for pt in LP do (
      f := minimalFace(P, pt);
      assert(2^(# select(pt, v -> v == 0)) == #f)
      )
  P2 = polar P
  for f in faceList P list (f, dualFace(P, f))
  faces P
  faceList P
  assert(set((flatten values faces P)/first) === set faceList P)
  for f in faceList P list (
      latticePointList(P, f)
      )
  
  -- annotatedFaces
  netList annotatedFaces P
  for p in annotatedFaces P do (
      f := p#1;
      assert(p#0 == dim(P,f));
      lps := latticePointList(P, f);
      assert(p#2 == lps);
      assert(p#3 == # interiorLatticePointList(P,f));
      assert(p#4 == genus(P,f))
      )  

-*
  TODO: BUG!! These give segfaults on my Apple M1.  
  debugLevel = 3  
  regularFineTriangulation vertices P
  regularStarTriangulation P -- it is a shame that this returns something different from regularFineTriangulation.
*-  
///

TEST ///
  -- Test of triangulation code
-*
  restart
  needsPackage "StringTorics"
*-  
  -- XXX
  topes = kreuzerSkarke 3;
  A = matrix topes_30
  Q = cyPolytope topes_30
  hh^(1,1) Q
  P1 = polytope(Q, "M")
  P2 = polytope(Q, "N")
  vertices P1
  vertices P2
  P2 == polytope Q
  findAllFRSTs P2
  Amat = transpose matrix latticePointList polytope Q
  assert(Amat == transpose matrix {{-1, -1, 0, 0}, {-1, -1, 0, 1}, {-1, -1, 2, 0}, {-1, 0, 0, 0}, {1, -1, -1, 1}, {1, 2, -1, -1}, {-1, -1, 1, 0}, {0, 0, 0, 0}})
///

TEST ///
  -- Test of triangulation code
-*
  restart
  needsPackage "Triangulations"
*-  
  Amat = transpose matrix {{-1, -1, 0, 0}, {-1, -1, 0, 1}, {-1, -1, 2, 0}, {-1, 0, 0, 0}, {1, -1, -1, 1}, {1, 2, -1, -1}, {-1, -1, 1, 0}, {0, 0, 0, 0}}
  regularSubdivision(Amat, matrix{{0,0,1,3,6,9,20,30}}) -- seems incorrect.
  TRI = regularFineTriangulation Amat -- is this including the origin automatically?
  wts = regularTriangulationWeights TRI
  -- check that this is a triangulation!
  TRI2 = regularSubdivision(Amat, matrix{{2, 4, 2, 0, 0, 0, 0, 0}}) -- good!
  TRI = TRI//max/sort//sort
  TRI2 = TRI2/sort//sort
  assert(TRI === TRI2) -- works!
  assert topcomIsTriangulation(Amat, TRI)
  assert naiveIsTriangulation(Amat, TRI)

  -- let's check 'affineCircuits'
  C = affineCircuits(Amat, TRI)  
  4! * volumeVector(Amat, TRI)
  bistellarFlip(TRI, C_0) === null
  bistellarFlip(TRI, C_1) === null

  TRI2 = bistellarFlip(TRI, C_2)
  wts2 = regularTriangulationWeights(Amat, TRI2)
  assert(wts2 == {1, 1, 2, 0, 0, 0, 0, 0}) -- doesn't really need to be the same.
  TRI2' = regularSubdivision(Amat, matrix {wts2}) -- good!
  assert(TRI2 == TRI2')

  TRI3 = bistellarFlip(TRI, C_3)
  wts3 = regularTriangulationWeights(Amat, TRI3)
  TRI3' = regularSubdivision(Amat, matrix {wts3}) -- good!
  assert(TRI3 == TRI3')

  bistellarFlip(TRI, C_4) === null
  
  bistellarFlip(TRI, C_5) === null

  TRI6 = bistellarFlip(TRI, C_6)
  wts6 = regularTriangulationWeights(Amat, TRI6)
  assert(wts6 == {-1, 3, 2, 0, 0, 0, 0, 0}) -- doesn't need to be the case.
  TRI6' = regularSubdivision(Amat, matrix{wts6}) -- good!
  assert(TRI6' == TRI6)

  4! * volumeVector(Amat, TRI)
  4! * volumeVector(Amat, bistellarFlip(TRI, C_2))
  4! * volumeVector(Amat, bistellarFlip(TRI, C_3))
  4! * volumeVector(Amat, bistellarFlip(TRI, C_6))
///

TEST ///
  -- Test functionality of triangulations, part 2. reflexive polytope in 4D
  -- We try one with h11=3
  needsPackage "StringTorics"
  tope = "4 13  M:64 13 N:8 7 H:3,51 [-96]
   1   0   0   2  -2   0  -1  -1   0   2   2  -3  -3
   0   1   1   3  -5   2   0   0   2   4   4  -6  -6
   0   0   4   0  -4   5   4  -1   0   1   0  -4  -5
   0   0   0   4  -4   1  -1  -1   1   5   5  -5  -5"
  A = matrix first kreuzerSkarke tope
  P = convexHull A
  P2 = polar P
  V = reflexiveToSimplicialToricVariety P
  assert isFavorable P
  assert isCompact P
  assert isCompact P2
  -- regularFineTriangulation (from topcom)
  -- regularFineStarTriangulation
  -- isRegularTriangulation (from topcom)
  nonzeroLP = transpose matrix drop(latticePointList P2,-1)
  LP = transpose matrix latticePointList P2
  chirotope nonzeroLP

  regularFineTriangulation nonzeroLP
  regularFineTriangulation LP
  regularStarTriangulation P2
///

TEST ///
  -- this is an example of a polytope with 200 lattice points.
  tope = "4 12  M:72 12 N:276 12 H:200,50 [300]
          1    0    0    0  -14   -1  -17   -4   -5   -9  -15  -29
          0    1    0    0   -9   -1  -11   -2   -4   -6  -10  -20
          0    0    1    0   -3    0   -4   -2   -2   -4   -6  -10
          0    0    0    1   -1    1   -1    1    2    2    2    2"
  A = matrix first kreuzerSkarke tope
  P = convexHull A
  P2 = polar P
  elapsedTime faces P2;
  elapsedTime halfspaces P2
  LP = latticePoints P2
  # faces(1,P2)
  elapsedTime regularFineTriangulation matrix transpose latticePointList P2;
  V = reflexiveToSimplicialToricVariety P
  rays V
  max V
///

TEST ///
  -- testing triangulations of fans
  -- Our plan: start with example #26 from Kreuzer-Skarke with h11=5, h12=57
  --  this one has a number of triangulations.
  -- How do we create Amat?  This is the way:
  --      4 11  M:58 11 N:10 8 H:5,51 [-92]
  -- XXXXXXXXXX This test is failing May 2022.
-*
  restart
*-
  needsPackage "StringTorics"
  mat = "  1   1   1  -1   0   1   1  -1  -3  -1  -3
         0   2   0   0   0   0   2  -2  -2  -2  -4
         0   0   2  -2   0  -2   2   2  -2   4   2
         0   0   0   0   1  -2   0   2   0   2   2"
  A = matrixFromString mat
  P1 = convexHull A -- (extremal) vertices in RR^4
  P2 = polar P1
  LP = latticePoints P2 -- these will be the origin, the extremal vertices of P2 and possibly some more.
  Amat = matrix {select(LP, x -> x != 0)}
  elapsedTime   allTRIS = generateTriangulations(Amat, RegularOnly => true); -- removed in commit 33e77a592d2890c7ebf134e13b95d5915a624039

  -- now let's change these to sage indexing
    Bstr = "[ 1 -1 -1  1 -1 -1 -1  1  1  0]
            [ 0  1  0 -1  1  0  0 -1 -1  0]
            [-1  0  0  0  1  1  0  0  0  0]
            [ 2  0  0  1 -1 -1 -1 -1  0  0]"
  -- rays coming from M2
  Bmat = matrixFromString Bstr
  
  (fromM2, toM2) = matchNonZero(Amat, Bmat)

  Ts = readSageTriangulations sageTri
  -- checkFan is no longer available.
  --elapsedTime for T in Ts do time checkFan(Bmat, T) -- this takes a while (24 seconds), too long for testing
  --elapsedTime checkFan(Bmat, Ts_5)
  applyPermutation(fromM2, allTRIS/max)

  -- the following are all in this list
  time TRI = regularFineStarTriangulation Amat -- this appears to not be returning regular triangulations?
  Amat0 = Amat | transpose matrix{{0,0,0,0}};
  wts = regularTriangulationWeights(Amat0, TRI)
  regularSubdivision(Amat0, matrix{wts})
  assert(oo == TRI)

  -- make a toric variety from one of the triangulations: (fine regular star...)
  -- elapsedTime assert({(true, true, true)} === 
  --     unique for T in allTRIS list (
  --     X = normalToricVariety(entries transpose Amat, max T);
  --     time (isWellDefined X, isSimplicial X, isComplete X, isProjective X)
  --     )
  -- )
///

TEST ///
  -- analyze polytopes with h11=30, h12=50.
  -- grabbed 3000 examples, but there are more!
-*
  restart
  needsPackage "StringTorics"
*-
  -- one of the "favorable" examples
  --  4 7  M:51 7 N:45 7 H:30,50 [-40]
  mat = "    1    0    0    0    0   -4  -10
    0    1    0    0    2    2  -10
    0    0    1    0   -3   -7    5
    0    0    0    1    2    2   -4"

  A = matrixFromString mat
  P1 = convexHull A -- (extremal) vertices in RR^4
  P2 = polar P1
  LP = latticePoints P2 -- these will be the origin, the extremal vertices of P2 and possibly some more.
  Amat = matrix {select(LP, x -> x != 0)}

  elapsedTime regularStarTriangulation P2;
  elapsedTime   TRI = regularFineStarTriangulation Amat;
  assert(sort unique flatten TRI == toList(0..numcols Amat))
///

TEST ///
      -- 4 10  M:105 10 N:71 10 H:50,80 [-60]
     mat = "     1    0    0    0   -1   -1   -1   -5  -15  -15
          0    1    0    0    0   -2   -2   -6   -8  -16
          0    0    1    0    1    1    0   -2   -6   -6
          0    0    0    1   -1    1    2    4    0    8"

  A = matrixFromString mat
  P1 = convexHull A -- (extremal) vertices in RR^4
  P2 = polar P1
  LP = latticePoints P2 -- these will be the origin, the extremal vertices of P2 and possibly some more.
  Amat = matrix {select(LP, x -> x != 0)}

  elapsedTime   TRI = regularFineStarTriangulation Amat;
  assert(sort unique flatten TRI == toList(0..numcols Amat))
///

TEST ///
  -- creating simplicial toric varieties from data base
      -- 4 10  M:105 10 N:71 10 H:50,80 [-60]
     mat = "     1    0    0    0   -1   -1   -1   -5  -15  -15
          0    1    0    0    0   -2   -2   -6   -8  -16
          0    0    1    0    1    1    0   -2   -6   -6
          0    0    0    1   -1    1    2    4    0    8"

  A = matrixFromString mat
  P1 = convexHull A
  P2 = polar P1
  (LP,tri) = regularStarTriangulation P2
  V = normalToricVariety(LP,tri)  
  assert isSimplicial V
  assert not isSmooth V
  --assert elapsedTime isWellDefined V -- this currently takes some time

  (LP,tri) = regularStarTriangulation(2,P2)
  V = normalToricVariety(LP,tri)  
  assert isSimplicial V
  assert not isSmooth V
  --assert isWellDefined V -- this currently takes some time
  
  V1 = reflexiveToSimplicialToricVariety P1
  assert isSimplicial V
  assert not isSmooth V
///

TEST ///
  -- toric complete intersection cohomology code
  -- YYY
-*
  restart
  needsPackage "StringTorics"
*-
  mat = "  1   1   1  -1   0   1   1  -1  -3  -1  -3
         0   2   0   0   0   0   2  -2  -2  -2  -4
         0   0   2  -2   0  -2   2   2  -2   4   2
         0   0   0   0   1  -2   0   2   0   2   2"
  A = matrixFromString mat
  P = convexHull A
  V = reflexiveToSimplicialToricVariety P
  KV = toricDivisor V
  X = completeIntersection(V, {-KV})
  hodgeDiamond X
  assert(h11OfCY P == 5)
  assert(h21OfCY P == 51)
  assert(cohomologyVector X == {1,0,0,1})
  cohomologyVector(X, degree(V_0+V_1+V_2))
  allsums = drop(subsets for i from 0 to #rays V - 1 list V_i, 1);
  allsums = allsums/sum;
  elapsedTime for i from 0 to 100 list cohomologyVector(X, degree allsums_i)
  for i from 0 to 10 list cohomologyVector(X, 2 * degree allsums_i)
///

-- example: regular star triangulations
///
-*
  restart
*-
  needsPackage "StringTorics"
  mat = "  1   1   1  -1   0   1   1  -1  -3  -1  -3
         0   2   0   0   0   0   2  -2  -2  -2  -4
         0   0   2  -2   0  -2   2   2  -2   4   2
         0   0   0   0   1  -2   0   2   0   2   2"
  A = matrixFromString mat
  P = convexHull A
  V = reflexiveToSimplicialToricVariety P
  pts = transpose matrix rays V
  pts = pts || matrix{{numColumns pts : 1}}
  T = max V -- triangulation
  -- is this list correct, or do we need to "homogenize 'pts'?
  annotatedFaces polar P
  ac = select(affineCircuits(pts,T), x -> #x#0 > 1 and #x#1 > 1)
  ac = unique(ac/sort//sort)
  netList oo
  volumeVector(pts, T)
  ac#0
  T1 = flip(T,ac#0)
  volumeVector(pts, T1)
  checkFan(pts, T1)
  checkFan(pts, T)
  isRegularTriangulation(pts,T)
  isRegularTriangulation(pts,T1)
  triS = new MutableHashTable from {T=>true}
  ac = select(affineCircuits(pts,T), x -> #x#0 > 1 and #x#1 > 1)
  newT = for a in ac list (t := flip(T,a); if t === null then continue else t)
  Ts = join({T},newT)
  Ts/(t -> volumeVector(pts,t))/sum
  Ts/(t -> elapsedTime checkFan(pts,t))
  Ts/(t -> sort unique flatten t)
  Ts/(t -> isRegularTriangulation(pts,t))
  unique oo
///

TEST ///
-*
  restart
*-
  -- from Kreuzer-Skarke database
  -- 4 9  M:32 9 N:11 8 H:6,30 [-48]
  polystr = "   1   0   1   1  -1   1   0  -1  -2
    0   1   0   0   0  -2  -2   2   2
    0   0   2   0  -2  -2  -2   2   2
    0   0   0   2  -2   2   1   0  -1"

  A = matrixFromString polystr
  P1 = convexHull A
  P2 = polar P1

  elapsedTime LP1 = latticePointList P1
  V1 = vertexList P1
  assert(take(LP1,#V1) == V1)

  elapsedTime LP2 = latticePointList P2
  V2 = vertexList P2
  assert(take(LP2,#V2) == V2)

  elapsedTime assert(faceList(0,P1) == for i from 0 to 8 list {i})
  assert(# faceList(1,P1) == 22)
  assert(# faceList(2,P1) == 21)
  assert(# faceList(3,P1) == 8)

  elapsedTime assert(faceList(0,P2) == for i from 0 to 7 list {i})
  assert(# faceList(1,P2) == 21)
  assert(# faceList(2,P2) == 22)
  assert(# faceList(3,P2) == 9)

-*
  faceList(4,P1) -- what should this do?
  assert(faceList(-1,P1) == {{}}) -- not correct yet
  faceList(4,P2) -- what should this do?
  assert(faceList(-1,P2) == {{}}) -- not correct yet
*-

  -- now test dual faces...    
  dualfaces = for f in faceList(1,P2) list dualFace(P2,f)
  origfaces = for g in dualfaces list dualFace(P1,g)
  assert(origfaces == faceList(1,P2))
  for g in dualfaces do assert(2 == dim(P1,g))
  
  -- now test lattice point containment
  -- for each face, want the lattice points on that face
  faceList(1,P2)  
  for f in faceList(1,P2) list f => latticePointList(P2,f)
  for f in faceList(2,P2) list f => latticePointList(P2,f)
  -- now test interior lattice points
  -- need to run through all lattice points, and find max face containing it
  -- check this against Polyhedra code
  H = latticePointHash P2

  (vertexMatrix P2)_{0}
  Q = convexHull oo
  vertexList Q
  vertexMatrix Q
  
  latticePoints Q
  -- TODO: fix the following bug.  This results when using these functions in cases when 
  -- the polytope isn't e.g. reflexive, or really, doesn't have the origin in the interior...
  -- latticePointList Q -- fails, since polar Q isn't really what should be used here...

  for f in (faceList P2)_{0} do (
      Q := convexHull (vertexMatrix P2)_f;
      lp := (latticePoints Q)/(m -> flatten entries m);
      lpi := latticePointList(P2,f);
      lpi2 := lp/(p -> H#p);
      assert(lpi == sort lpi2)
      )

  H = latticePointHash P1
  for f in faceList P1 do (
      Q := convexHull (vertexMatrix P1)_f;
      lp := (latticePoints Q)/(m -> flatten entries m);
      lpi := latticePointList(P1,f);
      lpi2 := lp/(p -> H#p);
      assert(lpi == sort lpi2)
      )
  
  -- test interior lattice point code
  lpi = (latticePointHash P1)#{-2, 2, 2, -1}
  assert(minimalFace(P1, {-2,2,2,-1}) == {0})
  
  assert(interiorLatticePointList(P1, {0,2,3,5,7}) == {24, 26, 28})

  hashTable for f in faceList(2,P2) list (
      f => {dim(P2,f), 
          latticePointList(P2,f), 
          # interiorLatticePointList(P2,f), 
          # interiorLatticePointList(P1, dualFace(P2,f))}
      )
  allinfo = sort for f in faceList P2 list (
      {dim(P2,f), 
          f, 
          latticePointList(P2,f), 
          # interiorLatticePointList(P2,f), 
          # interiorLatticePointList(P1, dualFace(P2,f))
          }
      )
  allinfo'ans = {
      {0, {0}, {0}, 1, 0}, 
      {0, {1}, {1}, 1, 0}, 
      {0, {2}, {2}, 1, 0}, 
      {0, {3}, {3}, 1, 0}, 
      {0, {4}, {4}, 1, 0}, 
      {0, {5}, {5}, 1, 0}, 
      {0, {6}, {6}, 1, 0}, 
      {0, {7}, {7}, 1, 0}, 
      {1, {0, 1}, {0, 1}, 0, 1}, 
      {1, {0, 2}, {0, 2}, 0, 0}, 
      {1, {0, 3}, {0, 3}, 0, 0}, 
      {1, {0, 5}, {0, 5}, 0, 0}, 
      {1, {0, 6}, {0, 6}, 0, 0}, 
      {1, {1, 2}, {1, 2}, 0, 0}, 
      {1, {1, 4}, {1, 4, 8}, 1, 0}, 
      {1, {1, 5}, {1, 5}, 0, 1}, 
      {1, {1, 6}, {1, 6}, 0, 0}, 
      {1, {1, 7}, {1, 7}, 0, 1}, 
      {1, {2, 3}, {2, 3}, 0, 0}, 
      {1, {2, 4}, {2, 4}, 0, 0}, 
      {1, {2, 5}, {2, 5}, 0, 1}, 
      {1, {3, 4}, {3, 4}, 0, 0}, 
      {1, {3, 6}, {3, 6}, 0, 0}, 
      {1, {4, 5}, {4, 5}, 0, 0}, 
      {1, {4, 6}, {4, 6}, 0, 0}, 
      {1, {4, 7}, {4, 7}, 0, 1}, 
      {1, {5, 6}, {5, 6, 9}, 1, 3}, 
      {1, {5, 7}, {5, 7}, 0, 1}, 
      {1, {6, 7}, {6, 7}, 0, 1}, 
      {2, {0, 1, 2}, {0, 1, 2}, 0, 0}, 
      {2, {0, 1, 3, 4}, {0, 1, 3, 4, 8}, 0, 1}, 
      {2, {0, 1, 5}, {0, 1, 5}, 0, 0}, 
      {2, {0, 1, 6}, {0, 1, 6}, 0, 1}, 
      {2, {0, 2, 3}, {0, 2, 3}, 0, 1}, 
      {2, {0, 2, 5}, {0, 2, 5}, 0, 0}, 
      {2, {0, 3, 6}, {0, 3, 6}, 0, 1}, 
      {2, {0, 5, 6}, {0, 5, 6, 9}, 0, 1}, 
      {2, {1, 2, 4}, {1, 2, 4, 8}, 0, 1}, 
      {2, {1, 2, 5}, {1, 2, 5}, 0, 0}, 
      {2, {1, 4, 7}, {1, 4, 7, 8}, 0, 1}, 
      {2, {1, 5, 6}, {1, 5, 6, 9}, 0, 0}, 
      {2, {1, 5, 7}, {1, 5, 7}, 0, 0}, 
      {2, {1, 6, 7}, {1, 6, 7}, 0, 0}, 
      {2, {2, 3, 4}, {2, 3, 4}, 0, 1}, 
      {2, {2, 3, 5, 6}, {2, 3, 5, 6, 9}, 0, 1}, 
      {2, {2, 4, 5}, {2, 4, 5}, 0, 1}, 
      {2, {3, 4, 6}, {3, 4, 6}, 0, 1}, 
      {2, {4, 5, 6}, {4, 5, 6, 9}, 0, 0}, 
      {2, {4, 5, 7}, {4, 5, 7}, 0, 0}, 
      {2, {4, 6, 7}, {4, 6, 7}, 0, 0}, 
      {2, {5, 6, 7}, {5, 6, 7, 9}, 0, 1}, 
      {3, {0, 1, 2, 3, 4}, {0, 1, 2, 3, 4, 8}, 0, 1}, 
      {3, {0, 1, 2, 5}, {0, 1, 2, 5}, 0, 1}, 
      {3, {0, 1, 3, 4, 6, 7}, {0, 1, 3, 4, 6, 7, 8}, 0, 1}, 
      {3, {0, 1, 5, 6}, {0, 1, 5, 6, 9}, 0, 1}, 
      {3, {0, 2, 3, 5, 6}, {0, 2, 3, 5, 6, 9}, 0, 1}, 
      {3, {1, 2, 4, 5, 7}, {1, 2, 4, 5, 7, 8}, 0, 1}, 
      {3, {1, 5, 6, 7}, {1, 5, 6, 7, 9}, 0, 1}, 
      {3, {2, 3, 4, 5, 6}, {2, 3, 4, 5, 6, 9}, 0, 1}, 
      {3, {4, 5, 6, 7}, {4, 5, 6, 7, 9}, 0, 1},
      {4, {0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 1, 0}
      }
  assert(allinfo == allinfo'ans)
  assert(allinfo == annotatedFaces P2)  

  hodgeOfCYToricDivisors P1
  hodgeOfCYToricDivisors P2
  elapsedTime hashTable for p in latticePointList P2 list p => hodgeCY(P1, p)
  elapsedTime hashTable for p in latticePointList P1 list p => hodgeCY(P2, p)
  
  netList annotatedFaces(0,P1)
  netList annotatedFaces(3,P2)
  netList annotatedFaces(1,P1)
  netList annotatedFaces(2,P1)
  netList annotatedFaces(3,P1)
///

TEST ///
  -- Test of h11 and h21 formulae for an example from Kreuzer-Skarke database.
  -- as well as minimalFace, dim, genus, isFavorable.
-*
  restart
*-
  -- from Kreuzer-Skarke database
  -- 4 9  M:32 9 N:11 8 H:6,30 [-48]
  polystr = "   1   0   1   1  -1   1   0  -1  -2
    0   1   0   0   0  -2  -2   2   2
    0   0   2   0  -2  -2  -2   2   2
    0   0   0   2  -2   2   1   0  -1"

  A = matrixFromString polystr
  P1 = convexHull A
  P2 = polar P1

  elapsedTime assert(h11OfCY P1 == 6)
  elapsedTime assert(h11OfCY P2 == 30)
  elapsedTime assert(h21OfCY P1 == 30)
  elapsedTime assert(h21OfCY P2 == 6)

  assert isFavorable P1  
  assert not isFavorable P2 
  
  for f in faceList P2 list {dim(P2,f), genus(P2,f)}
  minfaces = for lp in latticePointList P2 list dim(P2,minimalFace(P2,lp))
  assert(minfaces == {0,0,0,0,0,0,0,0,1,1,4})
///

TEST ///
  -- We work on one example in 4 dimensions, where we know the answers (or have computed them elsewhere).
  -- Second polytope (index 1) on h11=3 Kreuzer-Skarke list of 4d reflexive polytopes for h11=3.
  -- 4 5  M:48 5 N:8 5 H:3,45 [-84]
  str = "  1   0   2   4  -8
         0   1   5   3  -9
         0   0   6   0  -6
         0   0   0   6  -6
  "
  M = matrixFromString str
  assert(M == matrix {{1, 0, 2, 4, -8}, {0, 1, 5, 3, -9}, {0, 0, 6, 0, -6}, {0, 0, 0, 6, -6}})

  P = convexHull M  
  assert(# latticePointList P == 48)
  assert(# latticePointList polar P == 8)
  assert(h11OfCY P == 3)
  assert(h21OfCY P == 45)
  assert(h11OfCY polar P == 45)
  assert(h21OfCY polar P == 3)
  
  -- Now compute all of the cohomologies of the (irreducible) toric divisors 
  LP = latticePointList polar P
  assert(#LP == 8)
  LP = drop(LP, -1)
  cohoms = for v in LP list hodgeOfCYToricDivisor(P, v)
  cohomH = hodgeOfCYToricDivisors P
  cohoms1 = for v from 0 to #LP-1 list cohomH#v -- last LP is the origin, which isn't one of these divisors
  assert(cohoms == cohoms1)
  assert isFavorable P

  -- Now compute all of the cohomologies of the (irreducible) toric divisors for the polar dual
  LP = latticePointList P
  LP = drop(LP, -1)
  cohoms = for v in LP list hodgeOfCYToricDivisor(polar P, v)
  cohomH = hodgeOfCYToricDivisors polar P
  cohoms1 = for v from 0 to #LP-1 list cohomH#v
  assert(cohoms == cohoms1)
  assert not isFavorable polar P
///

TEST ///
  -- id=0 h11=11
  -- 4 13  M:23 13 N:16 13 H:11,18 [-14]
  str = "    1    0    0    0   -1    1    0    0    0    1   -1    1   -2
    0    1    0    0    1   -1    0    1    1   -1    1   -2    0
    0    0    1    0    1   -1    0    1    0    0   -1   -2    2
    0    0    0    1   -1    1   -1   -1   -1    1    1    0   -1
    "

  P = convexHull matrixFromString str -- last first eg11
  P2 = polar P
  assert(h11OfCY P == 11)
  assert(h21OfCY P == 18)
  assert(h11OfCY polar P == 18)
  assert(h21OfCY polar P == 11)

  -- Now compute all of the cohomologies of the (irreducible) toric divisors 
  LP = latticePointList polar P  
  LP = drop(LP, -1)
  assert(#LP == 15)
  
  --elapsedTime cohoms = for v in LP list hodgeOfCYToricDivisor(P, v)
  cohomH = hodgeOfCYToricDivisors P
  cohoms1 = for v from 0 to #LP-2 list cohomH#v

  -- the following line is too slow, and 
  cohoms = for v in drop(LP,-1) list hodgeOfCYToricDivisor(P, v)
  assert(cohoms == cohoms1)
  assert isFavorable P

  -- Now compute all of the cohomologies of the (irreducible) toric divisors for the polar dual
  LP = latticePointList P
  LP = drop(LP, -1)
  assert(#LP == 22)
  elapsedTime cohoms = for v in LP list hodgeOfCYToricDivisor(polar P, v)
  cohomH = hodgeOfCYToricDivisors polar P
  cohoms1 = for v from 0 to #LP-1 list cohomH#v
  assert(cohoms == cohoms1)
  -- the following line is too slow, and 
--  cohoms = for v in LP list hodgeOfCYToricDivisor(P, v) -- actually, gives an error here...
--  assert(cohoms == cohoms1)
  assert isFavorable polar P
///

TEST ///
-*
restart
*-
str = "    1    0    0    0  -11   -3   -1   -3
    0    1    0    0   -6   -2   -2   -6
    0    0    1    0   -2   -2   -2   -2
    0    0    0    1   -2    2    4    6 "
M = matrixFromString str
P = convexHull M
vertexMatrix P
lp = latticePointList polar P
sort for v in vertexList P list (latticePointHash P)#v

assert try (hodgeOfCYToricDivisor(P,{0,0,0,1}); false) else true -- {0,0,0,1} is not a lattice point
for p in lp list p => hodgeOfCYToricDivisor(P,p)
hodgeOfCYToricDivisors P

assert(h11OfCY(P) == 90)
assert(h21OfCY(P) == 50)
assert(h11OfCY(polar P) == 50)
assert(h21OfCY(polar P) == 90)
///

----------------------------
-- tests from MyPolyhedra --
----------------------------

TEST ///  
  debug needsPackage "StringTorics"
  A = transpose matrix {{-1,-1,2},{-1,0,1},{-1,1,1},{0,-1,2},{0,1,1},{1,-1,3},{1,0,-1},{1,1,-2}}
  tri = regularFineTriangulation A
  volumeVector(augmentWithOrigin A, max tri)
  P = convexHull A

  A = transpose matrix {{-1, 0, -1, -1}, {-1, 0, 0, -1}, {-1, 1, 2, -1}, {-1, 1, 2, 0}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, 0, -1, 2}, {1, 0, 1, 2}}
  C = transpose matrix latticePointList polar convexHull A
  tri = regularFineTriangulation C

  elapsedTime delaunaySubdivision C -- takes 20 seconds?! (now 8 seconds)
  isRegularTriangulation tri
  P2 = polar convexHull A
  regularStarTriangulation P2
  volume P2
///


TEST ///
  A = transpose matrix {{1, 0, 0, 0}, {1, 2, 0, 0}, {1, 0, 2, 0}, {0, 0, 0, 1}, {0, 4, 0, 1}, {0, 0, 4, 1}, {-2, -4, -6, -3}, {-2, -6, -6, -3}, {-2, -6, -4, -3}}
  P = convexHull A
  P2 = polar P

  -- triangulations
  A1 = transpose matrix latticePointList P
  tri = regularFineTriangulation A1
  -- isRegularTriangulation tri -- TODO: this takes forever, is that new?

  debugLevel = 3
  tri = regularStarTriangulation P
  normalToricVariety(drop(latticePointList P,-1), max tri)
  -- elapsedTime isWellDefined oo -- ouch, pretty long (4/27/20) (96 sec)

  -- vertices
  vertexList P
  vertexMatrix P
  vertexList P2
  vertexMatrix P2

  -- latticePoints
  latticePointList P
  latticePointList P2
  
  -- interior lattice points
  debug StringTorics
  hash1 = hashTable for f in faceList P list f => interiorLatticePointList(P, f)
  hash2 = hashTable select(pairs hash1, (k,v) -> #v > 0)
  hash3 = P.cache.TCIInteriorLatticeHash
  assert(hash2 === hash3)

  -- faces
  faceDimensionHash P
  faceList P
  faceList(0, P)
  faceList(1, P)
  faceList(2, P)
  faceList(3, P)
  faceList(4, P)

  FL1 = for f in faceList P list dualFace(P,f)
  FL2 = for f in FL1 list dualFace(polar P,f)  
  FL3 = for f in FL2 list dualFace(P,f)  
  assert(FL1 == FL3)
  assert(faceList P == FL2)

  time for f in faceList P list dim(P,f)
  time for f in faceList P list latticePointList(P,f)
  time for lp in latticePointList P list minimalFace(P,lp)
  time for f in faceList P list genus(P,f)

  netList annotatedFaces P
  for i from 0 to 4 list netList annotatedFaces(i,P)
///



TEST /// -- medium size (h^11 = 15) example
-*
  restart
  needsPackage "StringTorics"
*-
  topes = kreuzerSkarke(15, Limit=>10, Access=>"wget")
  A1 = matrix topes_8
  P = convexHull A1
  P2 = polar P
  A = transpose matrix latticePointList P2  

  elapsedTime tri = regularStarTriangulation P2;
  A2 = transpose matrix first tri
  assert(A2 == submatrix(A, 0..numcols A-2))
  tri = tri_1/(x -> append(x, numcols A - 1))
  isFine(A, tri)
  isStar(A, tri)
  wts = regularTriangulationWeights(A, tri)
  elapsedTime regularSubdivision(A, matrix{wts}) -- this is slower than we would like
  assert(oo == tri) -- both oo, tri should be already sorted.

  circs = affineCircuits(A, tri)
  circs0 = select(circs, x -> not member(numcols A - 1, flatten x))  
  for c in circs0 list bistellarFlip(tri, c)
  bistellarFlip(tri, circs0_1)
  
  elapsedTime tris = generateTriangulations(A, tri, Limit => 50);
  tris/isFine_A//tally
  tris/isStar_A//tally
  elapsedTime(tris/regularTriangulationWeights_A);
  
  elapsedTime tri = regularFineTriangulation A; -- topcom, fast.
  --  elapsedTime tris = allTriangulations(A, Fine => true); -- pretty long, how many are there?

  ans = {{0, 2, 3, 5, 17}, {0, 2, 3, 9, 17}, {0, 2, 5, 15, 17}, {0, 2, 6, 9, 17}, {0, 2, 6, 13, 17}, 
      {0, 2, 7, 13, 17}, {0, 2, 7, 15, 17}, {0, 3, 5, 11, 17}, {0, 3, 9, 11, 17}, {0, 5, 6, 11, 17}, 
      {0, 5, 6, 13, 17}, {0, 5, 13, 15, 17}, {0, 6, 9, 11, 17}, {0, 7, 13, 15, 17}, {1, 2, 3, 9, 17}, 
      {1, 2, 3, 10, 17}, {1, 2, 4, 10, 17}, {1, 2, 4, 12, 17}, {1, 2, 6, 9, 17}, {1, 2, 6, 12, 17}, 
      {1, 3, 5, 10, 17}, {1, 3, 5, 11, 17}, {1, 3, 9, 11, 17}, {1, 4, 5, 10, 17}, {1, 4, 5, 12, 17}, 
      {1, 5, 6, 11, 17}, {1, 5, 6, 12, 17}, {1, 6, 9, 11, 17}, {2, 3, 5, 10, 17}, {2, 4, 5, 10, 17}, 
      {2, 4, 5, 16, 17}, {2, 4, 12, 16, 17}, {2, 5, 7, 15, 17}, {2, 5, 7, 16, 17}, {2, 6, 12, 13, 17}, 
      {2, 7, 12, 13, 17}, {2, 7, 12, 16, 17}, {4, 5, 12, 14, 17}, {4, 5, 14, 16, 17}, {4, 8, 12, 14, 17}, 
      {4, 8, 12, 16, 17}, {4, 8, 14, 16, 17}, {5, 6, 12, 14, 17}, {5, 6, 13, 14, 17}, {5, 7, 13, 14, 17}, 
      {5, 7, 13, 15, 17}, {5, 7, 14, 16, 17}, {6, 8, 12, 13, 17}, {6, 8, 12, 14, 17}, {6, 8, 13, 14, 17}, 
      {7, 8, 12, 13, 17}, {7, 8, 12, 16, 17}, {7, 8, 13, 14, 17}, {7, 8, 14, 16, 17}}
  startri = regularFineStarTriangulation(A, ConeIndex => 17) -- last column of A is the origin
  assert(ans == startri)
///

TEST ///
-- This test is failing: May 2022.  It isn't a complete test anyway...
-- Test of intersection number computations.
-- This requires that V be favorable?
-- Remove this test?  In any case, make sure intersection numbers are bombproof!
-*
  restart
  needsPackage "StringTorics"
*-
  debug StringTorics

  topes = kreuzerSkarke(3, Limit => 50);    
  A = matrix topes_30
  P = convexHull A
  (V, basisElems) = reflexiveToSimplicialToricVarietyCleanDegrees(P, CoefficientRing => ZZ/32003)
  basisElems -- for the moment, we ignore this, and write down all of the elements...
  GLSM = transpose matrix degrees ring V
  X = completeIntersection(V, {-toricDivisor V})
  Xa = abstractVariety(X, base())
  IX = intersectionRing Xa
  elemsToConsider = toList(0..numcols GLSM-1);
  triples = (subsets(elemsToConsider, 3))/sort//sort;
  Htriples = hashTable for a in triples list (
      (i,j,k) := toSequence a;
      val := integral(IX_i * IX_j * IX_k);
      if val == 0 then continue else {i,j,k} => val
    )


  -- Now we try the code above
  (singles, doubles, triples) = toSequence possibleNonZeros V

  Htriples = hashTable for a in join(singles, doubles, triples) list (
      (i,j,k) := toSequence a;
      val := integral(IX_i * IX_j * IX_k);
      if val == 0 then continue else {i,j,k} => val
    )

  -- assert(triples === (keys Htriples)/sort//sort) -- failing.  And it shouldn't be correct anyway...?
  
  elapsedTime tripleProductsCY V
  elapsedTime CY3NonzeroMultiplicities V -- much slower for small h11...
  
  ans = toSequence topologyOfCY3(V, basisElems)
  (h11, h21, C, L) = ans

  -- What was this supposed to do?
  -- hashTable for x in keys H3 list (
  --     if isSubset(x, basisElems) then x => H3#x else continue
  --     )

  (h11, h21, C, L) = toSequence topologyOfCY3(V, basisElems)  
  (h11', h21', C', L') = toSequence topologyOfCY3(V, {0,1,2}, Ring => ring C)  
  A = ring C;
  M = GLSM_{0,1,2}  
  phi1 = map(A, A, flatten entries((M) * transpose vars A))
  L' == phi1 L
  C' == phi1 C
  netList {L, L', phi1 L, phi1 L'}
  
///


///
  -- Favorable h11=5 polytope.
  -- This is too long for a test
-*
  restart
*-  
  needsPackage "StringTorics"
  topes = kreuzerSkarke(5, Limit => 20);
  
  A = matrix topes_3
  P = convexHull A  
  assert isReflexive P
  h11OfCY P == 5
  h11OfCY polar P == 29
  h21OfCY P == 29
  assert isFavorable P
  assert not isFavorable polar P

  P2 = polar P
  vertices P2
  transpose matrix latticePointList P2
  netList annotatedFaces P2  

  V = reflexiveToSimplicialToricVariety P
  classGroup V
  transpose matrix degrees ring V

--  elapsedTime Qs = for tope in topes list cyPolytope(tope);
  elapsedTime Qs = for i from 0 to #topes -1 list cyPolytope(topes#i, ID => i);
  Qs/isFavorable -- takes .1 - .2 seconds per polytope.  Why so long?
  favorables = positions(Qs, isFavorable)

  RZ = ZZ[a,b,c,d,e]

  elapsedTime Xs = flatten for i in favorables list elapsedTime findAllCYs Qs#i;
  XH = hashTable for X in Xs list label X => X;
  topXH = hashTable for k in keys XH list elapsedTime k => topologicalData(XH#k, RZ);
  UtopXH = partition(k -> topXH#k, keys topXH)
  hashTable select(pairs UtopXH, k -> hh^(1,2) k#0 == 39)

///

///
  -- THIS TEST CURRENTLY FAILS (Aug 2022).
  -- Non favorable example.
  -- Either implement functionality for this situation, or give reasonable error messages!
  -- XXX start here Aug 2022.
  -- this is an h11=5 polytope.  Let's make sure everything seems ok with it 
  -- reason: it is seemingly becoming an h11=4 polytope?
  -- Actually: it is a torsion grading.
-*
  restart
*-  
  needsPackage "StringTorics"
  topes = kreuzerSkarke(5, Limit => 50);
  A = matrix topes_1
  P = convexHull A  
  assert isReflexive P
  h11OfCY P == 5
  h11OfCY polar P == 29
  h21OfCY P == 29

  P2 = polar P
  vertices P2
  latticePoints P2
  netList annotatedFaces P2  
  methods reflexivePolytope

  V = reflexiveToSimplicialToricVariety P
  degrees ring V -- fails, 
  classGroup V
  picardGroup V
  h11OfCY P

  Q = cyPolytope topes_1
  -- Q = reflexivePolytope A -- really the dual of A.
  vertices polytope Q
  netList annotatedFaces Q -- annotated faces of the dual of A.
  peek Q.cache
  assert(h11OfCY Q == 5)
  assert(h21OfCY Q == 29)
  assert(dim Q == 4)
  assert not isFavorable Q -- i.e. whether the dual has any points interior to a 2-face, whose dual does too.
  isFavorable polar Q
  X = makeCY Q -- BUG/TODO: should allow CoefficientRing at least...
  V = ambient X -- TODO: need a way to make this directly from Q...
  Xs = findAllFRSTs Q -- only one here, not surprisingly...
  ring V -- BUG: get inscrutable error.  Probably due to a placed GLSM charge matrix...
  
  RZ = ZZ[a,b,c,d]
  topologicalData(X, RZ) -- fails with bad error message.
  abstractVariety X -- fails for similar reason... (bad glsm matrix added...?)
///

///
  -- analyzing the triangulations related to one via bistellar flip.
  -- CURRENT WORK: grabbing all FRST's starting with one.  This test isn't really a test, and currently FAILS.
  -- The method below using bistellar flips seems to work better than topcom.
  -- Although, it is still much slower than it needs to be (don't need to use any circuit twice when moving from one to another...)

  -- my TODO:
  --  keep track of circuits used, and don't use one twice...?
  --  make sure that all but (2,2) flips do not change the 2-faces (I think this is clear, but check anyway).
  --  maybe: get the graph of all such.
-*  
  restart
  needsPackage "StringTorics"
*-
  topes = kreuzerSkarke(6, Limit => 10)  
  Q = cyPolytope(topes_7, ID => 7)
  X = makeCY Q
  assert isFavorable Q
  elapsedTime Xs = findAllCYs Q;
  assert(#Xs == 21)
  PXs = partition(X -> restrictTriangulation X, Xs)
  assert(#keys PXs == 4) -- at most 4 different topologies

  sampleXs = for k in keys PXs list PXs#k#0; -- a list of 4 CY's that have the 4 different topologies.
  sampleXsGV = for X in sampleXs list partitionGVConeByGV(X, DegreeLimit => 20)  
  for p in subsets({0,1,2,3}, 2) list p => findLinearMaps(sampleXsGV#(p#0), sampleXsGV#(p#1))

  -- this shows that 2 of the 4 are likely the same.  We next compute the topology of these 4 X's
  -- and see if that is in fact the case.
  
  -- So now we compute the topologies of the 4 potentially different CY's.
  RZ = ZZ[x_1..x_6]
  elapsedTime topOfXs = sampleXs/(X -> topologicalData(X, RZ));
  debug StringTorics -- FIXME: this should not be needed
  cubics = topOfXs/cubicForm
  c2s = topOfXs/c2

  -- each matrix in Ms determines a topological isomrphism of sampleXs#0 and sampleXs#2:
  Ms = findLinearMaps(sampleXsGV#1, sampleXsGV#3)
  Ms = Ms/(m -> lift(m, ZZ))
  assert all(Ms, m -> det m == 1 or det m == -1)
  phis = for m in Ms list map(RZ, RZ, m)
  for phi in phis do assert(phi(cubics_1) ==  cubics_3 and phi(c2s_1) == c2s_3)

  -- the others appear that they might be different.
  -- in fact, we can show that these are not the same, by looking at
  -- the jacobian locus of each cubic.
  RQ = QQ[gens RZ]
  for f in cubics list (fQ = sub(f, RQ); decompose saturate ideal jacobian fQ)
  netList oo -- shows that 0, 1, 3 are all unique (not related by invertible integral change of basis).




  
  Qs = for i from 0 to #topes-1 list cyPolytope(topes_i, ID => i)
  Qs/isFavorable

  topes = kreuzerSkarke(7, Limit => 10)  
  topes = kreuzerSkarke(8, Limit => 10)

  Qs = for i from 0 to #topes-1 list cyPolytope(topes_i, ID => i)
  Qs/isFavorable
  --A = transpose matrix rays Qs_0
  A = transpose matrix rays Qs_7
  --A = transpose matrix rays Qs_6
  A = A | transpose matrix{{0,0,0,0}}
  t0 = regularFineTriangulation A
  t1 = triangulation(A, fineStarTriangulation(A, max t0, ConeIndex => numcols A - 1))
  isWellDefined t0
  gkzVector t0
  volume convexHull A

  isWellDefined t1
  gkzVector t1
  volume convexHull A

  elapsedTime Ts = allTriangulations(A, Fine => true); -- 387 triangulations this takes quite a while.  Which example has 387?
  elapsedTime aaTs = allTriangulations(A, RegularOnly => false); -- 1278 triangulations
  Ts = Ts/(t -> triangulation(A, t))
  Ts/ isRegularTriangulation //tally -- all 387 are regular (as they should be).
  Ts/isStar//tally
  # (Ts/isFine)
  elapsedTime Xs = findAllCYs Qs_7;
  
  RZ = ZZ[x_1..x_6]
  Q = Qs_7  
  X = cyData(Q, max t1)
  elapsedTime findAllConnectedStarFine t1; -- for h11=7, Qs_7
  elapsedTime findStarFineGraph t1
  annotatedFaces Q
  restrictTriangulation X
  Xs = findAllCYs Q  
  PXs = partition(X -> restrictTriangulation X, Xs)
  X0s = PXs#((keys PXs)#0)
  X0s/(X -> topologicalData(X, RZ))//unique
  X1s = PXs#((keys PXs)#1)
  X1s/(X -> topologicalData(X, RZ))//unique
  X2s = PXs#((keys PXs)#2)
  X2s/(X -> topologicalData(X, RZ))//unique
  X3s = PXs#((keys PXs)#3)
  X3s/(X -> topologicalData(X, RZ))//unique

  gv0 = partitionGVConeByGV(X0s_0, DegreeLimit => 20)
  gv1 = partitionGVConeByGV(X1s_0, DegreeLimit => 20)
  gv2 = partitionGVConeByGV(X2s_0, DegreeLimit => 20)
  gv3 = partitionGVConeByGV(X3s_0, DegreeLimit => 20)  

  partitionGVConeByGV(X0s_0, DegreeLimit => 30)
  partitionGVConeByGV(X1s_0, DegreeLimit => 30)
  partitionGVConeByGV(X2s_0, DegreeLimit => 30)
  partitionGVConeByGV(X3s_0, DegreeLimit => 30)

  findLinearMaps(gv0, gv2)
  findLinearMaps(gv0, gv1)
  findLinearMaps(gv0, gv3)
  findLinearMaps(gv1, gv2)
  findLinearMaps(gv1, gv3)
  findLinearMaps(gv2, gv3)

  T = QQ[t_(1,1) .. t_(6,6)]
  M = genericMatrix(T, 6, 6)

  gv0, gv2
  id1 = trim ideal(M * transpose matrix{gv0#1#0} - transpose matrix{gv2#1#0})
  id2 = trim ideal(M * transpose matrix{gv0#2#0} - transpose matrix{gv2#2#0})
  id4 = (p) -> (
      trim ideal(M * transpose matrix gv0#4 - (transpose matrix gv2#4)_p)
      )
  id8 = (p) -> (
      trim ideal(M * transpose matrix gv0#8 - (transpose matrix gv2#8)_p)
      )

  id4s = for p in permutations 4 list (I := id4 p; if I != 1 then p => I else continue)
  id8s = for p in permutations 2 list (I := id8 p; if I != 1 then p => I else continue)
  ids = flatten for x in id4s list for y in id8s list trim(id1 + id2 + x#1 + y#1)
  ids = select(ids, i -> i != 1)
  Ms = ids/(i -> M % i) -- 4 matrices here!
  for m in Ms list det m
  Ms = for m in Ms list lift(m, ZZ)
  phis = for m in Ms list map(RZ, RZ, m)
  psis = for m in Ms list map(RZ, RZ, m^-1)
  Ms 
  trim(id1 + id2 + id4s#0#1 + id8s#0#1)
  
  T0 = topologicalData(X0s#0, RZ)
  T2 = topologicalData(X2s#0, RZ)
  F0 = cubicForm T0
  F2 = cubicForm T2
  L0 = c2 T0
  L2 = c2 T2

  phis_0 L0 - L2
  phis_0 F0 - F2

  phis_1 L0 - L2
  phis_1 F0 - F2

  phis_2 L0 - L2
  phis_2 F0 - F2

  phis_3 L0 - L2
  phis_3 F0 - F2

  m0 = Ms#0 -- order 2
  m1 = Ms#1 -- order 4
  m2 = Ms#2 -- order 2
  m3 = Ms#3 -- order 4
  for i from 1 to 8 list (m0*m1)^i -- order 2
  for i from 1 to 8 list (m0*m2)^i -- order 2
  for i from 1 to 8 list (m0*m3)^i -- order 2
  for i from 1 to 8 list (m1*m2)^i -- order 2  
  for i from 1 to 8 list (m1*m3)^i -- order 2  
  for i from 1 to 8 list (m2*m3)^i -- order 2  
  for i from 1 to 8 list (m1*m2*m3)^i -- order 2  
  for i from 1 to 8 list (m0*m2*m3)^i -- order 4
  for i from 1 to 8 list (m0*m1*m2*m3)^i -- order 1
  for i from 1 to 8 list m3^i
  
  -- it looks like X0, X2 are the same topology, X1, X3 are not, and so there arr 3 distinct topologies...  
  T1 = topologicalData(X1s#0, RZ)
  T3 = topologicalData(X3s#0, RZ)
  L1 = c2 T1
  F1 = cubicForm T1
  L3 = c2 T3
  F3 = cubicForm T3

  ideal gens gb saturate(ideal F0 + ideal jacobian F0)
  RQ = QQ[gens RZ]
  F0Q = sub(F0, RQ)
  saturate ideal jacobian F0Q
  decompose oo -- 4 singular points of the projective cubic

  F1Q = sub(F1, RQ)
  saturate ideal jacobian F1Q
  decompose oo -- singular locus is quadric in P3 and a point

  F3Q = sub(F3, RQ)
  saturate ideal jacobian F3Q
  decompose oo -- point union 2 points union conic in a plane

///

TEST ///
-*
  restart
  needsPackage "StringTorics"
*-  
  topes = kreuzerSkarke(5, Limit => 10);
  Qs = for i from 0 to #topes-1 list cyPolytope(topes#i, ID => i)
  for tope in topes list isFavorable convexHull matrix tope
  Q = cyPolytope(topes_8, ID => 8)
  Ts = findAllFRSTs Q  
  RZ = ZZ[a,b,c,d,e]
  Xs = findAllCYs(Q, Ring => RZ)

  assert((for X in Xs list label X) === {(8,0)}) -- (X#"polytope data".cache#"id", X.cache#"id") -- id of each example.
  for X in Xs list intersectionNumbers X

  for X in Xs list topologicalData X
  assert(# unique oo == 1)

  Vs = Xs/ambient
  assert all(Vs, isSimplicial)
///  

TEST ///
-- XXX
-*
  restart
  needsPackage "StringTorics"
*-  
  topes = kreuzerSkarke(3, Limit => 50);    
  Q = cyPolytope(topes_30, ID => 30)
  Ts = findAllFRSTs Q
  RZ = ZZ[a,b,c]
  Xs = for i from 0 to #Ts-1 list calabiYau(Q, Ts#i, ID => i, Ring => RZ)
  assert(#Xs == #Ts)
  vertices polytope Q
  label Q
  assert((for X in Xs list label X) === {(30, 0), (30, 1)})
  X = Xs#0
  V = ambient X
  assert isSimplicial V
  assert isProjective V
  intersectionNumbers X
  intersectionNumbersOfCY X
  oo === ooo
  intersectionNumbersOfCY(V, basisIndices Q)

  assert(hh^(1,1) X == 3)
  assert(hh^(1,2) X == 69)

  elapsedTime T = topologicalData X
  hh^(1,1) T
  hh^(1,2) T

  partitionGVConeByGV(X, DegreeLimit => 10)
  partitionGVConeByGV(X, DegreeLimit => 20)
  hilbertBasis gvCone(X, DegreeLimit => 20)
  gv = gvInvariants(X, DegreeLimit => 20);
///  


TEST ///
-- XXX
-*
  restart
  needsPackage "StringTorics"
*-  
  -- Test the routines of this package on the example X given here (h11=3, h12=69)
  topes = kreuzerSkarke(3, Limit => 50);    
  A = matrix topes_30
  P = cyPolytope(topes_30, ID => 30)
  hh^(1,1) P == 3
  hh^(1,2) P == 69
  X = makeCY(P, Ring => (RZ = ZZ[x,y,z]), ID => 0)
  -- findAllFRSTs P
  -- X = cyData(P, first oo, ID => 0)
 
  assert(hh^(1,1) X == 3)
  assert(hh^(1,2) X == 69)
  assert(dim X == 3)
  elapsedTime topologicalData X
  dump X
  dump cyPolytope X
  elapsedTime restrictTriangulation X

  assert(dim X  == 3)
  assert isFavorable X
  rays X
  max X
  V = ambient X -- give the normal toric variety.  Works now, sort of. Problems though: TODO: cache it, allow options? degrees might be different...
  assert(rays V === rays X)
  assert(max V === max X)
  
  intersectionNumbers X  
  toricIntersectionNumbers X
  c2 X
  cubicForm X
  c2Form X
      
  elapsedTime topologicalData X -- cache this result?
  
  ambient X -- give the normal toric variety.  Works now.
  aX = abstractVariety X -- give the abstract variety.  -- TODO: should stash the value...?
  abstractVariety(X, base(a,b,c)) -- give the abstract variety
  IX = intersectionRing aX
  
  -- TODO: How is this computed?
  rays toricMoriCone X
  hilbertBasis toricMoriCone X

  gvInvariants(X, DegreeLimit => 10)

  -- TODO: add tests for line bundles on X, and their cohomology.
///

/// -- NOT TESTED.
-*
  restart
  needsPackage "StringTorics"
*-  
  -- Testing interface for calabiYau, cyPolytope, in presence of databases.
  DB = (currentDirectory) | "Databases/cys-ntfe-h11-5.dbm"  
  RZ = ZZ[x_0..x_4]
  --elapsedTime (Qs, Xs) = readCYDatabase(DB, Ring => RZ); -- takes 13 seconds.
  --# sort keys Xs == 13635
  
  -- reading examples direcly from the database
  db = openDatabase DB
  # sort keys db === 18625
  X = calabiYau(db, (4782,0), Ring => RZ)
  close db
  Q = cyPolytope X
  peek Q.cache
  netList annotatedFaces cyPolytope X

  Q = cyPolytope(DB, 4782)

  db = openDatabase DB
  Xlabs = select(keys db, lab -> (lab = value lab; instance(lab, Sequence) and lab#0 == 4510))
  Xlabs = Xlabs/value
  close db
  X1 = calabiYau(DB, Xlabs#0, Ring => RZ)
  X1' = calabiYau(DB, Xlabs#0, Ring => RZ)
  X1 === X1'
  cyPolytope X1 === cyPolytope X1'

  Q = cyPolytope X1
  X2 = calabiYau(DB, Xlabs#1, Ring => RZ)
  Q2 = cyPolytope X2
  Q === Q2
  assert(Q.cache === Q2.cache)
  automorphisms Q
  netList restrictTriangulation X,  netList restrictTriangulation X2
///

/// -- TODO: this test isn't finding the database file.
-*
  restart
  needsPackage "StringTorics"
*-  
  -- Testing interface for calabiYau, cyPolytope, in presence of databases.
  debug needsPackage "StringTorics"
  -- TODO: this line doesn't work in tests, since the 
  DB = "./Databases/cys-ntfe-h11-3.dbm"
  RZ = ZZ[a,b,c]
  elapsedTime (Qs, Xs) = readCYDatabase(DB, Ring => RZ);
  Q = Qs#6
  Xs = findAllCYs(Q, NTFE => false)
  G = automorphisms Q
  tri1 = restrictTriangulation_2 Xs#0
  gtri1 = normalizeByAutomorphisms(Q, tri1)

  tri2 = restrictTriangulation_2 Xs#0
  gtri2 = normalizeByAutomorphisms(Q, tri2)
  gtri1 === tri1
  gtri2 === tri2
  gtri1 === gtri2
  #Xs
  # unique{gtri1, gtri2}
  # findAllCYs(Qs#31, NTFE => false, Automorphisms => false) == 6
  # findAllCYs(Qs#31, NTFE => false, Automorphisms => true) == 1
  # findAllCYs(Qs#31, NTFE => true, Automorphisms => false) == 6
  # findAllCYs(Qs#31, NTFE => true, Automorphisms => true) == 1
  for lab in sort keys Qs list (
      Q = Qs#lab;
      elapsedTime {# findAllCYs(Q, NTFE => false, Automorphisms => false),
      # findAllCYs(Q, NTFE => false, Automorphisms => true),
      # findAllCYs(Q, NTFE => true, Automorphisms => false),
      # findAllCYs(Q, NTFE => true, Automorphisms => true)
      })
  Q = Qs#31
///
