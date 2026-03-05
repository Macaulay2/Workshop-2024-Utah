TEST ///
-*
  restart
  needsPackage "StringTorics"
*-
  vs = {{-1, -1, -1, 0}, {-1, -1, 0, 0}, {-1, -1, 1, -1}, {-1, 0, -1, 0}, {0, -1, 2, -1}, {0, 0, -1, 0}, {0, 1, -1, 0}, {1, 1, -1, 1}, {1, 1, 0, 1}}
  cones4 = {{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 7}, {0, 1, 4, 7}, {0, 2, 3, 5}, {0, 2, 4, 5}, {0, 3, 5, 7}, {0, 4, 5, 7}, {1, 2, 3, 8}, {1, 2, 4, 8}, {1, 3, 7, 8}, {1, 4, 7, 8}, {2, 3, 5, 6}, {2, 3, 6, 8}, {2, 4, 5, 6}, {2, 4, 6, 8}, {3, 5, 6, 7}, {3, 6, 7, 8}, {4, 5, 6, 7}, {4, 6, 7, 8}}
  Q = reflexivePolytope(vs, ID => 1000)
  rays Q == vs
  X = calabiYau(Q, cones4, ID => 0)
  assert(rays X == vs)
  assert(max X == cones4)

  debug needsPackage "StringTorics" -- for toRingElement??  TODO: export that?
  elapsedTime intersectionNumbers X
  toRingElement(oo, picardRing X)
  elapsedTime toricIntersectionNumbers X
  assert(intersectionNumbers X === intersectionNumbersOfCY X)
  elapsedTime c2 X
  c2Form X
  cubicForm X

  elapsedTime intersectionNumbers X
  elapsedTime intersectionNumbersOfCY X

  elapsedTime topologicalData X
///

TEST ///
-- test of the (currently internal) routines: exponentsToProduct,
-- productToExponents, multinomial, toCOO, toRingElement.
  debug StringTorics
  assert(exponentToProduct {} == {})
  assert(exponentToProduct {3} == {0, 0, 0})
  assert(exponentToProduct {0, 3, 1} == {1, 1, 1, 2})
  assert(exponentToProduct {1, 2, 1, 0, 1, 0, 0, 0} == {0, 1, 1, 2, 4})
  assert(exponentToProduct {1, 1, 1, 1, 1} == {0, 1, 2, 3, 4})

  assert(productToExponents({}, 0) == {})
  assert(productToExponents({}, 3) == {0, 0, 0})
  assert(productToExponents({0, 0, 0}, 1) == {3})
  assert(productToExponents({1, 1, 1, 2}, 3) == {0, 3, 1})
  assert(productToExponents({0, 1, 2, 3, 4}, 5) == {1, 1, 1, 1, 1})
  assert(productToExponents({0, 1, 1, 2, 4}, 8) == {1, 2, 1, 0, 1, 0, 0, 0})

  assert(multinomial {3, 0, 0} == 1)
  assert(multinomial {1,0,2} == 3)
  assert(multinomial {1, 1, 1} == 6)

  RZ = ZZ[a,b,c]
  F = (a+2*b+3*c)^3
  G = toCOO F
  F' = toRingElement(G, RZ)
  assert(F == F')
  G' = toCOO F'
  assert(G === G')

  L = 3*a+c
  toCOO L
  assert(toRingElement(toCOO L, RZ) == L)

  L = 1_RZ
  toCOO L
  assert(toRingElement(toCOO L, RZ) == L)

  L = 0_RZ
  toCOO L
  assert(toRingElement(toCOO L, RZ) == L)
///

TEST ///
-- As it turns out, 'monoms' is much faster than first creating the basis,
-- and applying exponentToProduct to (the exponent vector of) every monomial
-- e.g. on MES's Apple M1 Max, 2022, doing nv = 81 the latter way gives .32 + 1.6 seconds
-- instead of .25 seconds.
  debug StringTorics
  elapsedTime assert(# monoms(3, 0, 10) == binomial(13,3))
  elapsedTime assert(# monoms(3, 0, 12) == binomial(15,3))
  elapsedTime assert(# monoms(3, 0, 20) == binomial(23,3))
  elapsedTime assert(# monoms(3, 0, 50) == binomial(53,3))
  elapsedTime assert(# monoms(3, 0, 80) == binomial(83,3)) -- .25 seconds
  elapsedTime assert(# monoms(3, 0, 200) == binomial(203,3)) -- 1.5 seconds

  -- commented out so 'check' doesn't take too long
  --elapsedTime assert(# monoms(3, 0, 300) == binomial(303,3)) -- 5.2 seconds
  --elapsedTime assert(# monoms(3, 0, 400) == binomial(403,3)) -- 14.2 seconds
  --elapsedTime assert(# monoms(3, 0, 495) == binomial(498,3)) -- 31 seconds
  --elapsedTime assert(# monoms(3, 0, 490) == binomial(493,3)) -- 36 seconds, why longer?

  RZ = ZZ[t_1..t_20]
  elapsedTime B = flatten entries basis(3, RZ);
  #B
  mons1 = B/(b -> exponentToProduct first exponents b)
  mons2 = monoms(3, 0, 19)
  mons1 === mons2 -- in the same order

  -- nv = 81 gives the timing above.
  -- the order should be the same,  For testing, we use a smaller value of nv.
  nv = 10
  RZ = ZZ[t_1..t_nv]
  elapsedTime B = flatten entries basis(3, RZ);
  assert(#B == binomial(nv+2, 3))
  elapsedTime mons1 = B/(b -> exponentToProduct first exponents b);
  elapsedTime mons2 = monoms(3, 0, nv-1);
  assert(mons1 === mons2) -- in the same order
///

TEST ///
  -- Let's test the basis intersection numbers code at slightly higher h11...
  -- TODO: This fails, as it uses old naming...
-*
  restart
  needsPackage "StringTorics"
*-
  -- ks was obtained via:
  --  topes = kreuzerSkarke(7, Limit => 50);
  --  assert(#topes == 50)
  --  ks = topes_30
  -- but this requires a network connection, so we don't do that here.
  -- Here it is:
  ks = KSEntry "4 10  M:33 10 N:12 8 H:7,29 [-44] id:30
   1   0   0   1  -1   1  -1  -1  -2   0
   0   1   1   0  -2   0  -2   2  -1   1
   0   0   2   0  -4   2  -2   2  -2   2
   0   0   0   2  -2   2  -2   0  -2   2
   "
  Q = cyPolytope ks
  X = makeCY Q

  -- TODO: reinstate these?
  --toricMoriConeCap X
  --classifyExtremalCurves X

  V = ambient X
  assert isWellDefined V
  assert isProjective V
  assert isSimplicial V

  debug StringTorics
  coo = intersectionNumbers X
  RZ = ZZ[t_0..t_6]
  F = toRingElement(coo, RZ)
  assert(sort coo === sort toCOO F)
///

TEST ///
  -- Let's see how high we can go with this simplistic routine.
-*
  restart
  debug needsPackage "StringTorics"
*-
  debug needsPackage "StringTorics" -- for toRingElement
  --  h11 = 20
  --  topes = kreuzerSkarke(h11, Limit => 50);
  --  assert(#topes == 50)
  --  ks = topes_25
  ks = KSEntry "4 8  M:16 8 N:27 11 H:20,12 [16] id:25
    1    0    0    2    0    0   -1   -2
    0    1    0    0    2    0   -1   -2
    0    0    1    1    1    0   -2   -1
    0    0    0    0    0    1   -1    0
  "

  A = matrix ks
  P1 = convexHull A
  P2 = polar P1
  annotatedFaces P2
  elapsedTime P = reflexivePolytope ks -- reflexivePolytope A
  assert isFavorable P
  assert(hh^(1,1) P == 20)
  assert(hh^(1,2) P == 12)

  elapsedTime X = makeCY P
  elapsedTime coo = intersectionNumbers X; -- now quite fast.  Previously: 4 seconds at h11=20.  3.2 seconds of this is computing the intersection ring.
  assert(#coo == 175)

  RZ = ZZ[t_0..t_(hh^(1,1) P - 1)]
  F = toRingElement(coo, RZ)
  assert(sort coo === sort toCOO F)
///
