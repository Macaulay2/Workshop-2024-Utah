TEST ///
-*
  restart
*-
  debug needsPackage "StringTorics" -- remove 'debug'
  -- Let's consider smooth toric surfaces.
  for i from 0 to 4 list (
      V = smoothFanoToricVariety(2, 1);
      picardGroup V
      )
  V = smoothFanoToricVariety(2, 2)
  transpose matrix degrees ring V
  X = completeIntersection(V, {-toricDivisor V}, Basis => {0,3}, Variables => {symbol a, symbol b})
  linearForm X
  intersectionForm X
  linearForm X
  hh^* OO_X(0,0) -- X is an elliptic curve
  hh^*(OO_X(-1,1))

  L = OO_X(-1,1)
  assert(hh^0(L) == 1)
  assert(hh^1(L) == 0)
  assert(hh^* L == {1,0})
  intersectionForm X
  chern_1 tangentBundle abstractVariety X
  c2Form X
///

TEST ///
-*
restart
*-
  debug needsPackage "StringTorics" -- remove 'debug'

  -*    -- code to generate this example
    topes = kreuzerSkarke(3, Limit => 50);
    A = matrix topes_30
    P = convexHull A
    (V, basisElems) = reflexiveToSimplicialToricVarietyCleanDegrees(P, CoefficientRing => ZZ/32003)
  *-

  verts = {{-1, -1, 0, 0}, {-1, -1, 0, 1}, {-1, -1, 2, 0}, {-1, 0, 0, 0}, {1, -1, -1, 1}, {1, 2, -1, -1}, {-1, -1, 1, 0}}
  maxcones = {{0, 1, 3, 4}, {0, 1, 3, 6}, {0, 1, 4, 6}, {0, 3, 4, 5}, {0, 3, 5, 6}, {0, 4, 5, 6}, {1, 2, 3, 5}, {1, 2, 3, 6}, {1, 2, 4, 5}, {1, 2, 4, 6}, {1, 3, 4, 5}, {2, 3, 5, 6}, {2, 4, 5, 6}}
  V = normalToricVariety(verts, maxcones, CoefficientRing => ZZ/32003)
  glsm = transpose matrix degrees ring V
  basiselems = {0, 5, 6}
  X = completeIntersection(V, {-toricDivisor V}, Basis => {0, 5, 6}, Variables => {symbol a, symbol b, symbol c})
  assert(dim X == 3)
  intersectionRing X
  h = linearForm X
  integral(h^3)
  assert(intersectionForm X == a^3-6*a^2*b+6*a*b^2+12*b^3-9*a^2*c+24*a*b*c+21*a*c^2-24*b*c^2-16*c^3)
  assert(c2Form X == 10*a + 60*b + 8*c)

  X = completeIntersection(V, {-toricDivisor V}, Basis => {0, 5, 6}, Variables => symbol a)

  D = completeIntersection(V, {-toricDivisor V, V_0}, Basis => {0, 5, 6}, Variables => {symbol a, symbol b, symbol c})
  dim D == 2
  saturate(ideal equations D, ideal V)
  intersectionRing D
  intersectionForm D
  c2Form D
  chern_2 tangentBundle abstractVariety D

  L = OO_X(1,1,2)
  assert(hh^* L == {8, 17, 0, 0}) -- TODO: recheck these numbers!
  assert(variety L === X)
  assert(degree L == {1,1,2})

  -- TODO: ADD BACK IN once DanilovKhovanskii more functional
  -- needsPackage "DanilovKhovanskii"
  -- computeHodgeDeligne(-toricDivisor V)
  -- oo#1
  -- matrix for i from 0 to dim X list for j from 0 to dim X list (-1)^(i+j) * oo#(i,j)

  -- want to be able to turn this into a CalabiYauInToric...
  -- then we can check computations against each other too.

///

TEST ///
-*
restart
*-
  debug needsPackage "StringTorics" -- remove 'debug'
  V = kleinschmidt(3, {2,1}, CoefficientRing => ZZ/101)
  rays V
  max V
  picardGroup V
  isSmooth V
  transpose matrix degrees ring V
  X = completeIntersection(V, {-toricDivisor V}, Basis => {3,0}, Variables => {a,b})
  dim X == 2
  hh^*(OO_X(0,0)) == {1, 0, 1}
  F = first equations X
  saturate(ideal F + ideal jacobian F, ideal V) -- X is smooth
  -- hh^*(OO_V(-1,3)) -- ouch!  needs to work...
///
