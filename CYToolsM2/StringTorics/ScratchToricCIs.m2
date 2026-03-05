-- Scratch code from ToricCompleteIntersections.m2
-- Exploring compatibility between CompleteIntersectionInToric and CalabiYauInToric

///
 -- how compatible with CalabiYauInToric is this?
  -- I guess we need to know if the triangulation comes from a triangulation of the polytope...

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

  -- we need a function that determines if this is a Batryev CY3.
  -- and to return the corresponding CalabiYauInToric...
  transpose matrix rays V
  P = convexHull oo
  vertices P
  Q = cyPolytope P
  rays Q
  X1 = first findAllCYs Q
  hodgeDiamond X1
  hodgeDiamond X
  hh^(1,2) X1 == 69

  -- TODO: add back in once DanilovKhovanskii is in the system
  -- needsPackage "DanilovKhovanskii"
  -- computeHodgeDeligne (-toricDivisor V) -- this is not the correct answer I think!
  -- oo#1
  -- matrix for i from 0 to dim X list for j from 0 to dim X list (-1)^(i+j) * oo#(i,j)
  -- assert(hh^(1,2) Q == 69)
  -- assert(hh^(1,1) Q == 3)
///
