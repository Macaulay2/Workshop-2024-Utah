-*
  restart
  needsPackage "StringTorics"
*-
/// -- test of: reflexivePolytope, polar.
  -- this example is from this:
  tope = KSEntry "4 6  M:57 6 N:9 8 H:4,52 [-96] id:40
    1   -1    1    1    1   -3
    0    0    2    2    2   -6
    0    0    0    4    0   -4
    0    0    0    0    4   -4
  "

  -- method #1 from a list
  verts = {{-1, 0, 0, 0}, {-1, 0, 0, 1}, {-1, 0, 1, 0},
      {-1, 2, -1, -1}, {1, -1, 0, 0}, {1, -1, 0, 1},
      {1, -1, 1, 0}, {1, 1, -1, -1}}
  Q = reflexivePolytope(verts, ID => 40) -- from h11=4, id=40 Kreuzer Skarke
  assert isReflexive polytope Q
  assert isWellDefined Q

  -- method #2 from a matrix
  Q2 = reflexivePolytope(transpose matrix verts, ID => 40)
  assert isWellDefined Q2
  assert(Q2 === Q)
  Q2.cache === Q.cache -- this should be false?  Here is the actual check:
  assert not (hashTable pairs Q2.cache === hashTable pairs Q.cache)
  assert(vertices Q === vertices Q2)

  -- method #3: using polar
  Q' = polar Q -- It should be that Q stashes Q' (and vice versa?)
  Q3 = polar Q'
  isWellDefined Q3
  assert(Q3 === Q)
  hashTable pairs Q.cache === hashTable pairs Q3.cache
  -- TODO: Q' only stashes the polytope object, not the ReflexivePolytope object.
  --   fix this.

  -- method #4: using a Polyhedron object
  P4 = convexHull transpose matrix verts
  Q4 = reflexivePolytope(P4, ID => 40)
  isWellDefined Q4
  assert(Q4 === Q)

  -- method #5: Using Kreuzer-Skarke database
  Q5 = reflexivePolytope tope -- takes ID from tope.
  isWellDefined Q5
  assert(Q5 === Q)
///

-*
  restart
  needsPackage "StringTorics"
*-
/// -- test of: dump, reflexivePolytope String (inverses of each other
  -- this example is from this:
  verts = {{-1, 0, 0, 0}, {-1, 0, 0, 1}, {-1, 0, 1, 0},
      {-1, 2, -1, -1}, {1, -1, 0, 0}, {1, -1, 0, 1},
      {1, -1, 1, 0}, {1, 1, -1, -1}}
  Q = reflexivePolytope(verts, ID => 40) -- from h11=4, id=40 Kreuzer Skarke
  assert isWellDefined Q

  computeBasics Q
  str = dump Q
  Q2 = reflexivePolytope str
  assert(Q === Q2)
  str2 = dump Q2
  assert(str === str2)
///  

-*
  restart
  needsPackage "StringTorics"
*-
/// -- test of: findAllFRVTs
  topes = kreuzerSkarke 5;
  tope = topes_101
  Q = reflexivePolytope tope
  tri = findOneFRST Q
  computeBasics Q
  elapsedTime # findAllFRSTs Q
  elapsedTime # findAllFRVTs Q
  tris = findAllFRVTs Q
  assert all(for t in tris list (
      V := normalToricVariety(rays Q, t);
      result := isWellDefined V;
      << result << endl;
      result
      ))
///

-*
  restart
  needsPackage "StringTorics"
*-
///
  --topes = kreuzerSkarke(5, Limit => 100)
  --tope = topes_55
  tope = KSEntry "4 11  M:48 11 N:10 8 H:5,41 [-72] id:55
   1   0   0   0  -3   3   1  -3  -4  -1  -4
   0   1   0   1   2  -4   0   2   3   0   3
   0   0   1   0  -1   2   0  -2  -2  -2  -2
   0   0   0   3   3  -6   2   2   4   2   5
  "
  Q = reflexivePolytope tope

  tri = findOneFRST Q

  --findOneFRVT Q -- WRONG!!

  elapsedTime ptris = findAllFRSTs Q; -- relatively expensive
  assert(#ptris == 18)
  for tri in ptris list (
      V := normalToricVariety(rays Q, tri);
      isWellDefined V and isProjective V
      )
  assert all(oo, x  -> x)
  
  vtris = sort findAllFRVTs Q;
  for tri in vtris list (
      V := normalToricVariety(rays Q, tri);
      isWellDefined V and isProjective V
      )
  tally oo -- 18 false, 121 true.
  assert all(oo, x  -> x)
  assert(#vtris == 103)
  assert(#ptris == 18)

  ---------------------------------------------------------
  -- Testing with vector triangulations
  -- Now let's try to maneuver in the space of triangulations
  -- starting with tri above.
  cQ = chirotope(transpose matrix rays Q, Homogenize => false)
  "topcomfoo.in" << toString cQ << endl << close
  ctris = for L in lines get ("!chiro2finetriangs <topcomfoo.in") list (
      matches := regex("\\{\\{[0-9,\\{\\}]*", L);
      if matches === null or #matches != 1 then error "my logic is wrong";
      value substring(matches_0, L)
      );
  ctris = sort for x in ctris list (x/sort//sort);
  #ctris == 139
  toList(set vtris + set ptris - set ctris) -- every one in vtris and ptris is in ctris.
  extras = toList (set ctris - (set vtris + set ptris))
  for tri in extras list (
      V := normalToricVariety(rays Q, tri);
      isWellDefined V and isProjective V
      )

  get "!chiro2circuits <topcomfoo.in"
  circuits cQ
  circs = affineCircuits(transpose matrix rays Q, tri)

  bistellarFlip(tri, circs_3)
  T = triangulation(transpose matrix rays Q, tri)
  affineCircuits(transpose matrix rays Q, tri)
  debugLevel = 1
  newtris = generateTriangulations(T, Limit => 200, Homogenize => false)
  newtris = generateTriangulations(transpose matrix rays Q, tri, Limit => 200, Homogenize => false)
  T = triangulation(transpose matrix rays Q, vtris_0)
  newtris = generateTriangulations(T, Homogenize => false)
  set(newtris/max) - set ptris
  set(newtris/max) - (set vtris + set ptris)
  set(newtris/max) - set ctris
  -- why do the others not show up?
  o28 - vtris
  
  for c in circs list (
      newtri := bistellarFlip(tri, c);
      if newtri === null then continue;
      {c, newtri}
      )
  #ptris
  vtris
  ctris
  -- end of testing with vector triangulations, etc
  ----------------------------------------------------------

  
  vtris = findAllSimplicialFans(transpose matrix rays Q);
  #vtris == 121
  ptris = sort for t in vtris list if isTriangulationOfPolytope(Q,t) then max t else continue;
  assert(tris1 === ptris)

  alltris = findAllSimplicialFans(transpose matrix rays Q, Fine => false);
  #alltris == 369
  allptris = sort for t in alltris list if isTriangulationOfPolytope(Q,t) then max t else continue; -- 30
  
  -- findAllFRSTs works as follows currently:
  -- 1. find all fine regular triangulations of the point set, including origin
  -- 2. take those that are star.  In this case, these are the ones that
  --    all subsets contain the last index.
  -- 3. remove the origin from each simplex.

  -- Another way: compute all vector triangulations
  tris = topcomAllTriangulations(transpose matrix rays Q, Homogenize => false, Fine => true); -- 121 of these
  #tris == 121
  
  topcomAllTriangulations(A1, Fine => true)

  elapsedTime findAllFRSTs Q;
  options  allTriangulations
  options topcomAllTriangulations
  A = transpose matrix rays Q
  topcomAllTriangulations A 
  topcomAllTriangulations(A, Fine => true)
  topcomAllTriangulations(A, Homogenize => false, Fine => true)

  dump Q
  
///


-*
  restart
  needsPackage "StringTorics"
*-
TEST ///
  --topes = kreuzerSkarke(10, Limit => 1000);
  -- this is topes_57 (which has interior facet lattice points).
  tope = KSEntry "4 11  M:29 11 N:17 10 H:10,22 [-24] id:57
     1   0   0   0   0   3   3  -1  -2  -3  -4
     0   1   0   0   1  -2  -3   2   2   0   1
     0   0   1   0   0  -1  -1  -1  -1   0   0
     0   0   0   1   1  -1  -1  -1  -1   1   1
  "  

  Q = reflexivePolytope tope

  -- we need to exert some care that
  -- vertices, latticePoints, rays, annotatedFaces all refer to the same indices!
  elapsedTime assert(
      latticePoints Q
      ===
      {{-1, -1, -1, 0}, {-1, -1, -1, 1}, {-1, -1, 0, 0}, {-1, 0, -1, -1},
          {0, -1, -1, 0}, {0, 0, -1, -1}, {0, 0, -1, 2}, {0, 0, 2, -1},
          {0, 1, -1, -1}, {1, 1, -1, 2}, {0, 0, 0, -1}, {0, 0, 0, 1},
          {0, 0, 1, -1}, {0, 0, 1, 0}, {0, 0, -1, 0}, {0, 0, -1, 1},
          {0, 0, 0, 0}}
      )
  elapsedTime assert(
      rays Q
      ===
      (latticePoints Q)_{0..13}
      )
  elapsedTime assert(
      vertices Q
      ===
      (latticePoints Q)_{0..9}
      )
  assert(faceDimensions Q
      ===
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 3, 3, 4}
      )
  assert(faceDimension(Q, 10) == 1)

  -- the following are cached and computed together
  -- degrees, basisIndices, toricBasisIndices
  elapsedTime degrees Q -- currently a bit expensive
  elapsedTime transpose matrix degrees Q
  assert(
      basisIndices Q
      ===
      {0, 1, 2, 3, 4, 5, 6, 7, 10, 11}
      )

  -- the following are cached and computed together
  elapsedTime assert(hh^(1,1) Q == 10)
  assert(hh^(1,2) Q == 22)
  assert isFavorable Q

  -- the following are cached and computed together
  -- automorphisms, autPermutations
  -- TODO: add these back in
  elapsedTime automorphisms Q -- in this case, just the identity
  debug StringTorics
  elapsedTime automorphismsAsPermutations Q

  -- point and vector configurations
  -- TODO.
  -- elapsedTime findAllFRSTs Q; -- too many for this h11=9 example??
  -- options  allTriangulations
  -- options topcomAllTriangulations
  -- A = transpose matrix rays Q
  -- topcomAllTriangulations A -- 105
  -- topcomAllTriangulations(A, Fine => true) -- 75
  -- topcomAllTriangulations(A, Homogenize => false, Fine => true) -- 121
  -- dump Q
///
