-- also defined: 
--  (1) dim(P,f)
--  (2) genus(P,f)

-- functions from Polyhedra which we use here:
--   latticePoints P -- a different order than 'latticePoints P'
--   vertices P -- a different order than 'vertexList P' !
--   faces(d,P) -- only in faceDimensionHash
--   polar P -- used often
--   dim P -- used often

-- These are not exported yet...  They should appear elsewhere?
boundedLatticePoints = method(Options => {
        Bound => null,
        Limit => null,
        NodeLimit => 100000 -- set higher if you get an error.
        })
boundedLatticePoints(Matrix, Matrix) := Matrix => opts -> (A, b) -> (
    -- Input: A (m x n) represents a cone {x : Ax <= b} in RR^n, A, b, are integral
    -- Ouput: B (n x N) N is the number of lattice points, B is over ZZ.
    matrix map(ZZ, rawLatticePoints(raw A, raw b, opts.Bound, opts.Limit, opts.NodeLimit))
    )

normalizLatticePoints = method()
-- (A, b) represents the polytope Ax <= b, I believe (that is, we translate this to what normaliz uses.
normalizLatticePoints(Matrix, Matrix) := Matrix => (A, b) -> (
    matrix map(ZZ, rawLatticePointsNormaliz(raw A, raw b))
    )



-- Our plan: stash into a Polyhedron, the information here
protect TCILatticePointList
protect TCIVertexList
protect TCIVertexMatrix
protect TCIFaceDimensionHash
protect TCIInteriorLatticeHash
protect TCILatticePointHash

vertexList = method()
vertexList Polyhedron := (cacheValue symbol TCIVertexList) (P -> (
    verts := vertices P;
    verts = try lift(verts,ZZ) else verts;
    --if liftable(verts,ZZ) then verts = lift(verts,ZZ);
    sort entries transpose verts
    ))

vertexMatrix = method()
vertexMatrix Polyhedron := (cacheValue symbol TCIVertexMatrix) (P -> (
    transpose matrix vertexList P
    ))

faceDimensionHash = method()
faceDimensionHash Polyhedron := (cacheValue symbol TCIFaceDimensionHash) ((P) -> (
    L := vertexList P;
    M := vertices P; -- different ordering, possibly, and also a matrix over QQ
    M = try lift(M,ZZ) else M;
    --if liftable(M,ZZ) then M = lift(M,ZZ);
    verticesQ := entries transpose M;
    vertexHash := hashTable for i from 0 to #L - 1 list (L#i => i);
    hashTable flatten for i from 0 to dim P list for f in faces(dim P-i,P) list (
        -- USING INFO FROM POLYHEDRA
        verts := first f;
        if #(last f) > 0 then error "expected a polytope, but received a polyhedron";
        newverts := sort for v in verticesQ_verts list vertexHash#v;
        newverts => i
        )
    ))

dim(Polyhedron, List) := (P,f) -> (faceDimensionHash P)#f

faceList = method()
faceList(ZZ,Polyhedron) := (dimF,P) -> (
    H := faceDimensionHash P;
    sort select(keys H, f -> H#f === dimF)
    )
faceList Polyhedron := P -> (
    H := faceDimensionHash P;
    (pairs H)/((k,v) -> (v,k))//sort/last
    )

dualFace = method()
dualFace(Polyhedron, List) := (P,f) -> (
    -- f is a sorted list of integer indices of the vertices, giving a face of P
    -- returns a face of 'polar P', as a list of integer indices of the vertices of 'polar P'
    V1 := vertexMatrix P;
    V2 := transpose vertexMatrix polar P;
    vp := V1_f;
    g := positions(entries (V2 * vp), 
        x -> all(x, x1 -> x1 == -1));
    assert(sort g === g);
    g
    )

-- pts is a list of lattice points in P, or a single lattice point.
-- returns the minimal face of P (as sorted list of vertex indices) containing pts
minimalFace = method()
minimalFace(Polyhedron, List) := (P, pts) -> (
    if all(pts, x -> instance(x,ZZ)) then pts = {pts};
    V2 := transpose vertexMatrix polar P;
    dualfaceverts := positions(entries (V2 * transpose matrix pts), 
        x -> all(x, x1 -> x1 == -1));
    dualFace(polar P, dualfaceverts)
    )

latticePointList = method()
latticePointList Polyhedron := (cacheValue symbol TCILatticePointList) (P -> (
    -- this reorders the lattice points via dimension of smallest face containing them
    Q := P; -- use Q for calling functions in Polyhedra, just for doc...
    lp := latticePoints Q;
    -- the following is our putative list of lattice points.
    L := sort for p in lp list flatten entries lift(p,ZZ);
    -- now we reorder this list, and set:
    --   TCILatticePointList
    --   TCILatticePointHash
    --   TCIInteriorLatticeHash
    L1 := sort for i from 0 to #L-1 list (
        f := minimalFace(P, L#i);
        {dim(P,f), L#i, f}
        );
    --L1 := sort for i from 0 to #L-1 list {faceDim(P,minimalFace(P,L#i)),L#i};
    L2 := L1/(x -> x#1); -- select the actual lattice point
    -- now set the hash tables:
    -- lattice point => index
    H0 := hashTable for i from 0 to #L2-1 list L2#i => i;
    P.cache.TCILatticePointHash = H0;
    -- face => interior lattice points
    H1 := partition(x -> last x, L1); -- partition on minimal face
    H2 := applyPairs(H1, (k,v) -> (k,v/(v1 -> H0#(v1#1))));
    P.cache.TCIInteriorLatticeHash = H2;
    L2
    ))
latticePointList(Polyhedron, List) := (P,f) -> (
    -- f is a sorted list of integer indices of the vertices, giving a face of P
    -- returns the list of indices of lattice point on f.
    V1 := vertexMatrix P;
    LP1 := matrix latticePointList P;
    V2 := vertexMatrix polar P;
    g := dualFace(P,f);
    positions(entries(LP1 * V2_g), x -> all(x, x1 -> x1 == -1))
    )

latticePointHash = method()
latticePointHash Polyhedron := P -> (
    latticePointList P;
    P.cache.TCILatticePointHash
    )

-- f is a sorted list of integer indices of the vertices, giving a face of P
-- returns the list of indices of lattice points in the relative interior of f.
interiorLatticePointList = method()
interiorLatticePointList(Polyhedron, List) := (P,f) -> (
    L := latticePointList P;
    H := P.cache.TCIInteriorLatticeHash;
    if H#?f then H#f else {}
    )

-- number of interior points in the dual face
genus(Polyhedron, List) := (P,f) -> # interiorLatticePointList(polar P, dualFace(P,f))

annotatedFaces = method()
annotatedFaces Polyhedron := List => (P1) -> (
    P2 := polar P1;
    sort for f in faceList P1 list (
        {dim(P1,f), 
            f, 
            latticePointList(P1,f), 
            # interiorLatticePointList(P1,f), 
            # interiorLatticePointList(P2, dualFace(P1,f))
            }
      )
    )
-- Returns a list for each face of P1 of dimension i:
-- {faceIndices, all lattice pts, #interior lattice pts, #interior lattice pts in dual face of P2}
annotatedFaces(ZZ,Polyhedron) := List => (i,P1) -> (
    P2 := polar P1;
    sort for f in faceList(i,P1) list (
        {f, 
            latticePointList(P1,f), 
            # interiorLatticePointList(P1,f), 
            # interiorLatticePointList(P2, dualFace(P1,f))
            }
      )
    )

latticePointsAndDimensions = method()
latticePointsAndDimensions Polyhedron := P2 -> (
    LP := latticePointList P2;
    LPdim := for lp in LP list dim(P2, minimalFace(P2, lp));
    (LP, LPdim))


-- private function for isomorphisms
partialPermutations = (elems, num) -> (
    if num == 1 then return elems/(a -> {a});
    flatten for i from 0 to #elems-1 list for p in partialPermutations(drop(elems,{i,i}), num-1)
      list
        prepend(elems#i, p)
    )

isomorphisms = method()

isomorphisms(Polyhedron, Polyhedron, List, List) := (P, Q, annotatedFacesP, annotatedFacesQ) -> (
    nrows := numrows vertexMatrix P;
    if nrows != dim P or nrows != dim Q or nrows != numrows vertexMatrix Q
    then error "expected polytoeps to be full dimensional and same dimension";
    
    -- Step 1. Find a facet of P with the smallest size.
    facetsP := for f in annotatedFacesP list if f#0 != nrows-1 then continue else f#1;
    minsizeP := facetsP/length//min;
    facetsMinsizeP := select(facetsP, f -> #f === minsizeP);
    facetA := first facetsMinsizeP;

    -- Step 2. Find all facets of Q with this smallest size minsizeP, or return {}.
    facetsQ := for f in annotatedFacesQ list if f#0 != nrows-1 then continue else f#1;
    minsizeQ := facetsQ/length//min;
    if minsizeQ =!= minsizeP then return {};
    facetsMinsizeQ := select(facetsQ, f -> #f === minsizeP);

    -- Now find a subset of nrows elements if facetA which are full dimensional
    if #facetA > nrows then (
        -- we need to take a subset of these of size nrows that have full rank.
        -- we then call these facetA again.  We don't actually need facetA again,
        -- the only thing we use is Ainv.
        C := ((vertexMatrix P)_facetA) ** QQ;
        facetA = facetA _ (columnRankProfile mutableMatrix C);
        if #facetA != nrows then error "my logic is missing a case";
        );
    A := (vertexMatrix P)_facetA;
    Ainv := (A ** QQ)^-1;

    -- now we loop through all possible maps from facetA to other facets,
    -- and if it gives an integer matrix, we add it to the list.
    vertsP := (vertexList P)/(v -> transpose matrix {v});
    vertsQ := (vertexList Q)/(v -> transpose matrix {v});
    hashQ := hashTable for i from 0 to #vertsQ-1 list vertsQ#i => i;
    isos := flatten for f in facetsMinsizeQ list (
        for perm in partialPermutations(f, nrows) list (
            B := (vertexMatrix Q)_perm;
            M := B * Ainv;
            try (M = lift(M, ZZ)) else continue;
            if all(vertsP, v -> hashQ#?(M * v)) then M else continue
            )
        );
    isos
    )

isomorphisms(Polyhedron, Polyhedron) := (P, Q) ->
    isomorphisms(P, Q, annotatedFaces P, annotatedFaces Q)

automorphisms = method()
automorphisms Polyhedron := P -> isomorphisms(P, P)

///
-- This is code to test where things stand for speed.
-- TODO: update this to include the latest lattice point code.
  restart
  debug needsPackage "StringTorics"
  topes = kreuzerSkarke 3;

  P1s = elapsedTime for tope in topes list elapsedTime (
        convexHull matrix tope
        ); -- .25s

  P2s = elapsedTime for P1 in P1s list elapsedTime (
        polar P1
        ); -- .7s

  LPs = elapsedTime for P2 in P2s list elapsedTime (
        latticePointList P2
        ); -- 29.6s

  LPs = elapsedTime for P2 in P2s list elapsedTime (
        latticePoints P2
        ); -- 9 sec

  iLP2s = elapsedTime for P2 in P2s list elapsedTime (
        interiorLatticePoints P2
        ); -- 6 sec
    
  LP1s = elapsedTime for P1 in P1s list elapsedTime (
        latticePoints P1
        ); -- 20 sec.

  iLP1s = elapsedTime for P1 in P1s list elapsedTime (
        interiorLatticePoints P1 -- this is expensive!
        ); -- 84.9 sec
    
    
  Qs = elapsedTime for P2 in P2s list elapsedTime (
        cyPolytope P2
        ); -- .45 s

  -- here we put these together
  Qs = elapsedTime for tope in topes list elapsedTime (
        cyPolytope tope
        ); -- 33.2s
    

  debug Polyhedra    
  peek Qs_0 .cache
  peek Qs_0 .cache#"N polytope".cache
  peek Qs_0 .cache#"N polytope".cache.computedPolar.cache

  
  Ps = elapsedTime for Q in Qs list elapsedTime (
        polytope(Q, "N")
        ); -- .23s

  Ps = elapsedTime for Q in Qs list elapsedTime (
        latticePointList polytope(Q, "N")
        ); -- .2s
    
  aF = elapsedTime for Q in Qs list elapsedTime (
        annotatedFaces Q
        ); -- 76.4s, or 57.4s  STILL PRETTY LONG...

  autQs = elapsedTime for Q in Qs list elapsedTime (
        P := polytope(Q, "N");
        automorphisms P
        ); -- 62.0s after annotatedFaces.  Why so slow?
    
  for tope in topes list elapsedTime (
      Q := cyPolytope tope;
      P := polytope(Q, "N");
      ans := isomorphisms(P, P);
      if ans === null then (
        << "--- tope: " << label Q << " TOO LARGE FOR NOW" << endl
        )
      else
        << "--- tope: " << label Q << " #aut=" << #ans << " auts: " << netList ans << endl;
      ans
      );

  Q = cyPolytope topes_20
  Q = cyPolytope topes_0
  P = polytope(Q, "N")
  vertexList P
  annotatedFaces P
  isomorphisms(Q,Q)
///

-*
restart
needsPackage "StringTorics"
*-
TEST /// -- lattice points of a polytope, Nate MacFadden's code.
  H = matrix{{1, 1}, {-1, 0}, {0, -1}}
  rhs = transpose matrix {{1,1,1}}
  LP0 = boundedLatticePoints(H, rhs, Bound => 10, Limit => 10000)

  -- check this against Polyhedra
  P = polyhedronFromHData(H, rhs)
  LPofP = matrix{latticePoints P}
  assert(set entries transpose LP0 === set entries transpose LPofP)
///  

-*
restart
needsPackage "StringTorics"
*-
TEST ///
  ks = KSEntry "4 15  M:41 15 N:11 9 H:6,36 [-60] id:100
         1   0   0   0   0   0   2  -1  -1  -2   2   1  -2  -2   0
         0   1   0   0   0   0  -1   2  -1   2  -2  -2   1   0  -2
         0   0   1   1  -1  -1  -1  -1   2  -1  -1   0   0   1   1
         0   0   0   2  -2   0  -1  -1   2   0  -2  -2   2   2   0
         "
  Q = reflexivePolytope ks
  transpose matrix vertices Q
  Q' = polar Q

  -- objects to test and compare:
  -- Player #1: polytope polar to ks, with 11 vertices, 15 halfplane constraints
  P1 = polytope Q
  (first facets P1) * (vertices P1)
    assert(max flatten entries oo == 1)
  assert(numrows first facets P1 == 15)

  -- the following is tested against Polyhedra below.
  (A, b) = facets P1
  A = lift(A, ZZ)
  rhs = lift(b, ZZ)
  LP0 = boundedLatticePoints(A, rhs, Bound => 10, Limit => 10000)
  
  -- Player #2: A polytope constructed from facets directly.
  (A, b) = facets P1
  A = lift(A, ZZ)
  b = lift(b, ZZ)
  P2 = polyhedronFromHData(A, b)
  P1 == P2
  assert(vertices P2 === vertices P1)
  LP2 = entries transpose matrix{latticePoints P2}
  LP1 = entries transpose matrix{latticePoints P1}
  assert(LP2 === LP1)

  -- Now let's check that the lattice points of P2 match the new code.
  A * matrix LP0
    assert(max flatten entries(A * LP0) == 1)
  LP2 = matrix {latticePoints P1}
  assert(set entries transpose LP0 == set entries transpose LP2)

  -- d=4 dilation
  LP2d4 = matrix{latticePoints(4*P1)}
  LP0d4 = boundedLatticePoints(A, 4*rhs, Bound => 20, Limit => 10000)
  numcols LP0d4
  assert(set entries transpose LP2d4 === set entries transpose LP0d4)
  assert(numcols LP0d4 == 611)
///

-*
restart
needsPackage "StringTorics"
*-
TEST ///
  ks = KSEntry "4 9  M:11 9 N:58 14 H:50,6 [88] id:300
    1   -1    0    0    0    1    0   -1    2
    0    0    1    0   -1   -1    0    0   -1
    0    0    0    1   -1    0    0    1   -1
    0    0    0    0    0    0    1   -1    0
    "

  Q = reflexivePolytope ks
  P1 = polytope Q
  (first facets P1) * (vertices P1)
    assert(max flatten entries oo == 1)
  assert(numrows first facets P1 == 9)
  elapsedTime latticePoints Q

  -- Construct Ax <= b directly from P1.
  (A, b) = facets P1
  A = lift(A, ZZ)
  b = lift(b, ZZ)
  P2 = polyhedronFromHData(A, b)
  assert(P1 == P2)
  assert(vertices P2 === vertices P1)
  LP2 = entries transpose matrix{latticePoints P2}
  LP1 = entries transpose matrix{latticePoints P1}
  assert(LP2 === LP1)

  LP0 = boundedLatticePoints(A, b, Bound => 10, Limit => 10000)
  assert(max flatten entries(A * LP0) == 1)
  LP2 = matrix {latticePoints P1}
  assert(set entries transpose LP0 == set entries transpose LP2)

  elapsedTime LP1d = boundedLatticePoints(A, 5*b, Bound => 31*3, Limit => 100000);
  assert(numcols LP1d === 10026) -- the number is what the following line gives.
  -- elapsedTime matrix{latticePoints(5*P1)}; -- 11.7 sec too long! (2 June 2026).
///

-*
restart
needsPackage "StringTorics"
*-
TEST /// -- of rawFourierMotzkin, vs FourierMotzkin package.
  debug Core -- for rawFourierMotzkin, raw.
  needsPackage "FourierMotzkin"

  ks = KSEntry "4 9  M:11 9 N:58 14 H:50,6 [88] id:300
    1   -1    0    0    0    1    0   -1    2
    0    0    1    0   -1   -1    0    0   -1
    0    0    0    1   -1    0    0    1   -1
    0    0    0    0    0    0    1   -1    0
    "
  A = matrix ks
  A1 = A || matrix{{9:1}}

  B1 = rays dualCone posHull A1 -- dualCone is B1 s.t. (transpose A1) * A1 >= 0
  assert all(flatten entries ((transpose A1) * B1), a -> a >= 0)

  B1' = - first fourierMotzkin A1 -- columns of B1' are the same as columns of B1 (up to order)
  assert(set entries transpose B1 === set entries transpose B1')

  B1'' = - transpose map(ZZ, rawFourierMotzkin raw transpose A1)
  assert(set entries transpose B1 === set entries transpose B1'')

  A = matrix ks
  A1 = (A | A) || matrix{{18:1}}

  B1 = rays dualCone posHull A1 -- dualCone is B1 s.t. (transpose A1) * A1 >= 0
  (transpose A1) * B1

  B1' = - first fourierMotzkin A1 -- columns of B1' are the same as columns of B1 (up to order)
  assert(set entries transpose B1 === set entries transpose B1')

  B1'' = - transpose map(ZZ, rawFourierMotzkin raw transpose A1)
  assert(set entries transpose B1 === set entries transpose B1'')
///

end--

