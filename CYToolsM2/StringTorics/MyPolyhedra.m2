-- also defined: 
--  (1) dim(P,f)
--  (2) genus(P,f)

-- functions from Polyhedra which we use here:
--   latticePoints P -- a different order than 'latticePoints P'
--   vertices P -- a different order than 'vertexList P' !
--   faces(d,P) -- only in faceDimensionHash
--   polar P -- used often
--   dim P -- used often

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
if Polyhedra#Options#Version == "1.3" then ( -- version of Polyhedra in M2 <= 1.9.2
  faceDimensionHash Polyhedron := (cacheValue symbol TCIFaceDimensionHash) ((P) -> (
    L := vertexList P;
    vertexHash := hashTable for i from 0 to #L - 1 list (L#i => i);
    hashTable flatten for i from 0 to dim P list for f in faces(dim P-i,P) list (
        -- USING INFO FROM POLYHEDRA
        verts := vertices f;
        --if liftable(verts,ZZ) then verts = lift(verts,ZZ);
        verts = try lift(verts,ZZ) else verts;
        (sort for v in entries transpose verts list vertexHash#v) => i
        )
    ))
) else (
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
)

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

-- -- private function for `isomorphisms`
-- findCombinatorialData = (aP, i) -> (
--     -- aP: List, coming from annotated faces.
--     -- i: integer index: for a given vertex.
--     -- returns: list, of
--     --  {genus, edges: genus => number, 2faces: {#vertices, genus} => number}, 3faces: {#vertices, genus} => ZZ
--     -- these are counts for all faces containing i
--     -- or, maybe one hash table, {#vertices, genus} => count
--     -- #vertices can be 1,2,3.
--     sort for x in aP list if member(i, x#1) then {x#0, #x#1, x#3, x#4} else continue
--     )

-- -- private function for `isomorphisms`
-- findPossibleMatchings = (matchings) -> (
--     -- matchings: List of (alpha, beta), alpha and beta lists of integers of the same length >= 1.
--     -- returns a list of {i1 => j1, ..., ir => jr}.
--     -- the possible matchings.
--     if #matchings === 0 then error "incorrect logic on my part, apparently";
--     hd := matchings#0;
--     hdmatchings0 := permutations(#hd#0);
--     hdmatchings := for p in hdmatchings0 list for i from 0 to #hd#0-1 list hd#0#i => hd#1#(p#i);
--     if #matchings === 1 then return hdmatchings;
--     tl := drop(matchings, 1);
--     restmatchings := findPossibleMatchings tl;
--     Ps := permutations hd#1;
--     flatten for p in hdmatchings list for q in restmatchings list (
--         sort join(p,q)
--         )
--     )

-- ///
--   findPossibleMatchings{({0,1,2}, {0,1,2})}
--   findPossibleMatchings{({0,3}, {0,1})}
-- ///

-- checkPossibleMatching = (A, perm, verticesP, verticesQ) -> (
--     -- A is a n x n generic matrix over n^2 variables.
--     -- perm is a list {i1 => j1, ...} of indices into verticesP to indices into verticesQ
--     trim sum for ab in perm list (
--         a := transpose matrix{verticesP _ (first ab)};
--         b := transpose matrix{verticesP _ (last ab)};
--         ideal (A * a - b)
--         )
--     )

-- matchings = method()
-- matchings(List, List, ZZ) := (aP, aQ, nvertices) -> (
--     HP := partition(i -> findCombinatorialData(aP, i), toList(0..nvertices - 1));
--     HQ := partition(i -> findCombinatorialData(aQ, i), toList(0..nvertices - 1));
--     if sort keys HP =!= sort keys HQ then (
--         << "note: vertex data does not match" << endl;
--         return {};
--         );
--     if not all(keys HP, k -> #HP#k == #HQ#k) then (
--         << "note: vertex number data does not match" << endl;
--         return {};
--         );
--     matchings := for k in keys HP list (
--         HP#k, HQ#k
--         );
--     matchings
--     )

-- isomorphisms = method()
-- isomorphisms(Polyhedron, Polyhedron) := (P, Q) -> (
--     -- for now, we assume both are full dimensional?
--     vP := vertexList P;
--     vQ := vertexList Q;
--     n := #vP#0; -- TODO: check that all vP, vQ elements have the same length, n == dim P == dim Q
--     if #vP =!= #vQ then return {};
--     -- Step 1: get numerical invariants for each vertex.
--     aP := annotatedFaces P;
--     aQ := annotatedFaces Q;
--     HP := partition(i -> findCombinatorialData(aP, i), toList(0..#vP - 1));
--     HQ := partition(i -> findCombinatorialData(aQ, i), toList(0..#vQ - 1));
--     if sort keys HP =!= sort keys HQ then (
--         << "note: vertex data does not match" << endl;
--         return {};
--         );
--     if not all(keys HP, k -> #HP#k == #HQ#k) then (
--         << "note: vertex number data does not match" << endl;
--         return {};
--         );
--     matchings := for k in keys HP list (
--         HP#k, HQ#k
--         );
--     possibles := findPossibleMatchings matchings;
--     if #possibles > 1000 then (
--         << "#possibles == " << #possibles << endl;
--         return isomorphisms2(P, Q)
--         );
--     --return {possibles, matchings, HP, HQ};
--     t := getSymbol "t";
--     R := QQ[t_(0,0)..t_(n-1,n-1)];
--     A := genericMatrix(R, n, n);
--     As := for p in possibles list (
--         J := checkPossibleMatching(A, p,vP, vQ);
--         if J == 1 then continue; -- not an isomorphism!
--         A0 := A % J;
--         A0 = try lift(A0, ZZ) else null;
--         if A0 === null then continue;
--         (A0, p/last)
--         );
--     As
--     )

-- private function for isomorphisms2
partialPermutations = (elems, num) -> (
    if num == 1 then return elems/(a -> {a});
    flatten for i from 0 to #elems-1 list for p in partialPermutations(drop(elems,{i,i}), num-1)
      list
        prepend(elems#i, p)
    )

-- TODO: isomorphisms and isomorphisms2 should be combined in a smarter way:
-- use the known matches to restrict the possible maps.
-- then use this on only some of the vertices?
-- Vague description because I don't know how best to fix it yet.
-- isomorphisms2 = method()
-- isomorphisms2(Polyhedron, Polyhedron) := (P, Q) -> (
--     -- here we don't bother with matchings.
--     -- instead we first find a set of n vertices which do not lie on a hyperplane.
--     -- and then we compute all possible matrices 
--     VP := vertexMatrix P;
--     VQ := vertexMatrix Q;
--     m := numcols VP;
--     HP := hashTable for i from 0 to m-1 list (vertexList P)#i => i;
--     HQ := hashTable for i from 0 to m-1 list (vertexList Q)#i => i;
--     if m =!= numcols VQ then return {};
--     n := numrows VP;
--     if n =!= numrows VQ then error "expected polytopes in the same space";
--     indepset := for p in subsets(#vertexList P, n) list if det VP_p != 0 then break p;
--     t := getSymbol "t";
--     R := QQ[t_(0,0)..t_(n-1,n-1)];
--     A := genericMatrix(R, n, n);
--     VP0 := VP_indepset;
--     for q in partialPermutations(splice{0..#vertexList Q - 1}, n) list (
--         J := trim ideal(A * VP0 - VQ_q);
--         if J == 1 then continue;
--         A0 := A % J;
--         A0 = try lift(A0, ZZ) else null;
--         if A0 === null then continue;
--         if abs(det A0) != 1 then continue;
--         newverts := entries transpose(A0 * VP);
--         if any(newverts, v -> not HQ#?v) then continue;
--         perm := for v in newverts list HQ#v;
--         (A0, perm)
--         )
--     )

-- Based on code in CYTools.
-- Assumption: P, Q are full rank polyhedra (i.e. dim = #rows of vertices matrix).
-- Idea: find a facet of P with the smallest cardinality.
--       find a subset of these vertices that define a full dimensional set.
--       for each 

isomorphisms = method()
-- isomorphisms(Polyhedron, Polyhedron) := (P, Q) -> (
--     nrows := numrows vertexMatrix P;
--     if nrows != dim P or nrows != dim Q or nrows != numrows vertexMatrix Q
--     then error "expected polytoeps to be full dimensional and same dimension";
    
--     -- Step 1. Find a facet of P with the smallest size.
--     facetsP := (annotatedFaces(3, P))/first;
--     minsizeP := facetsP/length//min;
--     facetsMinsizeP := select(facetsP, f -> #f === minsizeP);
--     facetA := first facetsMinsizeP;

--     -- Step 2. Find all facets of Q with this smallest size minsizeP, or return {}.
--     facetsQ := (annotatedFaces(3, Q))/first;
--     minsizeQ := facetsQ/length//min;
--     if minsizeQ =!= minsizeP then return {};
--     facetsMinsizeQ := select(facetsQ, f -> #f === minsizeP);

--     -- Now find a subset of nrows elements if facetA which are full dimensional
--     if #facetA > nrows then (
--         -- we need to take a subset of these of size nrows that have full rank.
--         -- we then call these facetA again.  We don't actually need facetA again,
--         -- the only thing we use is Ainv.
--         C := ((vertexMatrix P)_facetA) ** QQ;
--         facetA = facetA _ (columnRankProfile mutableMatrix C);
--         if #facetA != nrows then error "my logic is missing a case";
--         );
--     A := (vertexMatrix P)_facetA;
--     Ainv := (A ** QQ)^-1;

--     -- now we loop through all possible maps from facetA to other facets,
--     -- and if it gives an integer matrix, we add it to the list.
--     vertsP := (vertexList P)/(v -> transpose matrix {v});
--     vertsQ := (vertexList Q)/(v -> transpose matrix {v});
--     hashQ := hashTable for i from 0 to #vertsQ-1 list vertsQ#i => i;
--     elapsedTime isos := flatten for f in facetsMinsizeQ list (
--         for perm in partialPermutations(f, nrows) list (
--             B := (vertexMatrix Q)_perm;
--             M := B * Ainv;
--             try (M = lift(M, ZZ)) else continue;
--             if all(vertsP, v -> hashQ#?(M * v)) then M else continue
--             )
--         );
--     isos
--     )

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
  restart
  debug needsPackage "StringTorics"
  topes = kreuzerSkarke 3;
  Q = cyPolytope topes_20
  Q = cyPolytope topes_0
  P = polytope(Q, "N")
  vertexList P
  annotatedFaces P

  isomorphisms3(Q,Q)
-- code I'm working on now
  elapsedTime isomorphisms3(P, P)
  elapsedTime isomorphisms(P,P)

  (Qs, Xs) = readCYDatabase "../Databases/cys-ntfe-h11-3.dbm";
  (Qs, Xs) = readCYDatabase "../Databases/cys-ntfe-h11-4.dbm";
  (Qs, Xs) = readCYDatabase "../Databases/cys-ntfe-h11-5.dbm";

  -- This one is long, I think because the annotated faces for P is not stashed.
  elapsedTime for lab in sort keys Qs list (
      automorphisms Qs#lab
      );

  elapsedTime for lab in sort keys Qs list (
      P := polytope Qs#lab;
      annotatedFaces P);


-- older code
  netList isomorphisms(P,P)
  isomorphisms2(P,P)
  first oo
  netList oo

  for tope in topes list (
      Q := cyPolytope tope;
      P := polytope(Q, "N");
      ans := isomorphisms(P, P);
      if ans == null then 
        << "--- tope: " << label Q << " TOO LARGE FOR NOW" << endl;
      else
        << "--- tope: " << label Q << " #aut=" << #ans << " auts: " << netList ans << endl;
      ans
      );
  
///
end--

