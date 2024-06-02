---------------------------------------
-- CYPolytope ---------------------
---------------------------------------
-- This type can be written to disk, and tries to retain computations computed already.
-- It does not retain Polyhedron objects, but hopefully it recreates these quickly.
-- 

CYPolytopeFields = {
    "rays" => {value, toString, List}
    }

-- These are the cache fields that we write to a string via 'dump'
CYPolytopeCache = {
    -- these fields do not need to exist.
    "face dimensions" => {value, toString, List},
    "id" => {value, toString, ZZ},
    "favorable" => {value, toString, Boolean},
    "h11" => {value, toString, ZZ},
    "h21" => {value, toString, ZZ},
    "basis indices" => {value, toString, List},
    "glsm" => {value, toString, List},
    "annotated faces" => {value, toString, List},
    "automorphisms" => {value, toString, List},
    "autPermutations" => {value, toString, List},
    "triangulations" => {value, toString, List}
    }

cyPolytope = method(Options => {ID => null, InteriorFacets => false})

-- This is the main creation function.  Other functions call this.
-- Goal: this function does NOT change vertices list
-- TODO: currently it is NOT this!  

-- vertices: A list of the integer coordinates (also a list)
cyPolytope List := CYPolytope => opts -> vertices -> (
    cyPolytope(transpose matrix vertices, opts)
    )

cyPolytope Polyhedron := opts -> P2 -> (
    topdim := if opts.InteriorFacets then dim P2 - 1 else dim P2 - 2;
    LP := latticePointList P2;
    LPdim := for lp in LP list dim(P2, minimalFace(P2, lp));
    -- now remove the ones that are in facets (or the origin):
    LP = for i from 0 to #LP-1 list if LPdim#i <= topdim then LP#i else continue;
    LPdim = for i from 0 to #LP-1 list if LPdim#i <= topdim then LPdim#i else continue;
    cyData := new CYPolytope from {
        symbol cache => new CacheTable,
        "rays" => LP,
        };
    cyData.cache#"face dimensions" = LPdim;
    if opts.ID =!= null then cyData.cache#"id" = opts.ID;
    cyData
    )

-- This version contains ALL lattice points
-- cyPolytope Polyhedron := opts -> P2 -> (    
--     LP := latticePointList P2;
--     LPdim := for lp in LP list dim(P2, minimalFace(P2, lp));
--     -- now remove the ones that are in facets (or the origin):
--     LP = for i from 0 to #LP-1 list if LPdim#i <= 3 then LP#i else continue;
--     LPdim = for i from 0 to #LP-1 list if LPdim#i <= 3 then LPdim#i else continue;
--     cyData := new CYPolytope from {
--         symbol cache => new CacheTable,
--         "rays" => LP,
--         };
--     cyData.cache#"face dimensions" = LPdim;
--     if opts.ID =!= null then cyData.cache#"id" = opts.ID;
--     cyData
--     )

-- vertices: Matrix whose columns are the vertices of the reflexive polytope.
cyPolytope Matrix := CYPolytope => opts -> vertices -> (
    P2 := convexHull vertices;
    cyPolytope(P2, opts)
    )
cyPolytope KSEntry := CYPolytope => opts -> tope -> (
    -- KSEntry is a Kreuzer-Skarke polytope entry, returned from
    --   ReflexivePolytopesDB functions.
    P1 := convexHull matrix tope;
    P2 := polar P1;
    cyPolytope(P2, opts)
    )

cyPolytope String := CYPolytope => opts -> str -> (
    L := lines str;
    if L#0 != "CYPolytopeData" then error "string is not in proper format";
    fields := hashTable for i from 1 to #L-1 list getKeyPair L#i;
    -- First get the main elements (these are required!):
    required := for field in CYPolytopeFields list (
        k := field#0;
        readFcn := field#1#0;
        if fields#?k then k => readFcn fields#k else error("expected key "|k)
        );
    cyData := new CYPolytope from prepend(symbol cache => new CacheTable, required);
    -- now read in the cache values (including "id" value, if any)
    for field in CYPolytopeCache do (
        k := field#0;
        readFcn := field#1#0;
        if fields#?k then cyData.cache#k = readFcn fields#k;
        );
    if opts.ID =!= null then cyData.cache#"id" = opts.ID; -- just for compatibility with other constructors...
    cyData
    )

-- todo: translation function: {1, 2, 3, 6} ==> "1 2 3 6" (and viceversa)
-- todo: translation function: {{1,3},{4,7},{6,7,8}} ==> "1 3;4 7;6 7 8;" or "1 3;4 7;6 7 8" (white space is not relevant after or before a ;)
-- Format
-- CYPolytope
--   rays: 1 0 0; 1 0 -1; 1 1 1
--   face dimensions: 0 0 0
--   id: 12
--   favorable: true
--   h11: 5
--   h21: 20
--   basis indices: 0 1 2 3
--   glsm: 1 1 1; 1 2 3

-- Then need to be able to set fields
--
-- Need a isWellFormed function.  Checks that the correct fields are
-- present, and the lengths of the various integer vectors and lists
-- are compatible.

dump CYPolytope := String => {} >> opts -> (Q) -> (
    s1 := "CYPolytopeData\n";
    strs := for field in CYPolytopeFields list (
        k := field#0;
        writerFunction := field#1#1;
        if not Q#?k then error("expected key: "|k#0);
        "  " | k | ":" | writerFunction(Q#k) | "\n"
        );
    strs2 := for field in CYPolytopeCache list (
        k := field#0;
        writerFunction := field#1#1;
        if not Q.cache#?k then continue;
        "  " | k | ":" | writerFunction(Q.cache#k) | "\n"
        );
    strs = join({s1}, strs, strs2);
    concatenate strs
    )

getKeyPair = method()
getKeyPair String := Sequence => str -> (
    str1 := replace("^ *", "", str);
    result := separate(" *: *", str1); -- separate at colon, ignoring white space around colon.
    if #result != 2 then error("expected a key and a value for "|str);
    toSequence result
    )

findTwoFaceInteriorDivisors = method()
findTwoFaceInteriorDivisors CYPolytope := List => Q -> (
    -- returns a list of:
    -- {i, {g, ind}}
    -- {i:nonfavorable divisor index, {g:genus of 2-face, ind:index of 2-face in annotatedFaces Q}}
    A := annotatedFaces Q;
    on1skeleton := set sort unique flatten for x in A list if x#0 <= 1 then x#2 else continue;
    A2 := positions(A, x -> x#0 == 2 and x#3 > 0 and x#4 > 0);
    flatten for a in A2 list (
        thisface := A#a; -- note thisface#2 is the list of all lattice point indices on this face.
                          --      thisface#4 is the genus of this face.
        nonfavs := sort toList(set thisface#2 - on1skeleton);
        for x in nonfavs list {x, {thisface#4, a}}
        )
    )

-- choosing basis indices: if any non-favorable rays, try to choose them!
-- then we can simply replace that generator with the g+1 that sum to it.

findSuitableSet = (setstotry, Z) -> (
    for g in setstotry do if abs det(Z_g) == 1 then return g;
    null
    )

-- This has been subsumed below?
-- computeBasis = method()
-- computeBasis CYPolytope := List => Q -> (
--     -- first find 2-face interiors with g>0.
--     -- our plan is to find a basis including these, so that we can 
--     -- easily just replace them with the divisors of the form (i, j), 0 <= j <= g, for i non-favorable.
--     nonfavsList := findTwoFaceInteriorDivisors Q;
--     nonfavs := for f in nonfavsList list f#0; -- list of indices of non-favorable divisors.
--     M := transpose matrix rays Q;
--     Z := transpose LLL syz M;
--     rest := sort toList(set(0..numcols Z-1) - set nonfavs);
--     setstotry := for f in subsets(rest, numrows Z - #nonfavs) list (f | nonfavs); -- really want to do these 1 by 1...?
--     good := findSuitableSet(setstotry, Z);
--     if good === null then error "rats: cannot find basis set including all the non-favorables";
--     H := hashTable nonfavsList;
--     flatten for i in good list if not H#?i then i else (
--         g := H#i#0; -- genus of the 2-face
--         for j from 0 to g list (i,j)
--         )
--     )

cySetGLSM = method()
-- Delete this version?
-- cySetGLSM CYPolytope := (cyData) -> (
--     if cyData.cache#?"glsm" then return;
--     mLP := transpose matrix cyData#"rays";
--     D := transpose syz mLP;
--     p := findFirstUnitVectors D; -- TODO: p,q computation can be slow!
--     q := findInvertibleSubmatrix(D, p);
--     if q === null then error ("oops, can't find a good GLSM matrix"); -- hasn't happened yet. HAS NOW!!
--     GLSM := (D_q)^-1 * D;
--     cyData.cache#"glsm" = entries transpose GLSM;
--     cyData.cache#"basis indices" = q
--     )
cySetGLSM CYPolytope := Q -> (
    if Q.cache#?"glsm" then return;
    mLP := transpose matrix rays Q;
    D := transpose syz mLP; -- use LLL?
    -- D := transpose LLL syz M; -- which line should we use?
    nonfavsList := findTwoFaceInteriorDivisors Q;
    -- TODO: should nonfavsList be stashed into Q?
    nonfavs := for f in nonfavsList list f#0; -- list of indices of non-favorable divisors.
    rest := sort toList(set(0..numcols D-1) - set nonfavs);
    setstotry := for f in subsets(rest, numrows D - #nonfavs) list (f | nonfavs); -- really want to do these 1 by 1...?
    good := findSuitableSet(setstotry, D);
    if good === null then error "rats: cannot find basis set including all the non-favorables";
    H := hashTable nonfavsList;
    basind := flatten for i in good list if not H#?i then i else (
        g := H#i#0; -- genus of the 2-face
        for j from 0 to g list (i,j)
        );
    GLSM := (D_good)^-1 * D;
    Q.cache#"basis indices" = basind;
    Q.cache#"toric basis indices" = good;
    Q.cache#"glsm" = entries transpose GLSM;
    )

cySetH11H21 = cyData -> (
    -- this version is only for CY 3-fold hypersurfaces...
    -- P:ReflexivePolytope
    -- P := polytope cyData;
    A := annotatedFaces cyData; -- polytope on N side.
    A0 := for x in A list if x#0 == 0 then drop(x,1) else continue; -- annotatedFaces(0, P);
    A1 := for x in A list if x#0 == 1 then drop(x,1) else continue; -- annotatedFaces(1, P);
    A2 := for x in A list if x#0 == 2 then drop(x,1) else continue; -- annotatedFaces(2, P);
    A3 := for x in A list if x#0 == 3 then drop(x,1) else continue; -- annotatedFaces(3, P);
    npM := A/(x -> x#4)//sum + 1;
    npN := A/(x -> x#3)//sum; -- origin is included in the dim 4 face.
    -- points in facets (on M side) -- this is part of h21
    -- points in facets (on N side) -- this is part of h11
    facetInteriorsM := A0/(v -> v#3)//sum;
    facetInteriorsN := A3/(v -> v#2)//sum;
    -- points interior to 2-faces (times their genus) (on M-side)
    -- points interior to 2-faces (times their genus) (on N-side)
    twoFacesM := A1/(v -> v#2 * v#3)//sum;
    twoFacesN := A2/(v -> v#2 * v#3)//sum;
    -- now set the h11, h21.
    h11 := npN - 5 - facetInteriorsN + twoFacesN;
    h21 := npM - 5 - facetInteriorsM + twoFacesM;
    cyData.cache#"h11" = h11;
    cyData.cache#"h21" = h21;
    cyData.cache#"favorable" = (twoFacesN == 0);
    (h11, h21)
    )

rays CYPolytope := List => {} >> opts -> cyData -> (
    cyData#"rays"
    )

dim CYPolytope := List => cyData -> dim polytope(cyData, "N")
degrees CYPolytope := List => cyData -> (
    if not cyData.cache#?"glsm" then cySetGLSM cyData;
    cyData.cache#"glsm"
    )
basisIndices = method()
basisIndices CYPolytope := List => cyData -> (
    if not cyData.cache#?"basis indices" then cySetGLSM cyData;
    cyData.cache#"basis indices"
    )
-- h11OfCY CYPolytope := ZZ => cyData -> (
--     if not cyData.cache#?"h11" then cySetH11H21 cyData;
--     cyData.cache#"h11"
--     )
-- h21OfCY CYPolytope := ZZ => cyData -> (
--     if not cyData.cache#?"h21" then cySetH11H21 cyData;
--     cyData.cache#"h21"
--     )
isFavorable CYPolytope := Boolean => cyData -> (
    if not cyData.cache#?"favorable" then cySetH11H21 cyData;
    cyData.cache#"favorable"
    )
annotatedFaces CYPolytope := cyData -> (
    if not cyData.cache#?"annotated faces" then
      cyData.cache#"annotated faces" = annotatedFaces polytope(cyData, "N");
    cyData.cache#"annotated faces"
    )
polytope(CYPolytope, String) := Polyhedron => (cyData, which) -> (
    if which === "N" then (
        if not cyData.cache#?"N polytope" then (
            LP := cyData#"rays";
            LPdim := cyData.cache#"face dimensions";
            verts := for i from 0 to #LP - 1 list if LPdim#0 == 0 then LP#i else continue;
            cyData.cache#"N polytope" = convexHull transpose matrix verts
            );
        cyData.cache#"N polytope"
        )
    else if which === "M" then (
        if not cyData.cache#?"M polytope" then (
            cyData.cache#"M polytope" = polar polytope(cyData, "N");
            );
        cyData.cache#"M polytope"
        )
    else
      error "expected second argument to be either \"M\" or \"N\""
    )
polytope CYPolytope := Polyhedron => cyData -> polytope(cyData, "N")

polar CYPolytope := cyData -> cyPolytope polytope(cyData, "M")

findAllFRSTs CYPolytope := List => cyData -> (
    if not cyData.cache#?"triangulations" then
        cyData.cache#"triangulations" = (findAllFRSTs(transpose matrix rays cyData))/last;
    cyData.cache#"triangulations"
    )

normalizeByAutomorphisms = method()
normalizeByAutomorphisms(List, List) := (gPerms, T) -> (
    -- gPerms should be a list of permutations of 0..#rays-1, for a CYPolytoe Q.
    -- T should be a list of list of integer indices into the rays of Q.
    first sort for g in gPerms list (
        sort for t in T list sort g_t
        )
    )

findAllCYs = method(Options => {Ring => null, NTFE => true, Automorphisms => true}) -- opts.Ring: ZZ[h11 variables].
findAllCYs CYPolytope := List => opts -> Q -> (
    Ts := findAllFRSTs Q;
    RZ := if opts#Ring === null then (
        a := getSymbol "a";
        h11 := hh^(1,1) Q;
        ZZ[a_1 .. a_h11]
        )
    else (
        opts#Ring
        );
    Xs := for i from 0 to #Ts - 1 list cyData(Q, Ts#i, Ring => RZ); -- we set the ID below.
    -- If NTFE and UseAutomorphisms:
    gPerms := if opts.Automorphisms then 
                 automorphismsAsPermutations Q
              else 
                 {splice{0..#rays Q - 1}}; -- only the identity permutation
    -- f is the function we use to partition the Xs.
    f := if opts.NTFE then 
             (X -> normalizeByAutomorphisms(gPerms, restrictTriangulation(2, X)))
         else 
             (X -> normalizeByAutomorphisms(gPerms, max X));
    H := partition(f, Xs);
    count := 0;
    Xs = for k in keys H list (
        X := H#k#0; -- take the first one
        X.cache#"id" = count;
        count = count+1;
        X);
    Xs
    )

hh(Sequence, CYPolytope) := (pq, Q) -> (
    cySetH11H21 Q;
    (p,q) := pq;
    if p > q then (p, q) = (q, p);
    if p == 0 then (
        if q == 3 or q == 0 then 1 else 0
        )
    else if p == 1 then (
        if q == 1 then Q.cache#"h11"
        else if q == 2 then Q.cache#"h21"
        else 0
        )
    else if p == 2 then (
        if q == 2 then Q.cache#"h11" else 0
        )
    else if p == 3 then (
        if q == 3 then 1
        else 0
        )
    )


isomorphisms(CYPolytope, CYPolytope) := (P, Q) -> (
    isomorphisms(polytope P, polytope Q, annotatedFaces P, annotatedFaces Q)
    )

automorphisms CYPolytope := Q -> (
    if not Q.cache#?"automorphisms" then Q.cache#"automorphisms" = (
        P := polytope(Q, "N");
        auts := isomorphisms(Q, Q);
        sort for x in auts list entries x
        );
    Q.cache#"automorphisms"
    )

automorphismsAsPermutations = method()
automorphismsAsPermutations CYPolytope := Q -> (
    if not Q.cache#?"autPermutations" then Q.cache#"autPermutations" = (
        G := automorphisms Q;
        raysQ := rays Q;
        raysMatrices := for v in rays Q list transpose matrix {v};
        raysHash := hashTable for i from 0 to #raysMatrices - 1 list raysMatrices#i => i;
        for g in G list (
            m := matrix g;
            for v in raysMatrices list raysHash#(m * v)
            )
        );
    Q.cache#"autPermutations"
    )

    -- nrows := numrows vertexMatrix P;
    -- if nrows != dim P or nrows != dim Q or nrows != numrows vertexMatrix Q
    -- then error "expected polytoeps to be full dimensional and same dimension";
    
    -- -- Step 1. Find a facet of P with the smallest size.
    -- annCYP := annotatedFaces CYP;
    -- facetsP := for f in annCYP list if f#0 != nrows-1 then continue else f#1;
    -- minsizeP := facetsP/length//min;
    -- facetsMinsizeP := select(facetsP, f -> #f === minsizeP);
    -- facetA := first facetsMinsizeP;

    -- -- Step 2. Find all facets of Q with this smallest size minsizeP, or return {}.
    -- annCYQ := annotatedFaces CYQ;
    -- facetsQ := for f in annCYQ list if f#0 != nrows-1 then continue else f#1;
    -- minsizeQ := facetsQ/length//min;
    -- if minsizeQ =!= minsizeP then (
    --     error "debug me";
    --     return {};
    --     );
    -- facetsMinsizeQ := select(facetsQ, f -> #f === minsizeP);

    -- -- Now find a subset of nrows elements if facetA which are full dimensional
    -- -- TODO: don't assume facetA is norws!
    -- if #facetA > nrows then (
    --     -- we need to take a subset of these of size nrows that have full rank.
    --     -- we then call these facetA again.  We don't actually need facetA again,
    --     -- the only thing we use is Ainv.
    --     C := ((vertexMatrix P)_facetA) ** QQ;
    --     facetA = facetA _ (columnRankProfile mutableMatrix C);
    --     if #facetA != nrows then error "my logic is missing a case";
    --     );
    -- A := (vertexMatrix P)_facetA;
    -- Ainv := (A ** QQ)^-1;

    -- -- now we loop through all possible maps from facetA to other facets,
    -- -- and if it gives an integer matrix, we add it to the list.
    -- vertsP := (vertexList P)/(v -> transpose matrix {v});
    -- vertsQ := (vertexList Q)/(v -> transpose matrix {v});
    -- hashQ := hashTable for i from 0 to #vertsQ-1 list vertsQ#i => i;
    -- elapsedTime isos := flatten for f in facetsMinsizeQ list (
    --     for perm in partialPermutations(f, nrows) list (
    --         B := (vertexMatrix Q)_perm;
    --         M := B * Ainv;
    --         try (M = lift(M, ZZ)) else continue;
    --         if all(vertsP, v -> hashQ#?(M * v)) then M else continue
    --         )
    --     );
    -- return isos;
    -- )                                                       
