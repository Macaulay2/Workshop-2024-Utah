-- TODO: I think changing the name to ReflexiveData is better?
--   data will contain: an ordered list of ALL lattice points.
--   I think I would like to make the restriction that the points
--   occur in some order (e.g. ordered by face dimension).

---------------------------------------
-- CYPolytope ---------------------
---------------------------------------
-- This type can be written to disk, and tries to retain computations computed already.
-- It does not retain Polyhedron objects, but hopefully it recreates these quickly.
--

-- A CYPolytope is fundamentally, a list of *all* the non-zero lattice points of
-- a reflexive polytope, together with cached information about the polytope,
-- which is useful for constructing Calabi-Yau hypersurfaces in simplicial
-- toric varieties arising from this polytope.
-- The information which is cached is designed to be easy to read and write to a disk
-- or database.  It is possible to supply some of this cached data, including
-- glsm degree matrix and annotated faces, but generally, the package will compute these
-- for you.
-- The above is not currently correct: (1) it is currently a list of all lattice points
--   not in facets, (2) you cannot set the rays yourself, otherwise the indices are all wrong
--   (in e.g.) annotated faces.
--   (3) you cannot easily set "glsm" and "basis indices".
-- Note that "basis indices", for nonfavorable, includes not just indices,
-- but items like (3,0), (3,1), ...(3,g), where 3 is the index of a lattice point
-- in the interior of a 2-face.
-- This will all change with the introduction of cyPolytopeFromRays, and a function
-- which makes it easy to collect the rays from the vertices of a polytope.

CYPolytopeFields = {
    "rays" => {value, toString, List}
    }

-- These are the cache fields that we write to a string via 'dump'
CYPolytopeCache = {
    -- these fields may or may not exist in a specific CYPolytope object.
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

cyPolytope = method(Options => {ID => null, InteriorFacets => false}) -- TODO: remove InteriorFaces.  That should be in construction of CY's.

-- This is the main creation function.  Other functions call this.
-- Goal: this function does NOT change vertices list
-- TODO: currently it is NOT this!  

-- vertices: A list of the integer coordinates (also a list)
cyPolytopeFromRays = method(Options => {ID => null, Degrees => null, BasisIndices => null})
cyPolytopeFromRays List := CYPolytope => opts -> allrays -> (
    -- we assume that the convex hull of all the rays is a reflexive polytope.
    -- we also assume that the origin is not an array (and so in any
    -- case, all rays are on the reflexive polytope.
    error "cyPolytopeFromRays function is not yet implemented";
    Q := new CYPolytope from {
        symbol cache => new CacheTable,
        "rays" => allrays
        };
    if opts.ID =!= null then Q.cache#"id" = opts.ID;
    Q
    )

cyPolytopeFromRays(List, List) := opts -> (verts, latticepoints) -> (
    error "cyPolytopeFromRays function is not yet implemented";
    -- WARNING: the order of lattice points given is used, the indices of verts given
    --  are not kept!
    -- make P2 from verts
    -- get lattice points
    -- get face dimensions for each lattice point (in that order!)
    -- maybe: stash the change of basis from latticePointList P2 to latticepoints?
    )

-- cyPolytope List := opts -> verts -> (
--     -- make P2 from verts
--     -- get lattice points
--     -- get face dimensions for each lattice point (in that order!)
--     )

-- When do we need:
-- M polytope
-- N polytope
-- lattice points
--
-- should use:
--  rays Q: all lattice points other than the origin.
--  vertices Q: returns the rays which are vertices of the polytope.


cyPolytope List := CYPolytope => opts -> vertices -> (
    cyPolytope(transpose matrix vertices, opts)
    )
cyPolytope Matrix := CYPolytope => opts -> vertices -> (
    P2 := convexHull vertices;
    cyPolytope(P2, opts)
    )
cyPolytopeWithGivenLatticePoints = method(Options => options cyPolytope)
cyPolytopeWithGivenLatticePoints(List, List) := opts ->  (latticePoints, faceDimensions) -> (
    error "cyPolytopeWithGivenLatticePoints function is not yet implemented";
    verts := for i from 0 to #latticePoints - 1 list if faceDimensions#i == 0 then latticePoints#i else continue;
    P2 := convexHull transpose matrix verts;
    (LP, LPdim) := latticePointsAndDimensions P2;
    -- now check that the give lattice points and face dimensions match:
    if #latticePoints =!= #LP then error("expected "|#LP|" lattice points");
    if #faceDimensions =!= #latticePoints then error "expected both arguments to have the same length";
    LPset1 := for i from 0 to #latticePoints - 1 list latticePoints#i => faceDimensions#i;
    LPset2 := for i from 0 to #LP - 1 list LP#i => LPdim#i;
    if sort LPset1 =!= sort LPset2 then error "expected a list of all lattice points";
    Q := new CYPolytope from {
        symbol cache => new CacheTable,
        "rays" => latticePoints,
        };
    Q.cache#"face dimensions" = faceDimensions;
    Q.cache#"N polytope" = P2;
    if opts.ID =!= null then Q.cache#"id" = opts.ID;
    Q
    )

    
-- TODO: remove once we are using InteriorFacets when constructing CY's
-- in construction: move InteriorFacets checks to elsewhere
-- cyPolytope Polyhedron := opts -> P2 -> (
--     topdim := if opts.InteriorFacets then dim P2 - 1 else dim P2 - 2;
--     LP := latticePointList P2;
--     LPdim := for lp in LP list dim(P2, minimalFace(P2, lp));
--     -- now remove the ones that are in facets (or the origin):
--     LP = for i from 0 to #LP-1 list if LPdim#i <= topdim then LP#i else continue;
--     LPdim = for i from 0 to #LP-1 list if LPdim#i <= topdim then LPdim#i else continue;
--     cyData := new CYPolytope from {
--         symbol cache => new CacheTable,
--         "rays" => LP,
--         };
--     cyData.cache#"face dimensions" = LPdim;
--     if opts.ID =!= null then cyData.cache#"id" = opts.ID;
--     cyData
--     )

-- in construction: move InteriorFacets checks to elsewhere?
cyPolytope Polyhedron := opts -> P2 -> (
    topdim := if opts.InteriorFacets then dim P2 - 1 else dim P2 - 2;
    (LP, LPdim) := latticePointsAndDimensions P2;
    LP = for i from 0 to #LP-1 list if LPdim#i <= topdim then LP#i else continue;
    LPdim = for i from 0 to #LP-1 list if LPdim#i <= topdim then LPdim#i else continue;
    Q := new CYPolytope from {
        symbol cache => new CacheTable,
        "rays" => LP,
        };
    Q.cache#"face dimensions" = LPdim;
    Q.cache#"N polytope" = P2;
    if opts.ID =!= null then Q.cache#"id" = opts.ID;
    Q
    )

faceDimensions = method()
faceDimensions CYPolytope := Q -> (
    if not Q.cache#?"face dimensions" then (
        P2 := polytope(Q, "N");
        pts := Q#"rays";
        -- LP := latticePointList P2;
        -- if sort pts != sort LP then (
        --     error "expected the given lattice points to be the same (except for their order) as the provided points";
        --     -- TODO: create the reorder permutation(s).
        --     );
        LPdim := for lp in pts list dim(P2, minimalFace(P2, lp));
        Q.cache#"face dimensions" = LPdim;
        );
    Q.cache#"face dimensions"
    )

vertices CYPolytope := Q -> (
    raysQ := rays Q;
    dims := faceDimensions Q;
    for i from 0 to #rays Q - 1 list if dims#i > 0 then continue else raysQ#i
    )

latticePoints CYPolytope := Q -> latticePointList polytope(Q, "N")

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
    Q := new CYPolytope from prepend(symbol cache => new CacheTable, required);
    -- now read in the cache values (including "id" value, if any)
    for field in CYPolytopeCache do (
        k := field#0;
        readFcn := field#1#0;
        if fields#?k then Q.cache#k = readFcn fields#k;
        );
    if opts.ID =!= null then Q.cache#"id" = opts.ID; -- just for compatibility with other constructors...
    Q
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
    dimQ := dim Q;
    rayIndices := for i from 0 to #rays Q - 1 list if (faceDimensions Q)#i <= dimQ-2 then i else continue;
    raysToKeep := for i in rayIndices list (rays Q)#i;
    mLP := transpose matrix raysToKeep;
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

rays CYPolytope := List => {} >> opts -> Q -> Q#"rays"

dim CYPolytope := ZZ => Q -> dim polytope(Q, "N")

degrees CYPolytope := List => Q -> (
    if not Q.cache#?"glsm" then cySetGLSM Q;
    Q.cache#"glsm"
    )

basisIndices = method()
basisIndices CYPolytope := List => Q -> (
    if not Q.cache#?"basis indices" then cySetGLSM Q;
    Q.cache#"basis indices"
    )

isFavorable CYPolytope := Boolean => Q -> (
    if not Q.cache#?"favorable" then cySetH11H21 Q;
    Q.cache#"favorable"
    )

-- faceDimensions CYPolytope := List => Q -> (
--     if not Q.cache#?"favorable" then cySetFaceDimensions Q;
--     Q.cache#"favorable"
--     )

annotatedFaces CYPolytope := Q -> (
    if not Q.cache#?"annotated faces" then (
        P2 := polytope(Q, "N");
        result := annotatedFaces P2;
        LPlist := latticePointList P2;
        if take(LPlist, #(rays Q)) =!= rays Q then (
            --if sort LPlist =!= sort rays Q then error "internal error in annotatedFaces: Q is not well defined";
            raysQ := rays Q;
            rayHash := hashTable for i from 0 to #raysQ - 1 list raysQ#i => i;
            lpHash := hashTable for i from 0 to #LPlist - 1 list LPlist#i => i;
            lp2rays := for i from 0 to #LPlist-1 list rayHash#(LPlist#i); -- a permutation
            result = sort for f in result list (
                -- each entry is of the form {dim of face, indices of vertices in face, indices of all lps in face, genus, genus}
                -- only items #1, #2 need to be recomputed.
                if #f =!= 5 then error "my logic is wrong";
                {f#0, sort for a in f#1 list lp2rays#a, sort for a in f#2 list lp2rays#a, f#3, f#4}
                );
            );
        Q.cache#"annotated faces" = result;
        );
    Q.cache#"annotated faces"
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

-- This function returns true if the given triangulation of the vector configuration
-- given by `rays Q` is also a (star) triangulation of the polytope.  This is the
-- case when every maximal cone in the triangulation is contained in a facet of the
-- polytope.
isTriangulationOfPolytope = method()
isTriangulationOfPolytope(CYPolytope, List) := (Q, T) -> (
    facetsQ := for x in annotatedFaces Q list if x#0 =!= dim Q - 1 then continue else set x#2;
    for t in T do (
        if any(facetsQ, f -> isSubset(t, f)) then continue else return false;
        );
    true
    )
isTriangulationOfPolytope(CYPolytope, Triangulation) := (Q, T) -> isTriangulationOfPolytope(Q, max T)

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
