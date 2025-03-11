ReflexivePolytope = new Type of HashTable
ReflexivePolytope.synonym = "reflexive polytope"
ReflexivePolytope.GlobalAssignHook = globalAssignFunction
ReflexivePolytope.GlobalReleaseHook = globalReleaseFunction
expression ReflexivePolytope := X -> if hasAttribute (X, ReverseDictionary) 
    then expression getAttribute (X, ReverseDictionary) else 
    (describe X)#0
describe ReflexivePolytope := X -> Describe (expression ReflexivePolytope) (
    expression X#"vertices")
net ReflexivePolytope := X -> net expression X

-- basic data:
--   list of all of the lattice points.
--   generally these should be sorted in increasing face dimension
--   i.e., vertices come first, the origin comes last.
--   we could relax this, but I think not for the moment...

ReflexivePolytopeFields = {
    "vertices" => {value, toString, List}
    }

-- These are the cache fields that we write to a string via 'dump'
ReflexivePolytopeCache = {
    -- these fields may or may not exist in a specific CYPolytope object.
    "latticePoints" => {value, toString, List},
    "faceDimensions" => {value, toString, List},
    "id" => {value, toString, ZZ},
    "favorable" => {value, toString, Boolean},
    "h11" => {value, toString, ZZ},
    "h21" => {value, toString, ZZ},
    "basisIndices" => {value, toString, List}, -- is this only valid for Batyrev CY's?
    "glsm" => {value, toString, List},
    "annotatedFaces" => {value, toString, List},
    "automorphisms" => {value, toString, List},
    "autPermutations" => {value, toString, List},
    "allFRVTs" => {value, toString, Boolean},
    "ptriangulations" => {value, toString, List},
    "vtriangulations" => {value, toString, List}  -- simplicial fans which are not induced by the polytope
    }

reflexivePolytope = method(Options => {
        ID => null,
        InteriorFacets => false
        }) -- TODO: remove InteriorFaces.  That should be in construction of CY's.

-- `reflexivePolytope` Polyhedron: create a ReflexivePolytope object from a Polyhedra Polyhedron object
reflexivePolytope Polyhedron := ReflexivePolytope => opts -> P2 -> (
    if not isReflexive P2 then error "expected a reflexive polytope";
    verts := vertexList P2;
    -- (LP, LPdim) := latticePointsAndDimensions P2;
    -- verts := for i from 0 to #LP - 1 list if LPdim#i === 0 then LP_i else continue;
    Q := new ReflexivePolytope from {
        symbol cache => new CacheTable,
        "vertices" => verts
        };
    Q.cache#"N polytope" = P2;
    -- Q.cache#"latticePoints" = LP;
    -- Q.cache#"faceDimensions" = LPdim;
    if opts.ID =!= null then Q.cache#"id" = opts.ID;
    Q
    )

dump ReflexivePolytope := String => {} >> opts -> (Q) -> (
    s1 := "ReflexiveData\n";
    strs := for field in ReflexivePolytopeFields list (
        k := field#0;
        writerFunction := field#1#1;
        if not Q#?k then error("expected key: "|k#0);
        "  " | k | ":" | writerFunction(Q#k) | "\n"
        );
    strs2 := for field in ReflexivePolytopeCache list (
        k := field#0;
        writerFunction := field#1#1;
        if not Q.cache#?k then continue;
        "  " | k | ":" | writerFunction(Q.cache#k) | "\n"
        );
    strs = join({s1}, strs, strs2);
    concatenate strs
    )

reflexivePolytope String := ReflexivePolytope => opts -> str -> (
    L := lines str;
    if L#0 != "ReflexiveData" then error ("string is not in proper format, received: "|L#0);
    fields := hashTable for i from 1 to #L-1 list getKeyPair L#i;
    -- First get the main elements (these are required!):
    required := for field in ReflexivePolytopeFields list (
        k := field#0;
        readFcn := field#1#0;
        if fields#?k then k => readFcn fields#k else error("expected key "|k)
        );
    Q1 := new ReflexivePolytope from prepend(symbol cache => new CacheTable, required);
    -- now read in the cache values (including "id" value, if any)
    for field in ReflexivePolytopeCache do (
        k := field#0;
        readFcn := field#1#0;
        if fields#?k then Q1.cache#k = readFcn fields#k;
        );
    if opts.ID =!= null then Q1.cache#"id" = opts.ID; -- just for compatibility with other constructors...
    Q1
    )

reflexivePolytope List := ReflexivePolytope => opts -> vertices -> (
    P2 := convexHull transpose matrix vertices;
    reflexivePolytope(P2, opts)
    )
reflexivePolytope Matrix := ReflexivePolytope => opts -> vertices -> (
    P2 := convexHull vertices;
    reflexivePolytope(P2, opts)
    )

cyPolytope List := CYPolytope => opts -> vertices -> (
    return reflexivePolytope(vertices, opts);
    error "calling cyPolytope List";
    cyPolytope(transpose matrix vertices, opts)
    )
cyPolytope Matrix := CYPolytope => opts -> vertices -> (
    return reflexivePolytope(vertices, opts);
    error "calling cyPolytope Matrix";
    P2 := convexHull vertices;
    cyPolytope(P2, opts)
    )
cyPolytope KSEntry := opts -> tope -> reflexivePolytope(tope, opts)
cyPolytope String := opts -> str -> reflexivePolytope(str, opts)
    


reflexivePolytope KSEntry := ReflexivePolytope => opts -> tope -> (
    -- KSEntry is a Kreuzer-Skarke polytope entry, returned from
    --   ReflexivePolytopesDB functions.
    P1 := convexHull matrix tope;
    P2 := polar P1;
    reflexivePolytope(P2, opts)
    )

vertices ReflexivePolytope := Q -> Q#"vertices"

latticePoints ReflexivePolytope := Q -> (
    if not Q.cache#?"latticePoints" then (
        computeRaysAndDimensions Q;
        );
    Q.cache#"latticePoints"
    )

faceDimensions ReflexivePolytope := Q -> (
    if not Q.cache#?"latticePoints" then (
        computeRaysAndDimensions Q;
        );
    Q.cache#"faceDimensions"
    )

faceDimension = method()
faceDimension(ReflexivePolytope, ZZ) := ZZ => (Q, indx) -> (
    LPdims := faceDimensions Q;
    if indx < 0 or indx >= #LPdims then
        error "index of lattice point out of range";
    LPdims#indx
    )

dim ReflexivePolytope := ZZ => Q -> dim polytope(Q, "N")

computeRaysAndDimensions = method(Options => { InteriorFacets => false })
computeRaysAndDimensions ReflexivePolytope := opts -> Q -> (
    -- sets faceDimensions too, 
    P2 := polytope Q;
    (LP, LPdim) := latticePointsAndDimensions P2;
    Q.cache#"latticePoints" = LP;
    Q.cache#"faceDimensions" = LPdim;
    )

rays ReflexivePolytope := List => { InteriorFacets => false } >> opts -> (Q -> (
    if not Q.cache#?"latticePoints" then (
        computeRaysAndDimensions Q; -- this keeps lattice points sorted in increasing face dimension, *always*
        );
    LP := Q.cache#"latticePoints";
    LPdim := Q.cache#"faceDimensions";
    dimQ := LPdim#-1; -- this is the dimension of the face containing the origin.
    topdim := if opts.InteriorFacets then dimQ - 1 else dimQ - 2;
    topval := position(LPdim, a -> a > topdim);
    if topval === null then error "internal error: ReflexivePolytope is missing an interior vertex!";
    LP_{0..topval-1}
    --for i from 0 to #LP-1 list if LPdim#i <= topdim then LP#i else continue
    ))

polytope(ReflexivePolytope, String) := Polyhedron => (Q, which) -> (
    if not Q.cache#?"N polytope" then (
        Q.cache#"N polytope" = convexHull transpose matrix vertices Q;
        );
    P2 := Q.cache#"N polytope";
    if which === "N" then (
        P2
        )
    else if which === "M" then (
        if not Q.cache#?"M polytope" then (
            Q.cache#"M polytope" = polar P2;
            );
        Q.cache#"M polytope"
        )
    else
      error "expected second argument to be either \"M\" or \"N\""
    )
polytope ReflexivePolytope := Polyhedron => cyData -> polytope(cyData, "N")

polar ReflexivePolytope := ReflexivePolytope => Q -> reflexivePolytope polytope(Q, "M")

annotatedFaces ReflexivePolytope := Q -> (
    if not Q.cache#?"annotatedFaces" then (
        -- note: annotatedFaces uses the list of lattice points of P2, but that is the same as
        -- the list of lattice points of Q
        P2 := polytope(Q, "N");
        Q.cache#"annotatedFaces" = annotatedFaces P2;
        );
    Q.cache#"annotatedFaces"
    )

findTwoFaceInteriorDivisors ReflexivePolytope := List => Q -> (
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

computeGLSM = Q -> (
    -- input: Q:ReflexivePolytope
    -- consequences: fields "glsm", "basisIndices", "toricBasisIndices" are set, if they have not been set yet.
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
    Q.cache#"basisIndices" = basind;
    Q.cache#"toricBasisIndices" = good;
    Q.cache#"glsm" = entries transpose GLSM;
    )

degrees ReflexivePolytope := List => Q -> (
    if not Q.cache#?"glsm" then computeGLSM Q;
    Q.cache#"glsm"
    )

basisIndices ReflexivePolytope := List => Q -> (
    if not Q.cache#?"basisIndices" then computeGLSM Q;
    Q.cache#"basisIndices"
    )

isFavorable ReflexivePolytope := Boolean => Q -> (
    if not Q.cache#?"favorable" then computeH11H21 Q;
    Q.cache#"favorable"
    )

computeH11H21 = Q -> (
    if Q#?"h11" then return;
    -- input: a 4D reflexive polytope
    -- consequences: h11, h21, favorable are all set in Q.cache
    A := annotatedFaces Q; -- polytope on N side.
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
    Q.cache#"h11" = h11;
    Q.cache#"h21" = h21;
    Q.cache#"favorable" = (twoFacesN == 0);
    )

hh(Sequence, ReflexivePolytope) := (pq, Q) -> (
    if dim Q =!= 4 then
        error "Hodge numbers for a ReflexivePolytope are only implemented currently for dim=4";
    computeH11H21 Q;
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

------------------------------------
-- Isomorphisms and automorphisms --
------------------------------------
-- The heavy lifting (such as it is) is in MyPolyhedra.m2
isomorphisms(ReflexivePolytope, ReflexivePolytope) := (P, Q) -> (
    isomorphisms(polytope P, polytope Q, annotatedFaces P, annotatedFaces Q)
    )

computeAutomorphisms = Q -> (
    -- Q is a ReflexivePolytope
    -- this functions sets "automorphisms" and "autPermutations" in Q.cache, if needed.
    -- We are setting both versions at once.
    -- This could be slower, if so, we should change that.
    if Q.cache#?"automorphisms" then return;
    Q.cache#"automorphisms" = (
        P := polytope(Q, "N");
        auts := isomorphisms(Q, Q);
        G := sort for x in auts list entries x;
        G
        );
    Q.cache#"autPermutations" = (
        raysQ := rays Q;
        raysMatrices := for v in rays Q list transpose matrix {v};
        raysHash := hashTable for i from 0 to #raysMatrices - 1 list raysMatrices#i => i;
        for g in G list (
            m := matrix g;
            for v in raysMatrices list raysHash#(m * v)
            )
        ); 
    )

automorphisms ReflexivePolytope := Q -> (
    computeAutomorphisms Q;
    Q.cache#"automorphisms"
    )
automorphismsAsPermutations ReflexivePolytope := Q -> (
    computeAutomorphisms Q;
    Q.cache#"autPermutations"
    )

-------------------------------
-- Triangulations -------------
-- Both point configurations --
-- and vector configurations --
-------------------------------
isTriangulationOfPolytope(ReflexivePolytope, List) := Boolean => (Q, T) -> (
    -- important assumption: T *is* a triangulation of the vector configuration given
    -- by `rays Q`.
    -- I could imagine this being slow, but so far it seems ok
    facetsQ := for x in annotatedFaces Q list if x#0 =!= dim Q - 1 then continue else set x#2;
    for t in T do (
        if any(facetsQ, f -> isSubset(t, f)) then continue else return false;
        );
    true
    )
isTriangulationOfPolytope(ReflexivePolytope, Triangulation) := (Q, T) -> isTriangulationOfPolytope(Q, max T)

computeTriangulations = Q -> (
    -- sets "ptriangulations", "vtriangulations" (these are the ones that are *not* in ptriangulations.
    -- also sets "allPtriangulations", "allVtriangulations"
    -- maybe: FRSTs, FRVTs
    --  where FRSTs are the fine regular star triangulations, that is, the triangulations of the
    --   point set given by the rays of Q.
    --  FRVTs: fine regular vector triangulations, that are not in the above list.
    -- If we have computed all triangulations, we also set allFRVTs to true.  If we have computed some,
    --  but are not sure if we have them all,we set it to false.  Not set at all: means we have not
    --  computed triangulations yet.
    if Q.cache#?"allFRVTs" and Q.cache#"allFRVTs" then return;
    -- these are currently all FINE, REGULAR triangulations.  TODO question: should we include non-fine
    -- triangulations too?
    -- TODO: maybe have a max number we stash?
    vtris := findAllSimplicialFans(transpose matrix rays Q);
    H := partition(t -> isTriangulationOfPolytope(Q,t), vtris, {true, false});
    Q.cache#"vtriangulations" = sort for t in H#false list t; -- there had better be at least one here (why!?)
    Q.cache#"ptriangulations" = sort for t in H#true list t; -- there had better be at least one here.
    Q.cache#"allFRVTs" = true; -- currently, we have not coded the partial computation of these fans..
    )

findAllFRSTs ReflexivePolytope := List => Q -> (
    if not Q.cache#?"allFRVTs" then computeTriangulations Q;
    Q.cache#"ptriangulations"
    )

findAllFRVTs = method()
findAllFRVTs ReflexivePolytope := List => Q -> (
    if not Q.cache#?"allFRVTs" then computeTriangulations Q;
    Q.cache#"vtriangulations"
    )

findOneFRST = method()
findOneFRST ReflexivePolytope := Q -> (
    -- TODO: this is NOT correct if any rays interior to facets 
    last regularStarTriangulation polytope Q -- the 'last' removes the list of vertices
    )

restrictTriangulation = method()
restrictTriangulation(ZZ, ReflexivePolytope, List) := List => (d, Q, tri) -> (
    -- given Q, we use its annotated faces and the given triangulation, to
    -- write down the triangulation of the dim d-faces of the
    -- corresponding reflexive polytope in the N lattice side.
    -- (as a list of lists of d+1 elements)
    F := annotatedFaces Q;
    dfaces := for x in F list if x#0 =!= d then continue else x#2;
    sort unique flatten for td in dfaces list (
        a := set td; -- these are the indices we want.
        atri := sort unique for t in tri list (
            b := sort toList(a * set t);
            if #b == d+1 then b else continue
            );
        atri
        )
    )

-- This one is messed up:
-- findOneFRVT = method() -- this might return a FRST?
-- findOneFRVT ReflexivePolytope := List => Q -> (
--     -- returns one FRVT (or maybe FRST).
--     A := transpose matrix rays Q;
--     topcomRegularFineTriangulation(A, Homogenize => false)
--     )


-- TODO to get this running:
--  1. calabiYau needs to accept a ReflexivePolytope.
--  2. restrictTriangulation needs to work on (Q, tri).
--  3. rename `restrictTriangulation CalabiYauInToric`.   This gives more data, useful for intersection rings...
--  4. CalabiYauInToric needs: better way to get toricMoriConeCap
-- This function creates all the Batyrev CY's, not those coming from FRVTs
-- this function needs:
--   calabiYau(ReflexivePolytope, ...)
--   normalizeByAutomorphisms: doesn't use CYPolytope, so is fine as is.
--   restrictTriangulation
partitionFRSTsByDFaceEquivalence = method(Options => {Automorphisms => true})
partitionFRSTsByDFaceEquivalence(ZZ, ReflexivePolytope) := HashTable => opts -> (d, Q) -> (
    Ts := findAllFRSTs Q;
    gPerms := if opts.Automorphisms then 
                 automorphismsAsPermutations Q
              else 
                 {splice{0..#rays Q - 1}}; -- only the identity permutation
    f := if d < dim Q then 
             (tri -> normalizeByAutomorphisms(gPerms, restrictTriangulation(d, Q, tri)))
         else 
             (tri -> normalizeByAutomorphisms(gPerms, tri));
    partition(f, Ts)
    )

-- TODO: working on this.  Use partitionFRSTsByDFaceEquivalence above to help here.
findAllCYs ReflexivePolytope := List => opts -> Q -> (
    RZ := if opts#Ring === null then (
        a := getSymbol "a";
        h11 := hh^(1,1) Q;
        ZZ[a_1 .. a_h11]
        )
    else (
        opts#Ring
        );
    H := partitionFRSTsByDFaceEquivalence(dim Q - 2, Q);
    count := 0;
    Xs := for k in sort keys H list (
        -- TODO: this is perhaps a good place to compute toricMoriConeCap...
        X := calabiYau(Q, H#k#0, Ring => RZ, ID => count); -- take the first one
        count = count+1;
        X);
    Xs
    )

-- REMOVE
-- findAllFRVTs ReflexivePolytope := List => Q -> (
--   cQ := chirotope(transpose matrix rays Q, Homogenize => false);
--   "foo-topcomfoo.in" << toString cQ << endl << close;
--   ctris := for L in lines get ("!chiro2finetriangs <foo-topcomfoo.in") list (
--       matches := regex("\\{\\{[0-9,\\{\\}]*", L);
--       if matches === null or #matches != 1 then error "my logic is wrong";
--       value substring(matches_0, L)
--       );
--   ctris)

-- REMOVE
-- findAllFRVTs CYPolytope := List => Q -> (
--     alltris := topcomAllTriangulations(transpose matrix rays Q, ConnectedToRegular => false, Fine => true, Homogenize => false, RegularOnly => false);
--     -- now we need to choose the ones that are regular...
--     -- this seems to be a good way to do it.
--     select(alltris, t -> isProjective normalToricVariety(rays Q, t))
--     )

-*
-- FIX
generateTriangulations Triangulation := opts -> T -> (
    allT := new MutableHashTable;
    allT#T = true;
    TODO := {T};
    while #TODO > 0 and #(keys allT) < opts.Limit do (
        nextTRI := TODO#0;
        TODO = drop(TODO,1);
        --flips := affineCircuits nextTRI;
        --flips := select(affineCircuits nextTRI, z -> #z#0 > 1 and #z#1 > 1);
        flips := select(affineCircuits nextTRI, z -> #z#0 > 1 or #z#1 > 1);
        fliptris := for f in flips list bistellarFlip(nextTRI, f);
        newT := select(fliptris, x -> x =!= null);
        for T in newT do (
            if not allT#?T then (
                --<< "new triangulation: " << T << endl;
                if not opts.RegularOnly or isRegularTriangulation T then (
                    allT#T = true;
                    TODO = append(TODO, T);
                    );
                ));
        if debugLevel > 0 then 
            << "todo = " << #TODO << " and #triang = " << #(keys allT) << endl;
        );
    keys allT
    )

generateTriangulations(Matrix, List) := opts -> (Amat, tri) -> (
    allT := new MutableHashTable;
    allT#tri = true;
    TODO := {tri};
    while #TODO > 0 and #(keys allT) < opts.Limit do (
        nextTRI := TODO#0;
        TODO = drop(TODO,1);
        --flips := affineCircuits nextTRI;
        --flips := select(affineCircuits nextTRI, z -> #z#0 > 1 and #z#1 > 1);
        flips := select(affineCircuits(Amat, tri), z -> #z#0 > 1 and #z#1 > 1);
        fliptris := for f in flips list bistellarFlip(nextTRI, f);
        newT := select(fliptris, x -> x =!= null);
        for T in newT do (
            if not allT#?T then (
                --<< "new triangulation: " << T << endl;
                if not opts.RegularOnly or isRegularTriangulation(Amat, T, Homogenize => false) then (
                    allT#T = true;
                    TODO = append(TODO, T);
                    );
                ));
        if debugLevel > 0 then 
            << "todo = " << #TODO << " and #triang = " << #(keys allT) << endl;
        );
    keys allT
    )
*-

computeBasics ReflexivePolytope := Q -> (
    computeRaysAndDimensions Q;
    annotatedFaces Q; -- compute annotated faces
    computeGLSM Q;
    computeH11H21 Q;
    computeAutomorphisms Q;
    computeTriangulations Q;
    )

TEST ///
-*
  restart
  needsPackage "StringTorics"
*-
///

end--

restart
needsPackage "StringTorics"
ReflexivePolytope
topes = kreuzerSkarke 3;
A = matrix topes_40
P1 = convexHull A
P2 = polar P1
elapsedTime Q = reflexivePolytope(P2, ID => 40); -- .26 seconds. now .03 sec (as it doesn't get faces or lattice points...)
elapsedTime netList annotatedFaces Q
vertices Q
transpose matrix vertices Q
polytope(Q, "M")
polar Q
latticePoints Q
dump Q
reflexivePolytope oo
peek oo.cache
peek Q.cache
reflexivePolytope topes_55
