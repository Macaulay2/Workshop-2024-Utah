--------------------------------------------------------------
-- CalabiYauInToric (soon to change back to CalabiYauInToric? ----------
--------------------------------------------------------------

CalabiYauInToric.synonym = "Calabi-Yau hypersurface in a normal toric variety"
CalabiYauInToric.GlobalAssignHook = globalAssignFunction
CalabiYauInToric.GlobalReleaseHook = globalReleaseFunction
expression CalabiYauInToric := X -> if hasAttribute(X, ReverseDictionary) 
     then expression toString getAttribute(X, ReverseDictionary) else 
     (describe X)#0
net CalabiYauInToric := X -> net expression X     
describe CalabiYauInToric := X -> Describe (
    "A Calabi-Yau "|dim X|"-fold hypersurface with h11="|hh^(1,1) X|" and h21="|hh^(1,2) X |" in a "|(dim X + 1)|"-dimensional toric variety"
    )

CYPolytope.synonym = "Calabi-Yau reflexive polytope"
CYPolytope.GlobalAssignHook = globalAssignFunction
CYPolytope.GlobalReleaseHook = globalReleaseFunction
expression CYPolytope := X -> if hasAttribute (X, ReverseDictionary) 
    then expression getAttribute (X, ReverseDictionary) else 
    (describe X)#0
describe CYPolytope := X -> Describe (expression CYPolytope) (
    expression rays X, expression max X)


CYDataFields = {
    -- first entry: true means it must exist and be in the main hash table
    --   false: it might exist, and is in the cache table.
    "polytope data" => {value, Q -> toString Q.cache#"id", CYPolytope},
    "triangulation" => {value, toString, List}
    }

-- These are the cache fields that we write to a string via 'dump'
CYDataCache = {
    -- first entry: true means it must exist and be in the main hash table
    --   false: it might exist, and is in the cache table.
    "id" => {value, toString, ZZ},
    "c2" => {value, toString, List},
    "intersection numbers" => {value, toString, List},
    "toric intersection numbers" => {value, toString, List},
    "toric mori cone cap" => {value, toString, List}
    }

setCYIntersectionRing = (X, R) -> (
    -- X is a CalabiYauInToric
    -- R is a polynomial ring, or null (if not, an error is raised).
    n := hh^(1,1) X;
    if R =!= null then (
        if not instance(R, PolynomialRing) or numgens R != n then 
            error ("expected polynomial ring with "|n|" variables");
        X.cache.PicardRing = R;
        )
    else (
        a := getSymbol "a";
        X.cache.PicardRing = ZZ[a_1..a_n];
        );
    )

calabiYau = method(Options => {ID => null, Ring => null})
-- TODO, BUG!! The triang needs to indices in the Q rays.
calabiYau(CYPolytope, List) := CalabiYauInToric => opts -> (Q, triang) -> (
    X := new CalabiYauInToric from {
        symbol cache => new CacheTable,
        "polytope data" => Q,
        "triangulation" => triang
        };
    if opts.ID =!= null then X.cache#"id" = opts.ID;
    setCYIntersectionRing(X, opts#Ring);
    X
    )

cyData = method(Options => options calabiYau)
cyData(CYPolytope, List) := opts -> (Q, triang) -> calabiYau(Q, triang, opts)

picardRing = method()
picardRing CalabiYauInToric := X -> X.cache.PicardRing

cyData(String, Function) :=
calabiYau(String, Function) := CalabiYauInToric => opts -> (str, F) -> (
    -- F is a function which takes an id of a CYPolytope and returns the CYPolytope
    -- The string is the value taken from a CY database .
    L := lines str;
    if L#0 != "CYData" then error "string is not in proper format";
    fields := hashTable for i from 1 to #L-1 list getKeyPair L#i;
    -- First get the main elements (these are required!):
    polytopeid := value fields#"polytope data";
    required := for field in CYDataFields list (
        k := field#0;
        if k === "polytope data" then (
            "polytope data" => F polytopeid
            )
        else (
            readFcn := field#1#0;
            if fields#?k then k => readFcn fields#k else error("expected key "|k)
        ));
    X := new CalabiYauInToric from prepend(symbol cache => new CacheTable, required);
    -- now read in the cache values (including "id" value, if any)
    for field in CYDataCache do (
        k := field#0;
        readFcn := field#1#0;
        if fields#?k then X.cache#k = readFcn fields#k;
        );
    if opts.ID =!= null then X.cache#"id" = opts.ID; -- just for compatibility with other constructors...
    setCYIntersectionRing(X, opts#Ring);
    X
    )

dump CalabiYauInToric := String => {} >> opts -> X -> (
    s1 := "CYData\n";
    strs := for field in CYDataFields list (
        k := field#0;
        writerFunction := field#1#1;
        if not X#?k then error("expected key: "|k#0);
        "  " | k | ":" | writerFunction(X#k) | "\n"
        );
    strs2 := for field in CYDataCache list (
        k := field#0;
        writerFunction := field#1#1;
        if not X.cache#?k then continue;
        "  " | k | ":" | writerFunction(X.cache#k) | "\n"
        );
    strs = join({s1}, strs, strs2);
    concatenate strs
    )

makeCY = method(Options => {ID => null, Ring => null})
makeCY CYPolytope := CalabiYauInToric => opts -> Q -> (
    P2 := polytope Q;
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    if rays Q =!= LP then error "I have a lattice point mismatch";
    cyData(Q, tri, opts)
    )    

makeCY(List, List) := CalabiYauInToric =>  opts -> (pts, triangulation) -> (
    -- We keep the translation around?
    Q := cyPolytope pts;
    -- now we need the translation from old vertices to new.
    H := hashTable for i from 0 to #rays Q - 1 list (rays Q)#i => i;
    mapping := hashTable for i from 0 to #pts-1 list (
        p := pts#i;
        if all(p,a -> a == 0) then continue; -- leave out the origin
        if H#?p then i => H#p else 
            error("lattice point found which is likely interior to a facet: "|(toString p))
        );
    tri := sort for t in triangulation list (
        sort for t1 in drop(t,1) list mapping#t1
        );
    Q.cache#"vertex translation" = mapping;
    cyData(Q, tri, opts)
    )


normalToricVariety CalabiYauInToric := opts -> X -> (
    if not X.cache.?NormalToricVariety then X.cache.NormalToricVariety = (
        Q := X#"polytope data";
        T := X#"triangulation";
        GLSM := transpose matrix degrees Q;
        normalToricVariety(rays Q, T, opts, WeilToClass => matrix GLSM)
        );
    X.cache.NormalToricVariety
    -- TODO: this fails if the class group is torsion! (Fails: later it gives an inscrutable error...)
    )

rays CalabiYauInToric :=  List => {} >> o -> X -> rays cyPolytope X
max CalabiYauInToric := X -> X#"triangulation"

-- TODO: triangulation is used with 2 different pieces of data:
--  with, without cone point!  Change this to use only one point.
-- Also: there are 4 matrices one can imagine: A, A0 (A with origin), Ah, A0h...
-- We need to be consistent about these!
-- TODO: do we really need this?
triangulation CalabiYauInToric := Triangulation => opts -> X -> (
    if not opts.Homogenize then error "Homogenize flag is not used in this method";
    if not X.cache#?"triangulation" then (
        rys := X#"polytope data"#"rays";
        d := #rys#0;
        B := (transpose matrix rys) | matrix{d:{0}};
        X.cache#"triangulation" = triangulation(B, for t in X#"triangulation" list append(t, #rys)); -- TODO: BUG?? where is "triangulation" key? In cache??
        );
    X.cache#"triangulation"
    )

cyPolytope CalabiYauInToric := opts -> X -> X#"polytope data"
dim CalabiYauInToric := X -> dim ambient X - 1
polytope CalabiYauInToric := X -> polytope cyPolytope X
polytope(CalabiYauInToric, String) := (X, which) -> polytope(cyPolytope X, which)
basisIndices CalabiYauInToric := List => X -> basisIndices cyPolytope X
degrees CalabiYauInToric := List => X -> degrees cyPolytope X

ambient CalabiYauInToric := X -> normalToricVariety X

label = method()
label CYPolytope := Q -> if Q.cache#?"id" then Q.cache#"id" else ""
label CalabiYauInToric := X -> (label cyPolytope X, if X.cache#?"id" then X.cache#"id" else "")

hh(Sequence, CalabiYauInToric) := (pq, X) -> hh^pq cyPolytope X

isFavorable CalabiYauInToric := Boolean => X -> isFavorable cyPolytope X

abstractVariety CalabiYauInToric := opts -> X -> (
    -- Store this with X.
    V := ambient X;
    aX := completeIntersection(V, {-toricDivisor V});
    abstractVariety(aX, base())
    )
abstractVariety(CalabiYauInToric, AbstractVariety) := opts -> (X, pt) -> (
    -- Store this with X?
    -- Check: pt is of dimension zero?
    V := ambient X;
    aX := completeIntersection(V, {-toricDivisor V});
    abstractVariety(aX, pt)
    )

-- restrictTriangulation: returns a List of
--   {2-face indices, 
--    all indices of points in in this 2-face, 
--    the triangles in this 2-face, 
--    genus of this face}
-- TODO: need also a function which returns just: triangles, genus information.
restrictTriangulation = method()
restrictTriangulation CalabiYauInToric := List => (X) -> (
    -- given X, we use its annotated faces and its triangulation, to write down the triangulations of the 2-faces
    -- of the corresponding reflexive polytope in the N lattice side.
    Q := cyPolytope X;
    F := annotatedFaces Q;
    twofaces := for x in F list if x#0 =!= 2 then continue else {x#1, x#2, x#4};
    T := max X; -- triangulation
    for t2 in twofaces list (
        a := set t2#1; -- these are the indices we want.
        atri := sort unique for t in T list (
            b := sort toList(a * set t);
            if #b == 3 then b else continue
            );
        {t2#0, t2#1, atri, t2#2}
        )
    )

restrictTriangulation(ZZ, CalabiYauInToric) := List => (d, X) -> (
    -- given X, we use its annotated faces and its triangulation, to
    -- write down the triangulations of the dim d-faces of the
    -- corresponding reflexive polytope in the N lattice side.
    Q := cyPolytope X;
    F := annotatedFaces Q;
    dfaces := for x in F list if x#0 =!= d then continue else x#2;
    T := max X; -- triangulation
    sort unique flatten for td in dfaces list (
        a := set td; -- these are the indices we want.
        atri := sort unique for t in T list (
            b := sort toList(a * set t);
            if #b == d+1 then b else continue
            );
        atri
        )
    )

lineBundle(CalabiYauInToric, List) := (X, deg) -> (
    if not all(deg, x -> instance(x, ZZ)) or #deg =!= degreeLength ring ambient X
    then error("expected multidegree of length "|degreeLength ring ambient X);
    new LineBundle from {
        symbol cache => new CacheTable,
        symbol variety => X,
        symbol degree => deg
        }
    )

degree LineBundle := L -> L.degree
variety LineBundle := L -> L.variety

installMethod(symbol _, OO, CalabiYauInToric, LineBundle => 
     (OO,X) -> lineBundle(X, (degree 1_(ring ambient X)))
     )

LineBundle Sequence := (L, deg) -> (
    lineBundle(variety L, degree L + toList deg)
    )

equations CalabiYauInToric := List => X -> (
    if not X.cache.?Equations then X.cache.Equations = (
        V := ambient X;
        {random(degree(-toricDivisor V), ring V)} -- TODO: (1) allow tuned equations, (2) do the random call more efficiently.
        );
    X.cache.Equations
    )

toricMoriCone = method()
toricMoriConeCap = method()

setToricMoriConeCap = method()
setToricMoriConeCap(CalabiYauInToric, List) := List => (Y, Xs) -> (
    -- does nothing if Y is not favorable
    -- Xs are all of the CY3's equivalent to Y (including Y), but NO others.
    --myNTFE := restrictTriangulation Y;
    --myXs := select(Xs, X0 -> restrictTriangulation X0 === myNTFE);
    if not isFavorable Y then null
    else
        Y.cache#"toric mori cone cap" = sort entries transpose rays dualCone posHull matrix{for X in Xs list rays dualCone toricMoriCone X}
    )
setToricMoriConeCap CalabiYauInToric := List => Y -> (
    -- Xs are all of the CY3's equivalent to Y (including Y), possibly includes others too?
    if not isFavorable Y then return null;
    Q := cyPolytope Y;
    Xs := findAllCYs(Q, NTFE => false, Automorphisms => false);
    myNTFE := restrictTriangulation Y;
    myXs := select(Xs, X0 -> restrictTriangulation X0 === myNTFE);
    Y.cache#"toric mori cone cap" = sort entries transpose rays dualCone posHull matrix{for X in myXs list rays dualCone toricMoriCone X};
    )
    
toricMoriConeCap CalabiYauInToric := List => Y -> (
    if not isFavorable Y then return null;
    if not Y.cache#?"toric mori cone cap" then setToricMoriConeCap Y;
    Y.cache#"toric mori cone cap"
    )

///
  -- Let's test restrictTriangulation, automorphisms and getting only triangulations
  -- unique up to linear automorphism of the lattice polytope.
  topes = kreuzerSkarke(3, Limit => 20);
  topes_15
  tope = KSEntry "4 7  M:74 7 N:8 7 H:3,63 [-120] id:15
   1   0   0   0  -3  -3   3
   0   1   1   1   1  -2  -2
   0   0   3   0   3   0  -6
   0   0   0   3  -3  -6   6
  "
  Q = cyPolytope tope
  Xs = findAllCYs Q
  #Xs
  auts = for e in automorphisms Q list matrix e
  rayvecs = for v in rays Q list transpose matrix {v}
  rayHash = hashTable for i from 0 to #rayvecs - 1 list rayvecs#i => i
  g = auts#0
  perms = for g in auts list
    for v in rayvecs list rayHash#(g*v)
  sort unique flatten restrictTriangulation(3, Xs#0)
  sort for x in oo list sort perms_0_x
///
