-- In this file, we have code to interface with CYTools as well as the CYTools databases for
-- e.g. CY3's with h11=4.

CYFromCYToolsDB = new Type of HashTable
cyFromCYToolsDB = method()
cyFromCYToolsDB(String, Sequence, Ring) := (DB, lab, RZ) -> (
    -- DIR should be a directory where to find the files containing computed info.
    result := new CYToolsCY3 from {
        "Directory" => DB,
        "Label" => lab,
        cache => new CacheTable
        };
    result.cache.PicardRing = RZ;
    result
    )

getCYFilePrefix = method()
getCYFilePrefix CYToolsCY3 := X -> (
     (h21,polynum,cynum) := label X;
     nm := X#"Directory" | "h21_"|h21|"/poly_"|polynum|"/CY_"|cynum
    )

getPolytopeFilePrefix = method()
getPolytopeFilePrefix CYToolsCY3 := X -> (
     (h21,polynum,cynum) := label X;
     nm := X#"Directory" | "h21_"|h21|"/poly_"|polynum
    )

findH21s = DIR -> (
    names := select(readDirectory DIR, s -> match("h21_", s));
    sort for n in names list (m := regex("h21_([0-9]+)", n); value substring(m_1, n))
    )
findPolys = (DIR, h12) -> (
    DIRNAME := DIR | "h21_" | h12 | "/";
    names := select(readDirectory DIRNAME, s -> match("poly_", s));
    sort for n in names list (m := regex("poly_([0-9]+)", n); value substring(m_1, n))
    )
findCYs = (DIR, h12, polyid) -> (
    DIRNAME := DIR | "h21_" | h12 | "/poly_" | polyid | "/" ;
    --print DIRNAME;
    --print readDirectory DIRNAME;
    names := select(readDirectory DIRNAME, s -> match("CY_", s));
    sort for n in names list (m := regex("CY_([0-9]+)", n); value substring(m_1, n))
    )

findAllCYToolsCY3s = method()
findAllCYToolsCY3s(String, Ring) := (DB, RZ) -> (
    -- RZ must be a ring over ZZ with #vars being h11 of the CY.
    -- all of these are expected to have the same h11.
    hashTable flatten flatten for h21 in findH21s DB list 
        for pol in findPolys(DB, h21) list
            for cy in findCYs(DB, h21, pol) list 
                (h21, pol, cy) => cyFromCYToolsDB(DB, (h21, pol, cy), RZ)
    )

label CYToolsCY3 := X -> X#"Label"

picardRing CYToolsCY3 := X -> X.cache.PicardRing

hh(Sequence, CYToolsCY3) := (pq, X) -> (
    lab := label X;
    if pq === (1,1) then 
        numgens picardRing X
    else if pq === (1,2) or pq === (2,1) then
        lab#0
    else
        error "hh only implemented for hh^(1,1) and hh^(1,2), hh^(2,1)"
    )

c2Form CYToolsCY3 := X -> (
    filename := (getCYFilePrefix X) | "/c2.dat";
    if not fileExists filename then error "c2.dat doesn't exist for this CY";
    coeffs := value get filename; -- format is one line, with e.g. [28, 12, 10, -6]
    RZ := picardRing X;
    sum for i from 0 to numgens RZ - 1 list coeffs#i * RZ_i
    )    


cubicForm CYToolsCY3 := X -> (
    filename := (getCYFilePrefix X) | "/IntersectionVariety.dat";
    if not fileExists filename then error "IntersectionVariety.dat doesn't exist for this CY";
    use picardRing X;
    value get filename -- result is a polynomial
    )    

isFavorable CYToolsCY3 := X -> (
    filename := (getPolytopeFilePrefix X) | "/favorable.dat";
    if not fileExists filename then error "favorable.dat doesn't exist for this CY";
    contents := get filename;
    if match("True", contents) then true
    else if match("False", contents) then false
    else error "cannot parse data in favorable.dat"
    )

-- TODO: should we keep this?  Problem: some of the database files have it, some don't
-- isToric = method()
-- isToric CYToolsCY3 := X -> (
--     filename := (getCYFilePrefix X) | "/toric.dat";
--     if not fileExists filename then error "toric.dat doesn't exist for this CY";
--     contents := get filename;
--     if match("True", contents) then true
--     else if match("False", contents) then false
--     else error "cannot parse data in toric.dat"
--     )    

invariantsAll CYToolsCY3 := X -> invariantsAll(c2Form X, cubicForm X, hh^(1,1) X, hh^(1,2) X)
invariantsAllX CYToolsCY3 := X -> invariantsAll(c2Form X, cubicForm X, hh^(1,1) X, hh^(1,2) X)

hasGVInvariants = method()
hasGVInvariants CYToolsCY3 := X -> (
    filename := (getCYFilePrefix X) | "/GVs.dat";
    fileExists filename
    )

gvInvariants CYToolsCY3 := optsUnused -> X -> (
    -- this just grabs the ones that are in the data base.
    filename := (getCYFilePrefix X) | "/GVs.dat";
    if not fileExists filename then error "GVs.dat doesn't exist for this CY";
    Ls := lines get filename;
    Ls/value/(x -> toList drop(x, -1) => x#-1)//hashTable
    )

partitionGVConeByGV CYToolsCY3 := HashTable => opts -> X -> (
    -- return null if we cannot computr GV invariants (i.e. if non-favorable).
    if not hasGVInvariants X then return null;
    gv := gvInvariants(X, opts); -- TODO: stash this?
    C := posHull transpose matrix ((keys gv)/toList);
    gvX := entries transpose rays C;
    partition(f -> if gv#?(toSequence f) then gv#(toSequence f) else 0, gvX)
    )

heft CYToolsCY3 := List => X -> (
    filename := (getCYFilePrefix X) | "/gradingVec.dat";
    if not fileExists filename then error "gradingVec.dat doesn't exist for this CY";
    contents :=  get filename;
    contents = replace("\\[ +", "[", contents);
    toList value replace(" +", ",",contents)
    )

gvInvariantsAndCone(CYToolsCY3, ZZ) := Sequence => opts -> (X, D) -> (
    -- D is the degree bound to start with.  We could start with 5, or DegreeLimit/2 or DegreeLimit/4, or ...
    -- returns a hash table of computed GV invariants, and the cone generated by all curves with nonzero GV invariant.
    if not hasGVInvariants X then return null;
    degvec := heft X;
    gv := gvInvariants(X, opts);
    keysgv := keys gv;
    H := hashTable for k in keysgv list k => dotProduct(k, degvec);
    firstSet := select(keys H, k -> H#k <= D);
    if debugLevel > 0 then << "The number of curves in the first set: " << #firstSet << endl;
    C := posHull transpose matrix (firstSet);
    Cdual := dualCone C;
    HC := transpose rays Cdual;
    curves := for k in keys H list if H#k > D and H#k <= opts.DegreeLimit then transpose matrix {k} else continue;
    set2 := select(curves, c  -> any(flatten entries (HC * c), a -> a < 0));
    if debugLevel > 0 then << "The number of curves not in the first cone: " << #set2 << endl;
    C2 := if #set2 == 0 then C else posHull (rays C | matrix{set2});
    if debugLevel > 0 and #set2 == 0 then (
        << "CY " << label X << " C = " << rays C  << endl
        )
    else
        << "*differs* CY " << label X << " C1 = " << rays C << " and C2 = " << rays C2 << endl;
    (gv, C2)
    )

partitionGVConeByGV(CYToolsCY3, ZZ) := HashTable => opts -> (X, D) -> (
    -- return null if we cannot compute GV invariants (i.e. if non-favorable).
    if not hasGVInvariants X then return null;
    (gv, C) := gvInvariantsAndCone(X, D, opts);
    gvX := entries transpose rays C;
    partition(f -> if gv#?f then gv#f else 0, gvX)
    )

partitionGVConeByGV(CYToolsCY3) := HashTable => opts -> (X) -> (
    -- return null if we cannot computr GV invariants (i.e. if non-favorable).
    if not hasGVInvariants X then return null;
    deg := opts.DegreeLimit;
    D := if deg <= 10 then deg else D = 10;
    (gv, C) := gvInvariantsAndCone(X, D, opts);
    gvX := entries transpose rays C;
    partition(f -> if gv#?f then gv#f else 0, gvX)
    )

moriConeCap = method()
moriConeCap CYToolsCY3 := Cone => X -> (
    -- the columns are the extremal rays of the cone
    filename := (getCYFilePrefix X) | "/mori_gens.dat";
    if not fileExists filename then return "mori_gens.dat doesn't exist for this CY";
    contents :=  get filename;
    contents = replace("\\[ +", "[", contents);
    contents = replace(" *\\] *", "]", contents);
    contents = replace(" +", ",",contents);
    mori := (lines contents)/value/toList;
    posHull transpose matrix mori
    )

rays CYToolsCY3 :=  List => {} >> o -> X -> (
    filename := (getPolytopeFilePrefix X) | "/points.dat";
    if not fileExists filename then return "points.dat doesn't exist for this CY";
    rys := (lines get filename)/value/toList; -- format: one ray per line (with comma between values), first one is origin
    if not all(rys#0, a -> a == 0) then error "my logic is wrong: the first element should be the origin";
    drop(rys, 1) --remove the origin.
    )

max CYToolsCY3 := List => X -> (
    filename := (getCYFilePrefix X) | "/simplices.dat";
    if not fileExists filename then return "simplices.dat doesn't exist for this CY";
    rys := (lines get filename)/value/toList; -- format: one ray per line (with comma between values), first one is origin
    rys
    )

basisIndices CYToolsCY3 := List => X -> (
    filename := (getPolytopeFilePrefix X) | "/basis.dat";
    if not fileExists filename then return "basis.dat doesn't exist for this CY";
    basIndices := toList value get filename; -- format: one ray per line (with comma between values), first one is origin
    basIndices
    )
    
cyPolytope CYToolsCY3 := opts -> X -> (
    pts := rays X;
    P2 := convexHull transpose matrix pts;
    Q := cyPolytope(P2, ID => (label X)_1);
    X.cache.CYPolytope = Q;
    --assert(hh^(1,1) X == hh^(1,1) Q);
    --assert(hh^(1,2) X == hh^(1,2) Q);
    Q
    )

-- Database exploration and GV computation scratch code moved to ScratchCYTools.m2

