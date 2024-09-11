newPackage(
    "PALPInterface",
    Version => "0.1",
    Date => "",
    Headline => "An interface to parts of the PALP polyhedra software",
    Authors => {
        {Name => "Mike Stillman", 
            Email => "mike@math.cornell.edu", 
            HomePage => "http://pi.math.cornell.edu/~mike"}
        },
    AuxiliaryFiles => false,
    DebuggingMode => true,
    PackageImports => {"ReflexivePolytopesDB", "Polyhedra", "NormalToricVarieties"}
    )

export {
    "palpVertices","getVerticesFromWS", "getWSFromDim",
    "stringToList", "linesToMatrix", "matrixToLines", 
    "getPartitionInfo", "readPartitions", "formDivisor",
    "getPolyhedralInfo",
    "palpMatrix", "fromPalpMatrix", -- perhaps do not export these?
    "normalForm",
    "runPoly", -- TODO: change name, or remove this function, once we understand output of PALP better.
    "runNEF",
    "Hodge",
    "Normal"
    }

--programPaths#"PALP" = executableDir;
PALP = findProgram("PALP", "poly.x -h")

--programPaths#"poly.x" = executableDir;
POLYX = findProgram("poly.x", "poly.x -h")
NEFX = findProgram("nef.x", "nef.x -h")

-- translate Matrix into a String suitable for input to PALP
-- warning: often, the convex hull of the columns must be reflexive.
-- we need to document when this is necessary
palpMatrix = method()
palpMatrix Matrix := String => M -> (
    if ring M =!= ZZ then error "expected integer matrix";
    header := toString(numRows M | " " | numcols M);
    e := entries M;
    es := for e1 in e list for a in e1 list (toString a | " ");
    s := concatenate between("\n", es);
    header | "\n" | s | "\n"
    )

-- From output of palp functions, we might need to grab only some of the lines
-- to form a matrix.
fromPalpMatrix = method()
fromPalpMatrix String := Matrix => str -> (
    matrix KSEntry str -- this will give a somewhat inscrutable error if the format is not correct
    )

TEST ///
-*
  restart
  needsPackage "PALPInterface"
*-
  M = matrix{{1,1,1,1},{0,1,2,3}}
  str1 = palpMatrix M
  assert(M == fromPalpMatrix str1)

  str2 = palpMatrix transpose M
  assert(transpose M == fromPalpMatrix str2)

  -- from smoothFanoToricVariety(3, 5)
  M = transpose matrix {{1, 0, 0}, {-1, 0, 1}, {0, 1, 0}, {0, -1, 1}, {0, 0, 1}, {0, 0, -1}}
  assert(M == fromPalpMatrix palpMatrix M)
///

runPoly = method()
runPoly(Matrix, String) := String => (M, opts) -> (
    "foo" << palpMatrix M << close;
    cmd := " -" | opts | " foo";
    result := runProgram(POLYX, cmd);
    result#"output"
    )

runNEF = method()
runNEF(Matrix, String) := String => (M, opts) -> (
    "foo" << palpMatrix M << close;
    cmd := " -" | opts | " foo";
    result := runProgram(NEFX, cmd);
    result#"output"
    )

normalForm = method()
normalForm Matrix := Matrix => (M) -> (
    "foo" << palpMatrix M << close;
    cmd := " -N foo";
    result := runProgram(POLYX, cmd);
    print result#"output";
    matrix KSEntry result#"output"
    )

TEST ///
-*
  restart
  needsPackage "PALPInterface"
*-
  -- 3d reflexive example
  -- from smoothFanoToricVariety(3, 5)
  M = transpose matrix {{1, 0, 0}, {-1, 0, 1}, {0, 1, 0}, {0, -1, 1}, {0, 0, 1}, {0, 0, -1}}
  assert(M == fromPalpMatrix palpMatrix M)
  normalForm M

  needsPackage "Polyhedra"
  P = convexHull M
  P' = polar P
  M1 = lift(vertices P', ZZ)

  needsPackage "IntegerEquivalences"
  A = extendToMatrix{3,4,5}
  A = matrix {{1, -2, 1}, {2, 1, -2}, {-1, 1, 0}}
  assert(det A == 1)
  assert(normalForm (A * M1) == normalForm M1)

  runPoly(M1, "p")
  runPoly(M, "p")

  runPoly(M1, "v")
  runPoly(M, "v")

  runPoly(M, "e")
  runPoly(M1, "e")  

  runPoly(M, "m")
  runPoly(M1, "m")  

  runPoly(M, "g") -- what does output mean? -- M:7 6 N:29 8 Pic:17 Cor:0
  runPoly(M1, "g") -- M:29 8 N:7 6 Pic:3 Cor:0

  runPoly(M, "a")
  runPoly(M1, "a")

  runPoly(M, "i")
  runPoly(M1, "i")

  runPoly(M, "I")
  runPoly(M1, "I")
  runPoly(matrix{{1,1,1,1},{0,1,2,3}}, "I") == "No IP\n"
  runPoly(transpose matrix {{1, 0, 0}, {-1, 0, 1}, {0, 1, 0}, {0, -1, 1}, {0, 0, 1}, {0, 0, -1}}, "I")

  runPoly(M, "S")
  runPoly(M1, "S")

  runPoly(M, "vT")
  runPoly(M, "pT")
  runPoly(M1, "vT")
  runPoly(M1, "pT")

  runPoly(M, "V")
  runPoly(M1, "V")

  runPoly(M, "B")
  runPoly(M1, "B")

  runPoly(M, "F")
  runPoly(M1, "F")

  runPoly(M, "A")
  runPoly(M1, "A")

  runPoly(M, "G")
  runPoly(2*M1, "G")
  ///

nefPartitions = method()
nefPartitions Matrix := M -> (
    -- M should be the matrix whose columns are vertices of the reflexive polytope in N lattice.
    
    )
getPolyhedralInfo = method(Options => {Hodge => true, Normal => false})
getPolyhedralInfo Matrix := opts -> M -> (
    -- idea: create the matrix of M.
    -- create the call to PALP
    -- call PALP (poly.x here)
    -- parse the output and return the answer as a hash table of desired info.
    w := palpMatrix M;
    -- this needs to be changed depending on desired info
    str1 := "!poly.x -v  << FOO\n";
    --str1 := " -v << FOO\n";
    str2 := "\nFOO\n";
    str3 := str1|w|str2;
    PALPOutput := get str3;
    -- the following obtains the matrix remaining after removing the given lines.
    keeplines := select(lines PALPOutput, s -> (
            not match("^Degrees", s) and not match("^Type", s) and not match("or", s)));
    str := concatenate between("\n", keeplines);
    matrix KSEntry str
    )
getPolyhedralInfo Matrix := opts -> M -> (
    -- idea: create the matrix of M.
    -- create the call to PALP
    -- call PALP (poly.x here)
    -- parse the output and return the answer as a hash table of desired info.
    w := palpMatrix M;
    -- this needs to be changed depending on desired info
    str1 := "!poly.x -g  << FOO\n";
    --str1 := " -v << FOO\n";
    str2 := "\nFOO\n";
    str3 := str1|w|str2;
    PALPOutput := get str3;
    return PALPOutput;
    -- the following obtains the matrix remaining after removing the given lines.
    keeplines := select(lines PALPOutput, s -> (
            not match("^Degrees", s) and not match("^Type", s) and not match("or", s)));
    str := concatenate between("\n", keeplines);
    matrix KSEntry str
    )

-- str = get "!poly.x -v -r << FOO
-- 10 1 2 3 4
-- FOO
-- "

-- L = lines str
-- L = drop(drop(L, 3), -2)
-- netList L
-- L0 = separate(" +", L_0)
-- for x in drop(L0, 1) list value x
-- M = matrix for ell in L list (
--     L0 := separate(" +", ell);
--     for x in drop(L0, 1) list value x
--     )
-- needsPackage "Polyhedra"
-- needsPackage "StringTorics"
-- P = convexHull M
-- vertices P
-- isReflexive P

-- from Matrix to palp form of a matrix.
-- from palp form of a matrix to Matrix.
-- from ws or cws to matrix whose columns are the vertices of the corresponding reflexive polytope
-- from reflexive polytope vertices to ws or cws... (?)
-- optional: HodgeNumbers: h11, h12, ..., h1(dim X - 1).
-- nef partitions
-- 


getVerticesFromWS = method()

getVerticesFromWS List := Matrix => wlist ->(
    w1 := for i in wlist list(toString(i)|" ");
    w := concatenate(drop(w1,-1),toString(wlist_(-1)));
    
    str1 := "!poly.x -v -r << FOO\n";
    str2 := "\nFOO\n";
    str3 := str1|w|str2;
    PALPOutput := get str3;
    L := lines PALPOutput;
    L = drop(drop(L, 3), -2);
    M := for ell in L list (
        L0 := separate(" +", ell);
        for x in L0 list if x!="" then value x else continue
        );
    if #M == 0 then return null;
    matrix M)

getVerticesFromWS({10,1,2,3,4})

getWSFromDim = method(Options => {Degrees => null})

getWSFromDim ZZ := List => opts -> d -> (
    drange := opts.Degrees;
    str0 := toString(d);
    str1 := if opts.Degrees === null then "!cws.x -w"|str0
    else "!cws.x -w"|str0|" "|toString(drange_0)|" "|toString(drange_1);

    PALPOutput := get str1;
    L := lines PALPOutput;
    L1 := if opts.Degrees === null then L else drop(L,-1);
    M := for ell in L1 list(
	L0 := separate(" +",ell);
	L0Mod := take(L0, {0,d+1});
	for x in L0Mod list value x
	)
    )

///
-*
  restart
  needsPackage "PALPInterface"
*-
  wss = getWSFromDim 3
  assert(#wss == 95)
  ws = wss_90
  assert(ws == {44, 4, 5, 13, 22})  

  M = getVerticesFromWS ws
  wt = transpose matrix{drop(ws, 1)}
  M * wt
  M2 = transpose LLL syz transpose wt  

  normalForm M2

  
  getWSFromDim(4)

wss = getWSFromDim(5, Degrees => (20,20))
assert(#wss == 61)
V = getVerticesFromWS wss_58
needsPackage "Polyhedra"
P = convexHull V
isReflexive P
latticePoints P
isSimplicial polar P
///

--matrix oo
--getVerticesFromWS "10 1 2 3 4"
--oo_0    
--oo/class


stringToList = method()
stringToList String := List => stemp -> (
    L0 := separate(" +", stemp);
    (for x in L0 list if x!="" then value x else continue)
)
-- testwm8x5dw = stringToList("8 5")


linesToMatrix = method()
linesToMatrix String := Matrix => stemp -> (
--    << stemp;
    L := lines stemp;
    M := for ell in L list (
--	<< "line " << ell << "\n";
	if ell == "" then continue;
        stringToList(ell)
        );
    matrix M
)




matrixToLines = method()
matrixToLines List := String => ll -> (
-- matrixToLines := ll -> (
    res := "";
    for ell in ll do (
        respart := "\n";
        for x in ell do respart = respart|toString(x)|" ";
        -- << "eachline: " << respart << " \n";
        res = res|respart;
    );
    res
);



getPartitionInfo = method()
getPartitionInfo (ZZ, ZZ, ZZ) := String => (dimen, indexing, cod) -> (
    V := smoothFanoToricVariety(dimen, indexing);
    -- << rays V;
    -- << # (rays V);
    -- << dim V;
    nrowcol := toString(# (rays V))|" "|toString(dim V);
    -- << cod << " \n";
    entrylines := matrixToLines(rays V);


    -- cod, nrowcol, entrylines
--    commp1 := "!nef.x -N -c"|cod|" -p << FOO\n";
    commp1 := "!nef.x -N -c"|cod|" << FOO\n";
    commp3 := "\nFOO\n";
    -- << commp1|nrowcol|entrylines|commp3;
    test := get (commp1|nrowcol|entrylines|commp3);
    -- << "result: ";
    -- << test;

    test
)

readPartitions = (output1) -> (
    reg1 := "P:[0-9 ]+V:([0-9 ]+ {1})";
    for thisline in (lines output1) list(
        -- reg1 = "M:(.*)N(.*)codim(.*)part*";
        
        -- (for x in L0 list if x!="" then value x else continue)
        regres := regex(reg1, thisline);
        -- << "start: " << thisline;
        -- << "reg: ";
        -- << thisline;
        -- << "\n";
        -- -- if not (regres === null) then << regres;
        -- << regex(reg1, thisline);
        -- << "\n";
        -- << "\n";
        -- << output1_(219, 4);
        if not (regres === null) then (
            << thisline << " " << thisline_(regres_1) <<"\n";
            -- for eachparti in thisline_(regres_1) list (regex("^[0-9]$", eachparti))
            stringToList(thisline_(regres_1))
        )
        else continue
    )
)

-- TODO: What is this function doing?
formDivisor = (V, partlist) -> (
    raylist := rays V;
    -- << # raylist;
    -- -- << raylist_0;
    -- << raylist;
    -- << "\n";
    -- << partlist;
    -- << "\n";
    fulllist := toList (0..(# raylist - 1));
    partlist = reverse partlist;
    complist := fulllist;

    for thispart in partlist do complist = drop(complist, {thispart, thispart});
    -- << complist;
    D1list := for thispart in partlist list (V_thispart);
    -- << "D1list: " << for thispart in partlist list ("V_"|toString(thispart));
    D1 := sum(D1list);
    -- << D1;
    
    D2list := for thispart in complist list (V_thispart);
    -- << "D2list: " << for thispart in complist list ("V_"|toString(thispart));
    D2 := sum(D2list);
    -- << D2;

    -- << isNef(D1);
    -- << isNef(D2);
    if not isNef(D1) then (<< "D1 false");
    if not isNef(D2) then (<< "D2 false");
    (D1, D2)
    -- completeIntersection(V, {D1, D2})

    -- test5 = (for thispart in partlist list drop(fulllist, {thispart, thispart}));
    -- test5
    -- complist = 
    -- for thispart in partlist do << thispart;


)










-- runPoly(Matrix, String) := (M, opts) -> (
--     "foo" << palpMatrix M << close;
--     cmd := "poly.x -" | opts | " foo";
--     runProgram(PALP, cmd)
--     )

-- runPoly(Matrix, String) := (M, opts) -> (
--     "foo" << palpMatrix M << close;
--     cmd := " -" | opts | " foo";
--     result := runProgram(POLYX, cmd);
--     matrix KSEntry result#"output"
--     )



-* Documentation section *-
beginDocumentation()

doc ///
Key
  PALPInterface
Headline
  interface to the polyhedral program PALP
Description
  Text
    The program PALP was designed by ... to help compute all reflexive polytopes in dimensions 3 and 4.
    However, it has some generally useful functionality beyond that.
///

TEST ///
-*
  restart
  needsPackage "PALPInterface"
*-
  assert true
///

-- template for doc nodes for methods/functions
///
Key
Headline
Usage
Inputs
Outputs
Description
  Text
  Example
SeeAlso
///

TEST ///
-*
  restart
  needsPackage "PALPInterface"
*-
  M = matrix({{ -1,  -1,  -1,  -1,   1},
       { -1,  -1,   0,   3,  -1},
       { -1,  -1,   4,   0,  -1},
       { -1,   2,  -1,   2,  -1},
       { -1,   6,   0,  -1,  -1},
       {  1,   2,  -1,   2,  -1},
       {  1,   6,   0,  -1,  -1},
       { -1,  -1,  -1,  -1,  -1},
       { -1,  -1,  -1,   3,  -1},
       { -1,  -1,   5,  -1,  -1},
       { -1,   0,  -1,   3,  -1},
       { -1,   2,   3,  -1,  -1},
       { -1,   3,   1,   0,  -1},
       { -1,   7,  -1,  -1,  -1},
       {  1,   2,   3,  -1,  -1},
       {  5,   0,  -1,   3,  -1},
       {  7,  -1,   5,  -1,  -1},
       {  7,   7,  -1,  -1,  -1},
       { 23,  -1,  -1,   3,  -1},
       {151,  -1,  -1,  -1,  -1}})
  needsPackage "Polyhedra"
  P = convexHull transpose M
  isReflexive P
  runPoly(M, "gD")

  runNEF(M, "c2")

  wss = getWSFromDim(5, Degrees => (30,30))
  #wss == 354
  ws = wss#10
  for ws in wss list (
    M = getVerticesFromWS(ws);
    if M === null then continue;
    if not isReflexive convexHull M then << "note: " << ws << " does not give reflexive" << endl;
    nefs = runNEF(M, "c2");
    print nefs;
    ws => nefs)

  wss = getWSFromDim(5, Degrees => (10,10))
  #wss === 5
  for ws in wss list (
    M = getVerticesFromWS(ws);
    if M === null then continue;
    if not isReflexive convexHull M then << "note: " << ws << " does not give reflexive" << endl;
    nefs = runNEF(M, "c2");
    print nefs;
    ws => nefs)

///

TEST ///
  restart
  needsPackage "PALPInterface"


  ws = {10, 1, 1, 1, 1, 3, 3}  
  needsPackage "StringTorics"
  M = getVerticesFromWS ws
  P = convexHull M
  dim P
  V = reflexiveToSimplicialToricVariety P
  isSimplicial V
  isSmooth V
  transpose matrix rays V
  max V
  assert isWellDefined V
  normalForm transpose matrix rays V

  P2 = polar P
  vertices P2
  normalForm lift(vertices P2, ZZ)
  (LP, tri) = regularStarTriangulation(3, P2)
  V = normalToricVariety(LP, tri)
  isSimplicial V
  isSmooth V -- yes!

  m = lift(matrix vertices P2, ZZ)

  vertsP2 = entries transpose vertices P2
  raysV = rays V
  
  runNEF(m, "N -c2")
  runNEF(matrix rays V, "N -c2")
  -- this gives us the following nef partition of V
  D1 = V_2 + V_3
  D2 = -toricDivisor V - D1
  isNef D1
  isNef D2

  - degree toricDivisor V
  degree D1
  degree D2
  transpose matrix LP

  SR = dual monomialIdeal V
  S1 = ZZ/32003[(gens ring V)]
  minimalBetti sub(SR, S1)
  picardGroup V
  classGroup V

  X = completeIntersection(V, {D1, D2})
  hodgeDiamond X -- h11=3, although induced is only ZZ^2.
  
  -- if  ws = {30, 1, 2, 4, 4, 7, 12}
  -- these don't work so well...
  --elapsedTime HH^1(V, OO_V(0, 0, 1, 0, 1, 0, 1, 0, 0, 0))
  --elapsedTime cohomCalg(V, V_2 + V_4 + V_8)

  elapsedTime cohomCalg(V, V_2 + V_4 + V_5)
///

end--

-* Development section *-
restart
debug needsPackage "PALPInterface"
check "PALPInterface"

uninstallPackage "PALPInterface"
restart
installPackage "PALPInterface"
viewHelp "PALPInterface"

restart
debug needsPackage "PALPInterface"
M = matrix{{1,1,1,1},{0,1,2,3}}
palpMatrix transpose M
getPolyhedralInfo M
getPolyhedralInfo transpose M



"foo" << palpMatrix M << close
runPoly(M, "v")
normalForm M
viewHelp runProgram

needsPackage "StringTorics"
kss = kreuzerSkarke(7, Limit => 10)

L = vertices polar convexHull matrix kss_7
L = sub(L, ZZ)
L1 = normalForm L

first lines runPoly(L1, "g") -- informational line
runPoly(L1, "p")
runPoly(L1, "v")
runPoly(L1, "e")
runPoly(L1, "m")
runPoly(L1, "d")
runPoly(L1, "a")
runPoly(L1, "D")
runPoly(L1, "i")
runPoly(L1, "I")
runPoly(L1, "S")
runPoly(L1, "Tv")
runPoly(L1, "Tp")
runPoly(L1, "N")
runPoly(L1, "t")

-- Using nef.x
3 1 1 1 0 0 0 0 0  2 0 0 0 1 1 0 0 0  3 0 0 0 0 0 1 1 1
nef.x -h
nef.x -p -c2 -V << FOO
3 1 1 1 0 0 0 0 0  2 0 0 0 1 1 0 0 0  3 0 0 0 0 0 1 1 1
FOO

-- Using cws.x
cws.x -h

10 1 1 1 1 1 5 M:1128 6 N:8 6 H:1,0,976 [5910]

# points are the columns
poly.x -v << FOO
5 1 1 1 1 1
FOO

# points are the rows
poly.x -e << FOO
5 1 1 1 1 1
FOO

# points are the rows
poly.x -v << FOO
45 5 6 7 8 9 10
FOO

run "poly.x -e << FOO
3 1 1 1 0 0 0 3 0 0 0 1 1 1
FOO
"

run "poly.x -e << FOO
36 1 4 4 6 9 12
FOO
"

str = get "!poly.x -v -r << FOO
36 1 4 4 6 9 12
FOO
"

str = get "!poly.x -v -r << FOO
10 1 2 3 4
FOO
"

L = lines str
L = drop(drop(L, 3), -2)
netList L
L0 = separate(" +", L_0)
for x in drop(L0, 1) list value x
M = matrix for ell in L list (
    L0 := separate(" +", ell);
    for x in drop(L0, 1) list value x
    )
needsPackage "Polyhedra"
needsPackage "StringTorics"
P = convexHull M
vertices P
isReflexive P

get "!poly.x -e << FOO
42 2 3 5 5 6 21
FOO
"

get "!poly.x -e << FOO
143233 43 1651 3328 20226 47194 70791
FOO
"

get "!poly.x -v << FOO
143233 43 1651 3328 20226 47194 70791
FOO
"

lines get///!curl "http://rgc.itp.tuwien.ac.at/fourfolds/db/5d_reflexive,h11=100.txt"///;
#oo == 204046
o31_100000




restart
needsPackage"PALPInterface"
M = getVerticesFromWS("10 1 2 3 4")

M2 = getVerticesFromWS("3 1 1 1 0 0 0  3 0 0 0 1 1 1")

needsPackage "StringTorics"
P = convexHull(M2)
isReflexive P
#latticePoints P


---------- nef-partitions ---------
restart
needsPackage "PALPInterface"
needsPackage "StringTorics"
needsPackage "NormalToricVarieties"
viewHelp NormalToricVarieties

V = smoothFanoToricVariety(3, 4)
rays V
max V
V_0
V_1
V_4
toricDivisor V
S = ring V
describe S
2*V_3 + 10*V_2
degree V_3
degree V_1
degree(V_1 + V_3)
degree (-toricDivisor V)
basis({3,2}, S)
F = random({3,2}, S)
size F
isNef V_0
isNef V_1
isNef V_4

V = smoothFanoToricVariety(5, 6)
rays V
transpose matrix degrees ring V
dim V

nef.x -N -c2 -p << FOO
8 5
-1 0 0 0 0
0 -1 0 0 0
0 0 -1 0 0
0 0 0 -1 0
0 0 0 0 1
0 0 0 0 -1
0 0 0 1 -1
1 1 1 0 -3
FOO

D1 = V_4 + V_6 + V_7
D2 = V_0 + V_1 + V_2 + V_3 + V_5
isNef(V_4 + V_6 + V_7)
degree(V_4 + V_6 + V_7)
isNef(V_0 + V_1 + V_2 + V_3 + V_5)
degree(V_0 + V_1 + V_2 + V_3 + V_5)
for i from 0 to 7 list isNef V_i
for i from 0 to 7 list if i == 5 then continue else isNef(V_i + V_3 +  V_5)

X = completeIntersection(V, {D1, D2})
hodgeDiamond X
pt = base(a,b)
Xa = abstractVariety(X, pt)
IX = intersectionRing Xa
describe IX
basis(1, IX)
h = a * t_6 + b * t_7
integral(h^3)
integral((chern_2 tangentBundle Xa) * h)





---------- nef-partitions ---------
restart
needsPackage "PALPInterface"
needsPackage "StringTorics"
needsPackage "NormalToricVarieties"

testinfo = getPartitionInfo(5, 6, 2)
Y = smoothFanoToricVariety(5, 6)
testPartitionList = readPartitions(testinfo)
-- (testD1, testD2) = formDivisor(smoothFanoToricVariety(5, 6), {4, 6, 7})
divisorsList = for thispartcomb in testPartitionList list (
    (testD1, testD2) = formDivisor(Y, thispartcomb)
)
oo/(x -> (x/isNef))
-- X = completeIntersection(testV, {D1, D2});
X = completeIntersection(Y, toList divisorsList_0);
hodgeDiamond X
