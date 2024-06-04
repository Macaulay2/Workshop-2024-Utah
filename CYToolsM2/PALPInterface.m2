newPackage(
    "PALPInterface",
    Version => "0.1",
    Date => "",
    Headline => "An interface to parts of the PALP polyhedra software",
    Authors => {{ Name => "", Email => "", HomePage => ""}},
    AuxiliaryFiles => false,
    DebuggingMode => true,
    PackageImports => {"ReflexivePolytopesDB"}
    )

export {
    "palpVertices","getVerticesFromWS", "Drange"
    }

executableDir = "/Users/frenly/Desktop/palp-2.21/";

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

getVerticesFromWS = method()

getVerticesFromWS String := Matrix => w -> (
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
    matrix M)

getWSFromDim = method(Options => {Drange => null})

getWSFromDim ZZ := List => opts -> d -> (
    drange := opts.Drange;
    str0 := toString(d);
    str1 := toString(drange_0);
    str2 := toString(drange_1);
    str3 := if opts.Drange === null then "!cws.x -w"|str0
    else "!cws.x -w"|str0|" "|str1|" "|str2;

    PALPOutput := get str3;
    L := lines PALPOutput;
    M := for ell in L list(
	L0 := separate(" +",ell);
	L0Mod := take(L0, {0,d+1});
	for x in L0Mod list value x
	)
    )
///
getWSFromDim(3,Drange => {5,8})
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
testwm8x5dw = stringToList("8 5")


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

testwm8x5 = linesToMatrix(
"
-1 0 0 0 0
0 -1 0 0 0
0 0 -1 0 0
0 0 0 -1 0
0 0 0 0 1
0 0 0 0 -1
0 0 0 1 -1
1 1 1 0 -3
"
)



--getNOfPartitions = method()

--getNOfPartitions := w -> (
    
--)

--getNOfPartitions(2, testwm8x5dw, testwm8x5)










--programPaths#"PALP" = executableDir;
PALP = findProgram("PALP", "poly.x -h")

--programPaths#"poly.x" = executableDir;
polyx = findProgram("poly.x", "poly.x -h")

-* Code section *-
palpMatrix = method()
palpMatrix Matrix := String => M -> (
    if ring M =!= ZZ then error "expected integer matrix";
    header := toString(numRows M | " " | numcols M);
    e := entries M;
    es := for e1 in e list for a in e1 list (toString a | " ");
    s := concatenate between("\n", es);
    header | "\n" | s | "\n"
    )

runPoly = method()
-- runPoly(Matrix, String) := (M, opts) -> (
--     "foo" << palpMatrix M << close;
--     cmd := "poly.x -" | opts | " foo";
--     runProgram(PALP, cmd)
--     )

runPoly(Matrix, String) := (M, opts) -> (
    "foo" << palpMatrix M << close;
    cmd := " -" | opts | " foo";
    result := runProgram(polyx, cmd);
    matrix KSEntry result#"output"
    )


runPoly(Matrix, String) := String => (M, opts) -> (
    "foo" << palpMatrix M << close;
    cmd := " -" | opts | " foo";
    result := runProgram(polyx, cmd);
    result#"output"
    )

normalForm = method()
normalForm Matrix := Matrix => (M) -> (
    "foo" << palpMatrix M << close;
    cmd := " -N foo";
    result := runProgram(polyx, cmd);
    matrix KSEntry result#"output"
    )

-* Documentation section *-
beginDocumentation()

end--


doc ///
Key
  PALPInterface
Headline
Description
  Text
  Tree
  Example
  CannedExample
Acknowledgement
Contributors
References
Caveat
SeeAlso
Subnodes
///

doc ///
Key
Headline
Usage
Inputs
Outputs
Consequences
  Item
Description
  Text
  Example
  CannedExample
  Code
  Pre
ExampleFiles
Contributors
References
Caveat
SeeAlso
///

-* Test section *-
TEST /// -* [insert short title for this test] *-
-- test code and assertions here
-- may have as many TEST sections as needed
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
palpMatrix M
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


