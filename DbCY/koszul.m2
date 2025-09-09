restart
needs "./GlobalExt.m2"
nefGenerators ProjectiveVariety := X -> matrix{{1}}

R = ZZ/32003[x,y,z]/(x^2,x*y)
C = koszulComplex matrix{{x,y}}
apply(3, i -> prune sheaf HH_i C)
globalExt(C, C)
globalExt(0, C, C)
ExtTable(Proj R, {C})

---

S = ZZ/32003[x,y,z]
C = koszulComplex matrix{{x^2,x*y}}
minimize part_0 Hom(freeResolution truncate_2 S^1, C)

----

S = ZZ/32003[x,y,z,w]
C = naiveTruncation(koszulComplex vars S, (1,2))
prune HH C
globalExt(1, C, C)
ExtTable(Proj S, {C})


----

restart
needsPackage "Complexes"
load "GlobalExt.m2"
X = toricProjectiveSpace 3
S = ring X
I = monomialCurveIdeal(S, {1,3,4})
D = freeResolution I
C = complex(S^1);
globalExt(C, D)
globalExt(C, D[1])
globalExt(C[1], D)

for j in {-2,-1,0,1,2,3} do (
    D0 = D[-j];
    print globalExt(-j, C, D0) -- Ext^(-j)(C, D0)
    );
for j in {-2,-1,0,1,2,3} do (
    D0 = D[-j];
    print globalExt(1+j, C, D0)
    );
