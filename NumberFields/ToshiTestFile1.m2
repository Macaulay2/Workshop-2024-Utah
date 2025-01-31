----
isFieldAutomorphism(NF, F)
matrix(basis NF1)
print(flatten  entries (basis NF1))

------
loadPackage ("NumberFields", Reload=>true)
loadPackage("InvariantRing", Reload=>true)


R = QQ[w,v]/ ideal(w^3-2,v^2+v+1)
R = QQ[x]/ideal(x^2-2)
NF = numberField(R)
getGaloisGroup NF
--

getGaloisGroup= method(Options =>{});
getGaloisGroup(NumberField) := (nF) -> (
    R1 := nF[u];
    --We get all the roots and store them in rootList
    rootList := {};
    for i from 0 to (length gens coefficientRing R1)-1 do(
        r0 := (gens coefficientRing R1)_i;
        minPol := minimalPolynomial(r0);
        M := map(R1,ring minPol,{(gens R1)_0});
        rootList := append(rootList, getRoots(M(minPol)));
    );
    --We then generate the necessary permutations.
    perms := toList (0..(length rootList_0 - 1));
    for i from 1 to length rootList - 1 do (
        perms = perms ** toList (0..(length rootList_i - 1));
    );
    if length rootList == 1 then (
        perms = {};
        for i from 0 to length rootList_0 -1 do (
            perms = append(perms, {i})
        );
    );
    perms = perms / splice;
    --We then loop through, checking if they are field automorphisms, then adding them to all maps if they are.
    allMaps := {};
    for i from 0 to length perms -1 do (
        rootsImg := {}; 
        for j from 0 to length perms_i -1 do (
            rootsImg = append(rootsImg, substitute(rootList_j_(perms_i_j), nF));
        );
        M = map(nF, nF, rootsImg );
        F := matrixFromNumberFieldMap(M);
        if isFieldAutomorphism(nF, F) then (
            allMaps = append(allMaps, F);
        );
    );
    --We then create the group
    numVars := length flatten entries basis  nF;
    return group finiteAction(allMaps, QQ [x_1..x_numVars]);
)
--
R = QQ[w,v]/ ideal(w^3-2,v^2+v+1)

NF = numberField(R)
R1 = (NF)[u]
getRoots(u^3-2)
(basis NF)  

-------
R = QQ[w,v]/ ideal(w^3-2,v^2+v+1)
R = QQ[x]/ideal(x^2-2)
getGaloisGroup NF
NF = numberField(R)
R1 = (NF)[u]
r0 = (gens coefficientRing R1)_0
r1 = (gens coefficientRing R1)_1
------
 
---
length rootList_0 
numVars = length flatten entries basis  NF
length flatten entries basis  NF
R = QQ [x_0..x_numVars]

minPol = minimalPolynomial(r0)
getRoots(minPol)
rootList = {}
for i from 0 to (length gens coefficientRing R1)-1 do(
    r0 = (gens coefficientRing R1)_i;
    minPol = minimalPolynomial(r0);
    M = map(R1,ring minPol,{(gens R1)_0});
    toCoeffRing = map(coefficientRing R1, R1, {});
    rootList = append(rootList, getRoots(M(minPol)));
    --We now get the roots
)
print(rootList)
ring rootList_0_0
a1 = (gens NF)_0
a2 = (gens NF)_1
substitute(a1, NF)
rootList_0_0
allMaps = {}
for i from 0 to length rootList_0 -1 do (
    for j from 0 to length rootList_1 -1 do (
        M = map(NF, NF, {substitute(rootList_0_i, NF),substitute(rootList_1_j, NF)} );
        F = matrixFromNumberFieldMap(M);
        if isFieldAutomorphism(NF, F) then (
            allMaps = append(allMaps, F);
        );
    );
)
allMaps 
M = map(NF, NF, {substitute(rootList_0_1, NF),a2} )
F = matrixFromNumberFieldMap(M)
isFieldAutomorphism(NF, F)
NF 
NF 
G = finiteAction({allMaps_0,  allMaps_2} , QQ [a0,a1,a2,a3,a4,a5])
group G 
isAbelian G 
----------------------
F_2(a1)





F
for i from 0 to length rootList_0 do (
    for j from 0 to length rootList_1 do (
        M = pushFwd(map(NF, NF,  ))
    );
) 

for i from 0 to 0 do(
    r0 = (gens coefficientRing R1)_i;
    minPol = minimalPolynomial(r0);
    M = map(R1,ring minPol,{(gens R1)_0});
    toCoeffRing = map(coefficientRing R1, R1, {});
    rootList = getRoots(M(minPol));
    --We now get the roots 

    for j from 0 to ((length rootList)-1) do(
        print(j);
        print(rootList_j);
        a = toCoeffRing(rootList_j);
    
        vList = append(vList, NF#pushFwd#2(rootList_j));
    );
    --We now find our n linear independent vectors
    --We start by trying all the roots
    indepVctrs = {NF#pushFwd#2(1_NF)};
    for j from 0 to ((length rootList)-1) do (
        trialBasis = append(indepVctrs, NF#pushFwd#2(toCoeffRing(rootList_j)));
        print(trialBasis);
        if (checkLinIndep(trialBasis)) then(
            indepVctrs = trialBasis;
        );
    );
);

a = toCoeffRing(a)
t = terms toCoeffRing(a)
liftable (t_0 / (basis ring NF)_5_0, ring NF)
lift(t_0 / (basis ring NF)_5_0, ring NF)
(basis ring NF)_0_0
t_0
t_0
--We get roots in each variable
--We check each permutation of roots for variable to see if it's a field automorphism


--We expect a ring 

S = (QQ[x]/(x^2-2))[y]
f = y^2-2
getRoots(f)
R = (flattenRing S)#0
M = (flattenRing S)#1
Minv (flattenRing (S, Result=>3))#2
F = M(f) 
I = ideal F
primeFactors = decompose I
primeFactors_0_0 
(degree (primeFactors_0_0))#0 == 1
try lift (F / ((primeFactors_0_0)^2), R) then print("Hi") else print("hello")
m(F) 
decompose ideal (y^2-2)
z = (y^2-2) / (y-x)^2 


R = numberField(QQ[x]/(x^2-2))
R0 = ring R
a1 = (gens R0)#0  
R1 = R0[z]
f = z^2 - 2
decompose f 

factor f 
ring  f 
coefficientRing ring f 



gens R1 
()^3
RA = ring R 
x = (gens RA)_0
matrixFromRingEl(R, 1_RA)
f0 = pushFwd(map(RA^1, RA^1, matrix{{x+3*x^2-3}}))
ringElFromMatrix(R, f0)

 

--Matrix to field element
matrixFromRingEl = method();
matrixFromRingEl(NumberField, RingElement) := (nF, rEl) -> (
    R := ring nF;
    return pushFwd(map(R^1, R^1, matrix{{rEl}}));
)


ringElFromMatrix = method();
ringElFromMatrix(NumberField, Matrix) := (nF, mat) -> (
    --We basically turn the natural linear algebra basis of our number field into a matrix, then row reduce it to turn mat into an element in our number field.
    R0 := ring nF;
    R1 := coefficientRing R0;
    M0 := (pushFwd(map(R0,R1)))_1;
    vList := {};
    for i from 0 to ((numgens source M0)-1) do(
        M_i := matrixFromRingEl(R, (M0_i)_0);
        v_i := vector reshape(R1^((numgens target M_i)*(numgens source M_i)), R1^1, M_i);
        vList = append(vList, v_i);
    );
    vList = append(vList, vector reshape(R1^((numgens target mat)*(numgens source mat)), R1^1, mat));
    M := matrix(vList);
    RRM := reducedRowEchelonForm M;
    lastCol := RRM_{numColumns RRM-1};
    el := 0_R0;
    for i from 0 to ((numgens source M0)-1) do(
        el = el + (lastCol_0)_i * (M0_i)_0;
    );
    return el;
)














R = numberField(QQ[x]/(x^3-2))
R1 = ring R
R2 = coefficientRing R1
R3 = (pushFwd(map(R1,R2)))_1
P1 = (pushFwd(map(R1, R2)))#1
R = numberField(QQ[x]/(x^3-2))
x = (gens R1)_0
f0 = pushFwd(map(R1^1, R1^1, matrix{{x}}))
f1 = inverse f0
P1_0
(P1_1)_0
vList = {}
for i from 0 to ((numgens source P1)-1) do(
    print((P1_i)_0);
    M_i = matrixFromRingEl(R, (P1_i)_0);
    v_i = vector reshape(R2^((numgens target M_i)*(numgens source M_i)), R2^1, M_i);
    vList = append(vList, v_i);
)
vList = append(vList, vector reshape(R2^((numgens target f1)*(numgens source f1)), R2^1, f1))
M = matrix(vList)
reducedRowEchelonForm M
Matrix to field element





disc = method();
disc(NumberField) := (numField) ->(
    S = simpleExt(numField);
    R = ring S;
    print((presentation R));
);
R = numberField(QQ[x]/(x^3-2))
S = ring R
x = (gens S)_0
matrixFromRingEl(R, x)
-- describe S
gens S 
degree R 
f0 = pushFwd(map(S^1, S^1, matrix{{x}}))
f1 = pushFwd(map(S^1, S^1, matrix{{x^2}}))
f2 = pushFwd(map(S^1, S^1, matrix{{1_S}}))
f3 = pushFwd(map(S^1, S^1, matrix{{x^2+2*x}}))

inverse(f0)
inverse f1
inverse f2
inverse f3