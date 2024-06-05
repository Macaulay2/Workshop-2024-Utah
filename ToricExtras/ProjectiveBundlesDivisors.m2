needsPackage "NormalToricVarieties"

--this method takes as input a list of toric divisors on a normal toric variety 
--and outputs the projectivization of the associated vector bundle as a normal toric variety

projectivizationOfBundle = method()

--check that underlying varieties of given divisors are same?

projectivizationOfBundle List := (listOfDivisors) -> (

--get input variety
inputVariety := variety(listOfDivisors#0);

-- make list of coefficients of divisors 
listOfDivisorCoefficients := apply(listOfDivisors, i -> entries(i));

--get rank of vector bundle to projectivize minus 1, call this r
extendingDimension := #listOfDivisors - 1;

--get cones and rays of input variety
inputConeList := max(inputVariety);
inputRayList := rays(inputVariety);

-- get number of rays in the fan of the input variety
numberOfInputRays := #inputRayList;

-- get dimension of N_R
inputFanDimension := #(inputRayList#0);

-- padd with zeros to include the rays of the original fan in N_R x R^r
extendedInputRayList := {};

for oldRay in inputRayList do(

    extendedInputRayList = append(extendedInputRayList, (oldRay | toList (extendingDimension : 0)));

);


--build coefficients on standard basis vectors which will be added to input rays

rayLiftingCoefficients := {};

for i from 1 to extendingDimension do(

rayLiftingCoefficients = append(rayLiftingCoefficients , listOfDivisorCoefficients#i - listOfDivisorCoefficients#0);

);

-- make an identity matrix padded by 0's on top

extendedStandardBasisMatrix := map(ZZ^inputFanDimension, ZZ^extendingDimension, 0) || id_(ZZ^extendingDimension);

-- add lifted original rays to new ray list
newRayList := {};

for j from 0 to (#extendedInputRayList - 1) do(

    oldRay := extendedInputRayList#j;

    for i from 0 to (extendingDimension - 1) do(
        oldRay = oldRay + entries((rayLiftingCoefficients#i)#j*extendedStandardBasisMatrix_i);

    );
    
    newRayList = append(newRayList , oldRay);
    
    
    );

-- make e_0 = -e_1 -e_2 - ... - e_r and add to new ray list

eZero := entries((-extendedStandardBasisMatrix)*vector(toList (extendingDimension : 1)));

-- add e_0 to new ray list

newRayList = append(newRayList , eZero);

--add padded standard basis vectors to new ray list

for i from 0 to (extendingDimension - 1) do(
        newRayList = append(newRayList , entries(extendedStandardBasisMatrix_i));

    );

 


-- make a list which indexes the standard basis vectors at the end of the new ray list

shiftedEis := toList(0..extendingDimension) + toList ((extendingDimension + 1) : numberOfInputRays);


-- construct cone list of projectivization

newConeList := {};

for i from 0 to (#inputConeList - 1) do(

    for j from 0 to extendingDimension do(

        newConeList = append(newConeList , inputConeList#i | drop(shiftedEis, {j , j}) );

    );




);

-- make the projectivization!

normalToricVariety(newRayList , newConeList)
);

end ------------------------------

restart
needsPackage "NormalToricVarieties"
load("ProjectiveBundlesDivisors.m2")

rayListPP1 = {{1} , {-1}}
coneListPP1 = {{1}, {0}}
PP1 = normalToricVariety (rayListPP1, coneListPP1)
D0PP1 = toricDivisor ( { 0 , 0}, PP1)
D1PP1 = toricDivisor ( {0 , 7} , PP1)
H7 = projectivizationOfBundle({D0PP1, D1PP1})
rays H7
max H7

rayListPP2 = {{1 , 0}, {0 , 1}, {-1, -1}}
coneListPP2 = {{0, 1}, {1, 2}, {2 , 0}}
PP2 = normalToricVariety (rayListPP2, coneListPP2)
D0PP2 = toricDivisor ( { 9, 3 , 2}, PP2)
D1PP2 = toricDivisor ( {1 , 4, 7} , PP2)
Y = projectivizationOfBundle({D0PP2, D1PP2})
rays Y
max Y


dim X
rays X
X1 = toricProjectiveSpace 1
rays X1 == rays X
max X1
rayListH = { {1 , 0} , { 0 , 1} , {-1, 3}, {0 , -1}}
coneListH = { { 0 , 1} , { 1, 2} , { 2 , 3} , {3 , 0} }
Y = normalToricVariety (rayListH , coneListH)
Y1 = hirzebruchSurface 3
max Y1


variety(D1)

divisorsToVariety()

-----
PP2RayList = {
    {1,0}, {0,1}, {-1,-1}
}
PP2ConeList = {
    {0,1},
    {1,2},
    {2,0}
}
PP2 = normalToricVariety(PP2RayList, PP2ConeList)
D1 = toricDivisor ( {0 , 7, 9} , PP2)
D0 = toricDivisor ( { 8 , 0, 1}, PP2)

listOfDivisors = { D0, D1 }
inputVariety = variety listOfDivisors#0
r = #listOfDivisors - 1
oldDim = #entries(listOfDivisors#0)
-- build e_i
-- get original rays and then add zeros to the end
origRays = rays inputVariety
-- whichever side we add 0s to for the original rays, we need
-- to add to the complementary side for the e_is
origRaysInBigSpace = rays | { }