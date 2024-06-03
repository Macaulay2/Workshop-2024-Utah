needsPackage "NormalToricVarieties"

divisorsToVariety := method();

--check that underlying varieties of given divisors are same

divisorsToVariety List := (listOfDivisors) -> (

inputVariety := variety(listOfDivisors#0);


-- listOfDivisors = { {0 , 7, 9}, {8 , 0 , 1} }
-- extract coefficients
--  get r = numberOfDivisors
-- build e_i
-- set e_0 = -e_1 - ... - e_r
-- Construct Fi's: Fi = cone(e_0, ... , exclude e_i , ... , e_r )
-- 

listOfDivisorCoefficients := apply(listOfDivisors, i -> entries(i));

numberOfDivisors := #listOfDivisors - 1;

inputConeList := max(inputVariety);
inputRayList := rays(inputVariety);
inputFanDimension := #(inputRayList#0);

extendedInputRayList := {};

for oldRay in inputRayList do(

extendedInputRayList = append(extendedInputRayList, oldRay | {numberOfDivisors : 0});

)





-- for each cone sigma in inputConeList, for each ray rho in sigma
-- construct cone(inputRayList#rho +  )
-- 

-- extendedStandardBasisMatrix := map(ZZ^inputFanDimension, ZZ^numberOfDivisors, 0) || id_(ZZ^numberOfDivisors);

rayLiftingCoefficients := {};
-- build rayliftingCoeffs
for i from 1 to numberOfDivisors do(

rayLiftingCoefficients = append(rayLiftingCoefficients , listOfDivisorCoefficients#i - listOfDivisorCoefficients#0);

);

extendedRayLiftingMatrix = map(ZZ^inputFanDimension, ZZ^numberOfDivisors, 0) || matrix(rayLiftingCoefficients)

newRayList := {};

for oldRay in extendedInputRayList do(

    newRayList = append(newRayList , 
    
    
    );
    



    );
);

 

newConeList := {};
newVariety := normalToricVariety(newRayList , newConeList);

);

end ------------------------------

restart
needsPackage "NormalToricVarieties"
rayList = {{-1} , {1}}
coneList = {{1}, {0}}
X = normalToricVariety (rayList, coneList)
dim X
rays X
X1 = toricProjectiveSpace 1
X == X1
rays X1 == rays X
max X1
rayListH = { {1 , 0} , { 0 , 1} , {-1, 3}, {0 , -1}}
coneListH = { { 0 , 1} , { 1, 2} , { 2 , 3} , {3 , 0} }
Y = normalToricVariety (rayListH , coneListH)
Y1 = hirzebruchSurface 3
max Y1

D1 = toricDivisor ( {0 , 7} , X)
D0 = toricDivisor ( { 0 , 0}, X)
variety(D1)

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