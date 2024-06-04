loadPackage ("NumberFields", Reload=>true)

-- First example: compute splitting field of x^4+5*x^2+6
-- Use spaces instead of tabs
R = QQ[x]
f = x^4-5*x^2+6
L = decompose ideal f
currentEntry = entries gens L#0
S = R/(L#0)[y]

phi1 = map(S, R, {y})
f = phi1(f)
L = decompose ideal f
R = S
S = R/(L#2)[z]

phi1 = map(S, R, {z})
f = phi1(f)
L = decompose ideal f
length (entries gens L#0)#0
max degree (entries gens L#0)#0#0

T

f = x_1^4-5*x_1^2+6;

-- splittingField method
splittingField = method(Options => {})
splittingField(RingElement) := opts -> f1 -> (
    R1 = QQ[x_1];
    variableIndex = 2;
    finished = false;
    while not finished do (
        L = decompose ideal f1;
        finished = true;
        executeForLoop = true;
        for i from 0 to length(L)-1 do (
            if executeForLoop then (
                currentEntry = (entries gens L#i)#0;
                if not (length(currentEntry)==1 and max(degree(currentEntry#0))==1) then (
                    finished = false;
                    S1 = R1/(L#i)[x_variableIndex];
                    phi1 = map(S1, R1, {x_variableIndex});
                    f1 = phi1(f1);
                    R1 = S1;
                    variableIndex += 1;
                    executeForLoop = false;
                )
            )
        )
    )
    new NumberField from toField (R1)
)

decompose ideal(f)

polynomialSplits = f1 -> (
    L1 = decompose ideal f1;
    splits = true;
    for i from 0 to #L1-1 do (
        currentEntry := (entries gens L1#i)#0;
        if not (length(currentEntry)==1 and max(degree(currentEntry#0))==1) then (splits = false; break;)
    );
    splits
)

--*******************
-- To confirm that a polynomial splits in its splitting field
loadPackage ("NumberFields", Reload=>true)
R = QQ[x]
f1 = x^3+1
L1 = decompose ideal f1
polynomialSplits f1

S = splittingField f1
T = (ring target S)[y]
--T = (flattenRing ring S)#0[y]
f11 = (map(T, R, {y})) f1
L2 = decompose ideal f11
polynomialSplits f11

--TODO for Nicholas: make splittingField return a NumberField and then make it work when the parameter polynomial comes from 
--something like QQ[x]/(x^2-2). Actually, make splittingField return a NumberFieldExtension.
--TODO: make a test to assert that the degrees of some field 
--x^5+x+1 over QQ[x]
--TODO: instead of calling decompose ideal f every time, instead use synthetic division to divide f by the roots we already have.

loadPackage ("NumberFields", Reload=>true)
R = QQ[w]/(w^2+1)[x] -- R = QQ[i][x]
f1 = x^2-2
L1 = decompose ideal f1
polynomialSplits f1

S = splittingField f1
T = (ring target S)[y]
f11 = (map(T, R, {y})) f1
L2 = decompose ideal f11
polynomialSplits f11

--x^5-4*x+2 over QQ[x]
loadPackage ("NumberFields", Reload=>true)
R = QQ[x]
f1 = x^5-4*x+2
L1 = decompose ideal f1
polynomialSplits f1

S = splittingField f1
T = (ring target S)[y]
f11 = (map(T, R, {y})) f1
L2 = decompose ideal f11
polynomialSplits f11