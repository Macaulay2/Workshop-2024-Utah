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