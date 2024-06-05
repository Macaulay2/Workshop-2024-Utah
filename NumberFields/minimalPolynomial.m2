loadPackage ("NumberFields", Reload=>true)

-- Nicholas' minimal polynomial file.
-- Method to compute minimal polynomial via pushFwd.
-- Eventually extend to also compute minimal polynomials for a list of elements.
minimalPolynomialViaPushFwd = method(Options => {})
minimalPolynomialViaPushFwd(RingElement) := opts -> f1 -> (
    R1 := ring f1;
    S1 := (coefficientRing(R1))[y];
    P1 := pushFwd(map(R1, QQ));
    A1 := (P1#2)(1_R1);
    pow1 := 1;
    while not gens(kernel(A1))==0 do (
        A1 |= (P1#2)(f1^pow1);
        pow1 += 1;
    );
    M := matrix({{1_S1}});
    for i1 from 1 to pow1 do (
        M |= y^i1;
    );
    (entries (M*K))#0#0
)

loadPackage ("NumberFields", Reload=>true)
R1 = QQ[x]/(x^3-2)
S1 = (coefficientRing(R1))[y]
elt1 = x
--phi1 = map(R1, QQ)
P1 = pushFwd(map(R1, QQ))
g1 = (P1#2)(1_R1)
g2 = (P1#2)(elt1)
g3 = (P1#2)(elt1^2)
g4 = (P1#2)(elt1^3)
A = g1|g2|g3|g4
K = gens(kernel(A))
M = matrix({{1, y, y^2, y^3}})
(entries (M*K))#0#0

--****************
--Tests
--****************

loadPackage ("NumberFields", Reload=>true)
-- This one should return y^3 - 2
R1 = QQ[x]/(x^3-2)
f = x
minimalPolynomial f

-- This one should return y^3 - 4
f = x^2
minimalPolynomial f

-- This one should return y^3 - 6*y - 6
f = x^2+x
minimalPolynomial f

minimalPolynomial {x, x^2, x^2+x}

-- This one should return y^4 - 2*y^2 + 9
R1 = QQ[x,y]/(x^2+1,y^2-2)
f = x+y
minimalPolynomial f

simpleExt numberField R1