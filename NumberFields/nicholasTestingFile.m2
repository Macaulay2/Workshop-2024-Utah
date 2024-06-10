loadPackage ("NumberFields", Reload=>true)

R1 = QQ[x,z]/(x^2+1,z^2-2)
use R1
f1 = x
minimalPolynomial f1
NF1 = numberField R1
ring NF1
SE1 = simpleExt NF1
ring SE1

-- Running through simpleExt with the R1 example
K1 = ring NF1
describe K1
D1 = degree K1
d1 = 0
-- First pass
r1 = random(1, K1) -- Get a random homogeneous RingElement from K1 of degree 1
R1 = QQ[xx]
phi1 = map(K1, R1, {r1})
isPrime (kernel phi1)
I1 = kernel phi1 * sub(( 1/(((coefficients (first entries gens kernel phi1)_0)_1)_0)_0), R1)
simpleExt1 = numberField(R1/I1)
describe ring simpleExt1
d1 = degree simpleExt1 -- If d1==degree K1, we've found our simple extension!