restart
needsPackage "OrlovFunctors"

sheaf Complex := identity

-- Example 2.7
R = QQ[x_0,x_1,x_2] / sum(3, i -> x_i^3)
X = Proj R
a = x_0 + x_1
b = x_0^2-x_0*x_1+x_1^2
A = matrix {{x_2, a}, {b, -x_2^2}}
B = matrix {{x_2^2, a}, {b, -x_2}}
assert(A * B == 0 and B * A == 0)
A = map(R^{-1,0}, , A)
M = coker A
F = res(M, LengthLimit => 3)
F1 = orlovTruncateGeq(F, 1)
F1' = dual F1
G0 = canonicalTruncation(F1', -2, 0)
G' = dual freeResolution G0
G1' = orlovTruncateLess(G', 1)
G = sheaf G1'


-- Functorial example
N = x_0 * M
f = inducedMap(N, M, x_0 * id_M)
F = res(f, LengthLimit => 3)
F1 = orlovTruncateGeq(F, 1)
F1' = dual F1

G0s = canonicalTruncation(source F1', -2, 0)
G's = freeResolution G0s
G1's = orlovTruncateLess(dual G's, 1)
Gs = sheaf G1's

G0t = canonicalTruncation(target F1', -2, 0)
G't = freeResolution G0t
G1't = orlovTruncateLess(dual G't, 1)
Gt = sheaf G1't

phi = canonicalTruncation(F1', -2, 0)
-- c.g. https://github.com/Macaulay2/M2/issues/3865
phi = map(G0t, source phi, phi)
f = phi * G0s.cache.resolutionMap
g = G0t.cache.resolutionMap
h = f // g
f' = orlovTruncateLess(dual h, 1)

-- checking basics
assert(Gt == sheaf source f')
assert(Gs == sheaf target f')



-----
R = QQ[x, y, z] / (x*y - z^2)
M = R^1 / (x)
N = R^1 / (x,y,z)
FM = res(M, LengthLimit => 5)
FN = res(N, LengthLimit => 5)
A = map(complex N, complex M, { inducedMap(N, M) })
f = A * augmentationMap FM
g = augmentationMap FN
f' = f // g
assert(g * f' === f)
assert(source f' === FM)
assert(target f' === FN)
homotopyMap f'


-----
M = coker matrix {apply(R_*, g -> g^2)}
N = coker vars R
f = inducedMap(N, M)

end--
restart
needs "example.m2"

psi  = freeResolution(f, LengthLimit => 4)
psi' = orlovTruncateGeq(psi, 2)
phi  = Hom(F', psi')

source phi
target phi

F = source resolutionMap source phi
G = source resolutionMap target phi;

inducedMap(F, G)

prune HH target phi

HH^0 OO_X^{1}(>=0)



R = kk[x_0..x_4]/(x_0*x_1, x_2*x_3*x_4)
M = R^1/(x_0*x_2)

D = minimize singularityToDerived(M, 3, 7)

F = freeResolution(M, LengthLimit => 5)
orlovTruncateGeq(F,0)
supTruncate(F, 1)
