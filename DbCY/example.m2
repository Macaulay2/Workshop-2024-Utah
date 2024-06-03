needs "DbCY.m2"

R = QQ[x_0,x_1,x_2] / sum(3, i -> x_i^3)
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
