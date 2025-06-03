-- package is moved to OrlovFunctors.m2

restart
needsPackage "OrlovFunctors"
R = ZZ/101[x_0..x_4] / ideal(x_0*x_1, x_2*x_3*x_4)
X = Proj R
M = coker matrix{{x_0*x_2}}
D = minimize singularityToDerived(M, 3, 7)
D = sheaf D
errorDepth=2
debugLevel=2
C = OO_X^1
E = OO_X^{1}
D
C ++ C
elapsedTime RHom(C, D) -- 16.8s
elapsedTime RHom(C, D, 0) -- 92.2s

prune HH D
D.dd
elapsedTime prune HH RHom(D, OO_X^1) -- 22s -> 16s, ranks 113 103 69 7
elapsedTime RHom(OO_X^1, D)
debug Varieties

D = minimize singularityToModules(M, 3, 7)
D = sheaf D
errorDepth=2
debugLevel=2
RHom(OO_X^1, D)

prune HH D
D.dd
elapsedTime prune HH RHom(D, OO_X^1) -- 22s -> 16s, ranks 113 103 69 7
elapsedTime RHom(OO_X^1, D)
debug Varieties

restart
needsPackage "OrlovFunctors"
S = ZZ/101[x_0..x_4]
f = sum for i from 0 to 4 list x_i^5
R = S / ideal(f)
M = coker vars R
i = 0;
j = 7;
F = singularityToDerived(M, i, j)
dual((res M)[4]) -- this is RHom(k, R)[4], which is isomorphic to k(a) = k by the Gorenstein property
--This implies that \widetilde{F} = O_X[-3].
G = res M
K = coker G.dd_5
singularityToDerived(K, 0, 9)
prune ker oo.dd_(-1)
F = singularityToDerived(M, i, 6)
G = singularityToDerived(K, 0, 10)
G[-4] == F
F.dd_3
(G[-4]).dd_3
for i from -1 to 3 do print HH_i(o8) == 0

prune HH_0(o8) == 0
prune HH_(-1)(o8) == 0
prune HH_1(o8) == 0
prune HH_2(o8) == 0
prune HH_3(o8) == 0


--DEMO
--Michael Brown, Souvik Dey, Geoffrey Fatin, Alicia Lamarche,
--Guanyu Li, Mahrud Sayrafi, Tim Tribone, and Rachel Webb
restart
needsPackage "OrlovFunctors"

S = ZZ/101[x_0..x_4]
f = sum for i from 0 to 4 list x_i^5
R = S / ideal(f) --affine cone of the quintic
kk = coker vars R --the residue field of R
F = singularityToDerived(kk, 0, 7)
--This is \Phi(kk) before sheafifying.
--The complex has infinite length; it is unbounded on the right.
--The tail is a matrix factorization.
X = Proj(R)
Ftilde = sheaf F
--This is \Phi(kk).
--What is the homology of this complex?
naiveTruncation(prune HH Ftilde, -1, 3)
--Thus, \Phi(kk) = O_X[-3].
--(We're getting the shift by -3 because dim X = 3.)

--Future work: implement the functor going the other direction.
