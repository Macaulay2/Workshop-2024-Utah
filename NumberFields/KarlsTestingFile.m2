restart
loadPackage "NumberFields"

R = QQ[x]/ideal(x^3-2)
S = R[y]/ideal(y^2+y+1)
remakeField S
U = (flattenRing S)#0
pushFwd(map(U, QQ))
T = QQ[x,y, Degrees=>{0,0}]/ideal(x^3-2,y^2+y+1)

K = numberField R
L = numberField S

degree L
ring L

phi = map(T, QQ[], {})
pushFwd(phi)

loadPackage("PushForward", DebuggingMode => true, Reload=>true)



R = QQ[x,y]/ideal(y^2-x*(x-1)*(x-2), x^2-y*(y-1)*(y-2))
psi = map(R, QQ)
pushFwd(psi)



-----------------
--Karl testing isFieldAutomorphism
loadPackage "NumberFields"
R = QQ[t, z]/ideal(t^3-2, z^2+z*t+t^2)
K = numberField R
B = basis K
if ((B#1) == (gens(K))#0) and ((B#4)==(gens(K))#1) then (
    M=matrix{{1,0,0,0,0,0}, {0,0,0,0,1,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,1,0,0,0,0},{0,0,0,0,0,1}}**QQ;
    isFieldAutomorphism(K, M) -- this should be false
)

((gens(K))#1)^3
ringMapFromMatrix(K,M)
a1 = (gens(K))#0
a2 = (gens(K))#1
basis K
N = vector(1_(ring K), K) | vector(a2, K) | vector(a1*a2,K) | vector(a1^2*a2, K) | vector(a1, K) | vector(a1^2, K)
isFieldAutomorphism(K,N)