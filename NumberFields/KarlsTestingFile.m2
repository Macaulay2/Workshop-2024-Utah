restart

uninstallPackage "NumberFields"
loadPackage "NumberFields"
installPackage "NumberFields"
check NumberFields

restart
loadPackage "NumberFields"
debugLevel = 2
R = QQ[a,b]/ideal(a^3-2, b^2+b+1)
simpleExtension(R)
R = QQ[a,b,c]/ideal(a^2+1, b^4+b^3+b^2+b+1,c^5-2)
simpleExtension(R)

restart
loadPackage "NumberFields"

R = QQ[x]
f = x^3-2
K = (splittingField(f, Verbose=>true))
use R
g = x^5+x^4+x+3
K = (splittingField(g, Verbose=>true))

h = x^5-10*x+2
K = (splittingField(h, Verbose=>true))
h1 = x^5-10*x+2

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

restart
loadPackage("PushForward", DebuggingMode => true, Reload=>true)
loadPackage "NumberFields"
R = QQ[u,v]/ideal(u^3-2, v^2+v+1)
K = numberField R



R = QQ[x,y]/ideal(y^2-x*(x-1)*(x-2), x^2-y*(y-1)*(y-2))
psi = map(R, QQ)
pushFwd(psi)


break
restart

loadPackage "NumberFields"
S = QQ[a,b,c]/ideal(a^3-2, a^2+a*b+b^2, sum(apply(11, t->c^t)))
Ss = S[x]
f = x^3-2
time getRoots(f, Strategy=>decompose)
T = QQ[a,b]/ideal(a^3-2, a^2+a*b+b^2)
psi = map(S, T)
R = time numberField(S)
R2 = time numberField(T)
minimalPolynomial( sum gens S)
minimalPolynomial( (gens R)#0)
minimalPolynomial( sum gens R)
minimalPolynomial( c, psi)

time kappa = simpleExt(R2);
time inverse(kappa#1);
time inverseNumberFieldAutomorphism(kappa#1)

time simpleExt(R);
time simpleExt(R, Strategy=>kernel);
time omega = simpleExt(R, Strategy=>null);

elapsedTime inverse (omega#1);
time inverseNumberFieldAutomorphism(omega#1);

S = QQ[x]/ideal(x^2+1)
T = toField(QQ[y]/ideal(y^4+1))
psi = map(T, S, {y^2})
minimalPolynomial(y, psi)


-----------------
--Karl testing isFieldAutomorphism
loadPackage "NumberFields"
R = QQ[t, z]/ideal(t^3-2, z^2+z*t+t^2)
K = numberField R
B = basis K
a1 = (gens(K))#0
a2 = (gens(K))#1
phi = map(ring K, ring K, {a2,a1})
N = matrixFromRingMap(K, K, phi)
ringMapFromMatrix(K, N)
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
time isFieldAutomorphism(K,N)

restart
loadPackage "NumberFields"
    R = numberField(QQ[a]/ideal(a^4+a^3+a^2+a+1))
    b = (gens(R))#0
    h3 = map(R, R, {b^3})
    assert(isWellDefined h3)
    time inverse h3
    time inverseNumberFieldAutomorphism(h3)
