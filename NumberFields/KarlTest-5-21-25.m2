restart
uninstallPackage "NumberFields"
loadPackage "NumberFields"
R = QQ[t,i]/ideal(t^3-2,i^2+1)
numberField(R, Verbose=>true, usePari => false)

A = QQ[x]
f = x^3-2
splittingField (f, Verbose=>true)
B = QQ[y]
splittingField(y^4+y^3+y^2+y+1)

restart
loadPackage "NumberFields"
K = numberField(QQ[i]/(i^2+1), Verbose=>true)
A = K[x]
f = x^3-2
splittingField (f, Verbose=>true, UsePari =>false)

restart
loadPackage "NumberFields"
R = QQ[x];
f = x^2-2;
splittingField f
f = x^3-7;
splittingField f

restart
loadPackage "NumberFields"
needsPackage "PushForward"
k = QQ[x]/ideal(x^3-2)
L = k[y]/ideal(y^2-2)
nf = numberField L


S2 = (flattenRing(L))#0
describe S2
(myMod, myGens, myFun) = pushFwd(map(S2, coefficientRing S2))
prune myMod
A4 = (coefficientRing S2)[gens S2]
S4 = A4/sub(ideal(S2), A4)
S4 = newRing(S2, MonomialOrder=>GRevLex, Degrees=>{1,1})
(myMod, myGens, myFun) = pushFwd(map(S4, coefficientRing S4))

S3 = QQ[v,u, Degrees=>{{0,1},{1,0}}]/ideal(v^3-2,u^2-2)
(myMod3, myGens3, myFun3) = pushFwd(map(S3, coefficientRing S3))

nf = numberField L


restart
uninstallPackage "NumberFields"
loadPackage "NumberFields"
installPackage "NumberFields"
check NumberFields