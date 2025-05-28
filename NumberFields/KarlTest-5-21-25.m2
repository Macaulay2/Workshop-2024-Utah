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
uninstallPackage "NumberFields"
loadPackage "NumberFields"
installPackage "NumberFields"