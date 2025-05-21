restart
loadPackage "NumberFields"
R = QQ[t,i]/ideal(t^3-2,i^2+1)
numberField(R, Verbose=>true, usePari => false)

A = QQ[x]
f = x^3-2
splittingField (f, Verbose=>true)
B = QQ[y]
splittingField(y^4+y^3+y^2+y+1)
