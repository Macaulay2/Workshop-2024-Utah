restart
loadPackage "NumberFields"
R = QQ[t,i]/ideal(t^3-2,i^2+1)
numberField(R, Verbose=>true, usePari => false)
