restart
loadPackage "NumberFields"

R = QQ[x]/ideal(x^3-2)
S = R[y]/ideal(y^2+y+1)
remakeField S
U = (flattenRing)
T = QQ[x,y, Degrees=>{0,0}]/ideal(x^3-2,y^2+y+1)

K = numberField R
L = numberField S

degree L
ring L

phi = map(T, QQ[], {})
pushFwd(phi)

loadPackage("PushForward", DebuggingMode => true, Reload=>true)



R = QQ[x,y]/ideal(y^2-x*(x-1)*(x-2), x^2-y*(y-1)*(y-2))