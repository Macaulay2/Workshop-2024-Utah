restart
loadPackage "NumberFields"

R = QQ[x]/ideal(x^3-2)
S = R[y]/ideal(y^2+y+1)
remakeField S
T = QQ[x,y, Degrees=>{0,0}]/ideal(x^3-2,y^2+y+1)

K = numberField R
L = numberField S

degree L
ring L

phi = map(T, QQ[], {})
pushFwd(phi)

loadPackage("PushForward", DebuggingMode => true, Reload=>true)

