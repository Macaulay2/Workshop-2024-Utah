loadPackage ("NumberFields", Reload=>true)

R = QQ[]
S = QQ[x]/(x^6+x^5+x^4+x^3+x^2+x+1)
S = QQ[x]/(x^3-2)
S = QQ[x,y]/(x^3-2, y^2+y+1)

phi = map(S, R)
myList = pushFwd phi
elementOverR = myList#2

myMatrix = elementOverR(x^2)

A = map(S^1, S^1, {{x^3}}) -- ^1 turns a ring into a rank-one module, ^2 turns it into a rank-two module, etc.
pushFwd(map(S^1, S^1, {{x^2}}))
trace pushFwd(map(S^1, S^1, {{1_S}}))

degree(myList#0)