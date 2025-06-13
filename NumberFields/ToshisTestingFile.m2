restart
loadPackage ("NumberFields", Reload=>true)
K = QQ[x]/ideal(x^3-2)

L = K[y]/ideal(y^2-2)
NF = numberField L 


-- Discuss with Karl