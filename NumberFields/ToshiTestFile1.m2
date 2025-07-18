restart
loadPackage ("NumberFields", Reload=>true)
-- Compositums works for this simple case...

K= numberField (QQ[x]/(x^4+1)) 
R = K[z]
f = z^2+z+1 
K= numberField (QQ[x]/((x^4+1)))
splittingFieldPari(f)

-- Discuss with Karl about splitting field and no applicable strategy for minimalPrimes, ideal...

-- Testing for 
(gens P2)_0


-- GP FUNCTIONALITY BELOW
K = nfinit(y)
L = nfcompositum(K, x^3 - 2, x^2-2, 1)
L=nfcompositum(nfinit(y),x^2+1, x^3-2, 1)
for(d=0, poldegree(f), print(polcoeff(f,d),,))