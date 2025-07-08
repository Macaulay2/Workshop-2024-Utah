restart
loadPackage ("NumberFields", Reload=>true)
-- Compositums works for this simple case...
K0 = numberField (QQ[x]/ideal(x^3-2))
L = numberField (QQ[y]/ideal(y^2+1))
M = compositumPari (K0,L)
-- We expect the following to be two as M_1(x) should have "same behavior" as x in the compositum field.
M_1(x)^3 
-- ==============================Discuss this with Karl...=========================
-- K0 = numberField (QQ[x]/ideal(x^3-2))
-- K1 = numberField (K0[z]/(z^2+z+1))
-- L = numberField (QQ[y]/ideal(y^2+1))
M = compositumPari(K1, L)
-- ================================================================================


-- Discuss with Karl

-- Testing for 
(gens P2)_0


-- GP FUNCTIONALITY BELOW
K = nfinit(y)
L = nfcompositum(K, x^3 - 2, x^2-2, 1)
L=nfcompositum(nfinit(y),x^2+1, x^3-2, 1)
for(d=0, poldegree(f), print(polcoeff(f,d),,))