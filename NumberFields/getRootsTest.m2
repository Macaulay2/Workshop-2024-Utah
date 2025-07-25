restart
loadPackage "NumberFields"
NF = numberField( QQ[x]/(x^2-2))
R1 = NF[u]
p = minimalPolynomial primitive NF --This gives minimal polynomial over \mathbb{Q}[x]
M0 = map(R1,ring p,{(gens R1)_0}); --We use this to convert this min poly into one over NF[u] (relabel x into u)
M0(p)
getRoots(M0(p))
