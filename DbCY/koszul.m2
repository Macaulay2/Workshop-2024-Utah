restart
needs "./GlobalExt.m2"
S = ZZ/32003[x,y,z]/(x^2,x*y)
C = koszulComplex matrix{{x,y}}
prune sheaf HH C
globalExt(0, C, C)
ExtTable(Proj S, {C})

----

S = ZZ/32003[x,y,z,w]
C = naiveTruncation(koszulComplex vars S, (1,2))
prune HH C
globalExt(1, C, C)
ExtTable(Proj S, {C})
