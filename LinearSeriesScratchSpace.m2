restart
needsPackage "NormalToricVarieties"

P2 = toricProjectiveSpace 2

S = ring P2

basis (2,S)

m = flatten entries basis (3,S)
m#0
S = ring m#0

peek S

S.variety
B = ideal(S.variety)
I = radical ideal(m)
isSubset(B,I)



-- rays(T.variety)
-- T.variety.max