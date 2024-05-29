-- Note: to load the necessary packages from the repository, either:
-- first call installPackage in testbot.m2:
--  installPackage("Elimination", FileName => "/home/macaulay/Elimination.m2")
-- or specify the path to needsPackage:
--  needsPackage("Elimination", FileName => "/home/macaulay/Elimination.m2")
needsPackage "Elimination"

-- test code and assertions
R = ZZ/101[a..d]
I = monomialCurveIdeal(R, {1, 3, 4})
assert(eliminate(I, {b}) == ideal(c^4-a*d^3))
