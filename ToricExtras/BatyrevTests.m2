TEST ///
-- we know this construction works
V = Batyrev({1,1,1,1,1}, {}, {0})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V == 3)
///