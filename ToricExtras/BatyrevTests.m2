TEST ///
-- a 2-dimensional fan where the partition sizes of the primitive collections are equal
V = batyrevConstructor({1,1,1,1,1}, {0}, {})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V == 3)
///

TEST ///
-- a 3-dimensional fan where the partition sizes of the primitive collections are not the same
V = batyrevConstructor({1,1,1,2,1},{0,0},{})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V ==3)
///






