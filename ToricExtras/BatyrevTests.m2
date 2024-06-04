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
assert(rank picardGroup V == 3)
///

TEST ///
-- a 12-dimensional fan where the partition sizes of the primitive collection are all different
-- 15 rays with 105 maximal cones
V = batyrevConstructor({1,2,3,4,5},{0,0,0,0},{0,0})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V == 3)
///






