TEST ///
-- a 2-dimensional fan where the partition sizes of the primitive collections are equal
-- This constructs the unique smooth toric del Pezzo surface of Picard rank 3, which is the Bl_1(P^1 x P^1) = Bl_2(P^2)
V = batyrevConstructor({1,1,1,1,1}, {0}, {})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V == 3)
assert(isFano(V))
///

TEST ///
-- This constructs a smooth toric weak del Pezzo surface of Picard rank 3, which is the blowup of Bl_1 P^2 at the torus fixed-point on the exceptional divisor
V = batyrevConstructor({1,1,1,1,1}, {1}, {})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V == 3)
assert(not isFano(V)) -- not Fano
assert(isNef(-toricDivisor V)) -- but it is weak Fano
///

TEST ///
-- This constructs a smooth toric Fano of Picard rank 3 which is 3-26 on Fanography. It can be described as the blowup of P^3 at the disjoint union of a point and a line
V = batyrevConstructor({2,1,1,1,1}, {0}, {})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V == 3)
assert(isFano(V))

degreeOfThreefold = X -> (
    K := toricDivisor X;
    c := chern(1, OO (-K));
    integral (c*c*c)
)

assert(degreeOfThreefold(V) == 46)
///

TEST ///
-- This constructs a smooth toric Fano of Picard rank 3 which is 3-29 on Fanography. It can be described as the blowup of Bl_1 P^3 along a line in the exceptional divisor
V = batyrevConstructor({2,1,1,1,1}, {1}, {})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V == 3)
assert(isFano(V))

degreeOfThreefold = X -> (
    K := toricDivisor X;
    c := chern(1, OO (-K));
    integral (c*c*c)
)

assert(degreeOfThreefold(V) == 50)
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
-- a 4-dimensional fan where the partition sizes of the primitive collections vary
V = batyrevConstructor({1,2,1,2,1},{0,0},{})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V == 3)
///

-- GOOD TEST BUT TAKES TOO LONG :(
--TEST ///
-- a 12-dimensional fan where the partition sizes of the primitive collection are all different
-- 15 rays with 105 maximal cones
--V = batyrevConstructor({1,2,3,4,5},{0,0,0,0},{0,0})
--assert(isWellDefined V)
--assert(isSmooth V)
--assert(isProjective V)
--assert(rank picardGroup V == 3)
--///






