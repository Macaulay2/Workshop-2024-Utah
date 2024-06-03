TEST ///
-- a 2-dimensional fan where the partition sizes of the primitive collections are equal
V = Batyrev({1,1,1,1,1}, {}, {})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V == 3)
///

TEST ///
-- a 3-dimensional fan where the partition sizes of the primitive collections are not the same
V = Batyrev({1,1,1,2,1},{0,0},{})
assert(isWellDefined V)
assert(isSmooth V)
assert(isProjective V)
assert(rank picardGroup V ==3)
///

-----------------------------------------
------------- SCRATCH SPACE -------------
-----------------------------------------

restart
needsPackage "NormalToricVarieties";
load "Desktop/Projects/Workshop-2024-Utah/ToricExtras/BatyrevConstructions.m2";

allGood = (V) -> (
    i := isWellDefined V;
    j := isSmooth V;
    k := isProjective V;
    l := rank picardGroup V == 3;
    i and j and k and l
    )




