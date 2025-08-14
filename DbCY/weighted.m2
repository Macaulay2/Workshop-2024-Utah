-- Say X is the weighted projective space P(1,1,2).
-- Let M = S/(x,y). The sheaf on X associated to M(1) is zero (see Cox-Little-Schenck Example 5.3.11).
-- So obviously all of its twists have no cohomology. Over the associated weighted projective stack X',
-- the sheaf F associated to M(1) satisfies H^0(X’, F(i)) = k for i odd.

restart
needs "./GlobalExt.m2"

X = weightedProjectiveSpace {1,1,2}
S = ring X
M = S^{1}/(x_0,x_1)

assert all(3, i -> 0 == HH^i(X, sheaf M))
assert(0 == globalExt(complex S^1, complex M))
assert(0 == globalExt(sheaf S^1, sheaf M))
assert all(3, i -> 0 == globalExt(i, complex S^1, complex M))
assert all(3, i -> 0 == globalExt(i, sheaf S^1, sheaf M))

-- FIXME: NormalToricVarieties should prune this to zero:
0 == prune sheaf M

-- TODO: NormalToricVarieties doesn't support this yet
0 == hilbertPolynomial_X M
apply(0..10, i -> hilbertFunction(2*i, M))

-- manual test
X = weightedProjectiveSpace {1,1,2}
S = ring X
M = S^{1}/(x_0,x_1)
minimize part_0 Hom(freeResolution truncate_1 S^1, complex M)
