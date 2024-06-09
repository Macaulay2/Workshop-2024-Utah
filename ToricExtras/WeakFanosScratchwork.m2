--- just some stuff I was playing around with

-- Goal: product all examples of smooth weak Fano 3-folds of Picard number 3 w/ 5 primitive collections
-- After this, generalize to higher dimensions

-- needsPackage "NormalToricVarieties"
-- loadPackage "ToricExtras"

combos = {{2,1,1,1,1}, {1,2,1,1,1}, {1,1,2,1,1}, {1,1,1,2,1}, {1,1,1,1,2}}

-- Complete list of smooth toric Fano 3-folds of Picard rank 3 (w/ 5 primitive collections)

V26 = batyrevConstructor({2,1,1,1,1}, {0}, {});
V29 = batyrevConstructor({1,1,2,1,1}, {1}, {});
V26p = batyrevConstructor({1,1,2,1,1}, {0}, {}); -- isomorphic to V26 above

-- (These are all the possible inputs that generate Fano 3-folds)
-- We have one "redundancy"


-- Smooth toric weak Fano ( \ Fano) 3-folds of Picard rank 3 (w/ 5 primitive collections)

V1 = batyrevConstructor({2,1,1,1,1}, {2}, {}); -- deg = 58

V2 = batyrevConstructor({1,2,1,1,1}, {0}, {}); -- 48
V3 = batyrevConstructor({1,2,1,1,1}, {1}, {}); -- 54
V4 = batyrevConstructor({1,2,1,1,1}, {2}, {}); -- 64

V5 = batyrevConstructor({1,1,2,1,1}, {0}, {1}); -- 46
V6 = batyrevConstructor({1,1,2,1,1}, {1}, {0}); -- 46

V7 = batyrevConstructor({1,1,1,2,1}, {0,0}, {}); -- 46

V8 = batyrevConstructor({1,1,1,1,2}, {0}, {}); -- 46
V9 = batyrevConstructor({1,1,1,1,2}, {1}, {}); -- 48

l = {V1, V2, V3, V4, V5, V6, V7, V8, V9}

for X in l do (
    << degreeOfThreefold(X) << " " << isSimplicial X << endl;
);

-- Hodge number h^{1,2} (compute via HH^2(X, cotangentSheaf X)) not useful, seems to be always zero in toric case (why? some vanishing theorem?)

-- Are any entries in this list redundant? One idea (laborious) is to find explicit isomorphisms, since we have the rays and cones. But is there a better way? Other invariants besides degree to possible differentiate?
-- Question: is this the complete list?

V = V1
isFano(V)
isNef(-toricDivisor V)

degreeOfThreefold = X -> (
    K := toricDivisor X;
    c := chern(1, OO (-K));
    integral (c*c*c)
)


degreeOfThreefold(V)

-- Moreover, are V2 and V4 isomorphic?

for b from 0 to 10 list (
    V = batyrevConstructor({2,1,1,1,1}, {b}, {});
    if isWellDefined V and isSmooth V and isProjective V then (
	if (not isFano V) and (isNef(-toricDivisor V)) then (
	    << "{2,1,1,1,1}, {" << b << "}, {}" << endl;
	)
    )
);

for b from 0 to 10 list (
    V = batyrevConstructor({1,2,1,1,1}, {b}, {});
    if isWellDefined V and isSmooth V and isProjective V then (
	if (not isFano V) and (isNef(-toricDivisor V)) then (
	    << "{1,2,1,1,1}, {" << b << "}, {}" << endl;
	)
    )
);

for b from 0 to 10 list (
    for c from 0 to 10 list (
        V = batyrevConstructor({1,1,2,1,1}, {b}, {c});
        if isWellDefined V and isSmooth V and isProjective V then (
	    if (not isFano V) and (isNef(-toricDivisor V)) then (
	        << "{1,1,2,1,1}, {" << b << "}, {" << c << "}" << endl;
	    )
        )
    )
);

for b1 from 0 to 10 list (
    for b2 from 0 to 10 list (
        V = batyrevConstructor({1,1,1,2,1}, {b1, b2}, {});
        if isWellDefined V and isSmooth V and isProjective V then (
	    if (not isFano V) and (isNef(-toricDivisor V)) then (
	        << "{1,1,1,2,1}, {" << b1 << ", " << b2 << "}, {}" << endl;
	    )
        )
    )
);

for b from 0 to 10 list (
    V = batyrevConstructor({1,1,1,1,2}, {b}, {});
    if isWellDefined V and isSmooth V and isProjective V then (
	if (not isFano V) and (isNef(-toricDivisor V)) then (
	    << "{1,1,1,1,2}, {" << b << "}" << endl;
	)
    )
);


----


PP1 = toricProjectiveSpace 1
FF2 = hirzebruchSurface 2
prod = cartesianProduct (PP1, FF2)

degreeOfThreefold(prod)

-- write some code to compute primitive collections given a normalToricVariety (just loop over all subsets is a viable brute force approach)

primitiveCollections = X -> (
    n = #(rays X);
    maxConez = max X;
    isCone = c -> (for maxCone in maxConez do 
	(if isSubset(c, maxCone) then return true);
	return false;
    );
    pcollections = {};
    for subset in subsets(n) do (
	isValid = true;
	if isCone(subset) then (isValid = false)
	else for i from 0 to #subset-1 do (
	    a = subset_{0..i-1};
	    b = subset_{i+1..#subset-1};
	    trycone = join(a, b);
	    if not isCone(trycone) then (isValid = false; break)  
	);
    	if isValid then pcollections = append(pcollections, subset);
    );
    return pcollections;
);

for X in l do (<< primitiveCollections(X) << endl);


primitiveCollections(batyrevConstructor({1,1,1,1,1}, {0}, {})) -- 5 pc
primitiveCollections(prod) -- 3 pc
