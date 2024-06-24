--- just some stuff I was playing around with --Rohan

-- Goal: produce all examples of smooth weak Fano 3-folds of Picard number 3
-- After this, generalize to higher dimensions

-- needsPackage "NormalToricVarieties"
-- loadPackage "ToricExtras"

-- sometimes I may need to do this:
-- path = append(path, "~/Documents/PhD/Workshop-2024-Utah/"); stack path


-- "sameToricVariety" function written by the Projective Bundle group
-- just a way to see if two toric varieties are the same by a simple permutation
-- sufficient, but not necessary for isomorphism

sameToricVariety  = (X, Y) -> X === Y or fan X == fan Y or any(permutations(d := dim X), perm -> sameToricVariety'(X, Y, id_(ZZ^d)_perm))
sameToricVariety' = (X, Y, m) -> try (f := inducedMap map(X, Y, m)) then f ideal X == ideal Y else false

-- We will generate weak Fano varieties of Picard rank 3 with 5 primitive collections using the Batyrev constructor:

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

-- The following inputs produce Fano outputs (modifying the code above, to check for Fano instead of weak Fano):

-- Complete list of smooth toric Fano 3-folds of Picard rank 3 w/ 5 primitive collections

V26 = batyrevConstructor({2,1,1,1,1}, {0}, {}); -- 3-26 on Fanography
V29 = batyrevConstructor({2,1,1,1,1}, {1}, {}); -- 3-29 on Fanography
V26p = batyrevConstructor({1,1,2,1,1}, {0}, {0}); -- isomorphic to V26 above

sameToricVariety(V26, V26p) -- true!

-- These are all the possible inputs that generate Fano 3-folds
-- We have exactly one "redundancy"


----
-- Smooth toric weak Fano ( \ Fano) 3-folds of Picard rank 3 (w/ 5 primitive collections)

-- We get the following inputs producing weak Fano \ Fano outputs:

V1 = batyrevConstructor({2,1,1,1,1}, {2}, {}); -- deg = 58, inv = 38

V2 = batyrevConstructor({1,2,1,1,1}, {0}, {}); -- 48, 28
V3 = batyrevConstructor({1,2,1,1,1}, {1}, {}); -- 54, 34
V4 = batyrevConstructor({1,2,1,1,1}, {2}, {}); -- 64, 45

V5 = batyrevConstructor({1,1,2,1,1}, {0}, {1}); -- 46, 26
V6 = batyrevConstructor({1,1,2,1,1}, {1}, {0}); -- 46, 29

V7 = batyrevConstructor({1,1,1,2,1}, {0,0}, {}); -- 46, 26

V8 = batyrevConstructor({1,1,1,1,2}, {0}, {}); -- 46, 26
V9 = batyrevConstructor({1,1,1,1,2}, {1}, {}); -- 48, 31

l = {V1, V2, V3, V4, V5, V6, V7, V8, V9}

-- We will use the following invariants to differentiate: 
-- degree of a 3 fold: defined as (-K)^3
-- the cohomology group h^0(Omega^1(-K))
-- I tried different variations of the above cohomology group, this one seems to be the best
-- It also has a nice geometrical interpretation: Omega^1(-K) is T_X, the tangent sheaf. So it is the size of the space of global vector fields on X

degreeOfThreefold = X -> (
    K := toricDivisor X;
    c := chern(1, OO (-K));
    integral (c*c*c)
)

genusOfThreefold = X -> 1 + (degreeOfThreefold X)/2

randomInvariant = X -> HH^0(X, Omega1D(X, -toricDivisor X))

for X in l do (
    << degreeOfThreefold(X) << " " << randomInvariant(X) << endl;
);

sameToricVariety(V2, V9) -- false
sameToricVariety(V5, V6) -- false
sameToricVariety(V5, V7) -- false
sameToricVariety(V5, V8) -- false
sameToricVariety(V6, V7) -- false
sameToricVariety(V6, V8) -- false
sameToricVariety(V7, V8) -- true!
-- V7 and V8 are isomorphic, by the simple linear transformation (x,y,z) --> (z, y, x). :)

-- Thus, the only remaining question is, are V5 and V7 isomorphic? I cannot find a numerical invariant to distinguish them. So I will use a direct argument among relations in the fan. 

-- V5 and V7 are non-isomorphic, because V7 has two opposite rays ({0,0,-1} and {0,0,1}) while V5 does not. (test this with rays V5 and rays V7)

-- Thus we have 9-1 = 8 distinct weak Fano \ Fano varieties! This the complete list, we don't have to worry about trying larger values of bs and cs. This can be proved by Proposition 4 in [1], characteristic weak Fano torics (it basically says the input b / c can never be more than 2, for one. So going up to 10 is way overkill!


----
--misc


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





------------

--- 3 primitive collections

-- These are constructed in the following way.

-- Projective bundles of rank 2 vector bundles over smooth toric weak Fanos of Picard rank 2
-- There are only 3 smooth projective toric weak Fanos of Picard rank 2
-- P1 x P1, and Bl_1 P^2 and F_2 (i.e. Hirzebruch surfaces F_0, F_1, and F_2). Of these F_0 and F_1 are Fano.
-- How can we specify projective bundles on such varieties in a way as to reduce redundancy?

--  We are constructing threefolds, so it is a projectivization of a rank 2 vector bundle, which is given to be decomposable. So it is a direct sum of two line bundles. Since these are toric varieties, each line bundle can be written as O(D) where D is a toric divisor.
-- Since the base is Picard rank 2, we can write any toric divisors in terms of two generators. We will use the first two rays: on a Hirzebruch surface, these are always the ones that correspond to the rays (1,0) and (0,1). Our coefficients for the other rays will be zero. Call these d1 and d2. Call the coefficients c_1 and c_2. So our divisor is c_1d_1 + c_2d_2. 
-- P(O(D1) + O(D2)) is isomorphic to P(O + O(D2-D1)). So we can let one divisor be zero by default, and we only need to pick two numbers to be coefficients for the second divisor.
-- To further reduce redundancy, note that P(O + O(D)) is isomorphic to P(O + O(-D)). So we may flip D. Let us then flip D so that the first coefficient c1 is nonnegative.
-- To further reduce redundancy, note that if c1 is 0, then we can flip it so that c2 is nonnegative
-- Finally, in the case that the base is P1 x P1 (Hirzebruch surface F_0), we have a symmetry exchanging c1 and c2. So in this case we may insist that c1 >= c2 


-- Important fact to prove: if a projective bundle X --> B is weak Fano, B is weak Fano (probably true, idk how to prove it yet). Furthermore, if X is Fano, B is Fano. 

projectiveBundleConstructor = (a, l) -> (
    assert(0 <= a and a <= 2);
    assert(length l == 2);
    X = hirzebruchSurface a;
    zeroD = toricDivisor({0,0,0,0}, X);
    D = toricDivisor(join(l, {0,0}), X);
    projectivizationOfBundle({zeroD, D})
)

projectiveBundleConstructor(0, {0,0})
PP1 = toricProjectiveSpace 1;

sameToricVariety(projectiveBundleConstructor(0, {0,0}), cartesianProduct((PP1, PP1, PP1)))


-- Generating all Fano 3-folds that have 3 primitive collections, as projective bundles

-- As described by Fanography
-- Code	   | genus    | possible input	   | description      	 
-- 3-31	     27	        0, {1,1} 	
-- 3-30	     26	        1, {1,1}    	    1,1, because l ~ D1 + D2 in Pic
-- 3-28	     25	        1, {0,0}
-- 3-27	     25	        0, {0,0}
-- 3-25	     23	        0, {1,-1}    	    because we can shift O(1,0)+O(0,1) to O(0,0)+O(1,-1))


for a from 0 to 1 do (
    for c1 from 0 to 10 do (
	for c2 from -10 to 10 do (
	    if not ((c1 == 0 and c2 < 0) or (a == 0 and (c2 > c1))) then (
                X = projectiveBundleConstructor(a, {c1, c2});
	    	if (isFano X) then (
	    	    << "FF" << a << ". D = " << c1 << " d1 + " << c2 << " d2";
	            << "\t degree = " << degreeOfThreefold(X);
		    << "\tgenus = " << genusOfThreefold(X);
		    <<"\tinvariant = " << randomInvariant(X) << endl;
		)
	    )
       )
   )
);

-- The only redundancy generated by this code is it generates 3-28 twice (F_1 x P^1)
-- Once as     projectiveBundleConstructor(0, {1, 0}) -- to see this is F_1 x P^1, look at list of rays (and cones)
-- and once as projectiveBundleConstructor(1, {0, 0})

-- Question: is projectiveBundleConstructor(0, {r, 0}) isomorphic to F_r x P^1 ??



------
-- Now, onto weak Fanos!


counter = 0;
for a from 0 to 2 do (
    for c1 from 0 to 10 do (
	for c2 from -10 to 10 do (
	    if not ((c1 == 0 and c2 < 0) or (a == 0 and (c2 > c1))) then (
	    	X = projectiveBundleConstructor(a, {c1, c2});
	    	if (not isFano X and isNef(-toricDivisor X)) then (
		    counter = counter + 1;
		    << counter << ". ";
	    	    << "FF" << a << ". D = " << c1 << "d1 + " << c2 << "d2";
	            << "\t degree = " << degreeOfThreefold(X);
		    << "\tgenus = " << genusOfThreefold(X);
		    << "\t inv = " << randomInvariant(X) << endl;
		)
	    )
       )
    )
);

-- Outputs 15 varieties, with 8 distinct genuses: 17, 21, 24, 25, 27, 28, 29, 33
-- I.e. degrees 32, 40, 46, 48, 52, 54, 56, 64

-- degree 32: (0, {2,-2})
-- degree 40: unique since (0, {1,-2}) and (0, {2,-1}) are obv isomorphic
-- degree 46: (1, {0,1})
-- degree 48: inv = 30: (0, {2,0}), (2, {0,0})
--    	      inv = 28: (1, {1,0})
--    	      inv = 31: (1, {1,2})
-- degree 52: (2, {2,1})
-- degree 54: (1, {2,1})
-- degree 56: (0, {2,1}), (1, {2,2})
-- degree 64: inv = 45: (0, {2,2}), (1, {3,2})
--    	      inv = 47: (2, {4,2})

-- The only open questions are: 
-- are the two degree 48's with inv=30 the same?
-- are the two degree 56's the same?
-- are the two degree 64's with inv=45 the same?

-- sameToricVariety did not help here

V48v1 = projectiveBundleConstructor(0, {2,0})
V48v2 = projectiveBundleConstructor(2, {0,0}) -- this is F_2 x P^1
randomInvariantP(V48v1)
randomInvariantP(V48v2)

V56v1 = projectiveBundleConstructor(0, {2,1})
V56v2 = projectiveBundleConstructor(1, {2,2})
sameToricVariety(V56v1, V56v2)
randomInvariantP(V56v1)
randomInvariantP(V56v2)

V64v1 = projectiveBundleConstructor(0, {2,2})
V64v2 = projectiveBundleConstructor(1, {3,2})
randomInvariantP(V64v1)
randomInvariantP(V64v2)

-- So there are somewhere between 12 to 15 varieties here


--- 

-- According to paper of Sato [1], the "weakened" Fano-3 folds (a subset of weak Fano \ Fano) of Picard rank 3 are precisely F_2 x P^1, and a variety he calls X^0_3 (which we are told has degree 52

-- F_2 x P^1 is generated by (2, {0,0})
-- The degree 52 surface by above can only generated exactly one way, by (2, {2,1}).


-- Some notes on Sato's degree 52 surface, X^0_3. 

-- It is a blowup of P_P^1( O + O + O(2)).
-- We can construct the latter variety via

D = toricDivisor({2,0}, PP1)
zeroD = toricDivisor({0,0}, PP1)
satoBase = projectivizationOfBundle({zeroD, zeroD, D})

rays satoBase

X03 = projectiveBundleConstructor(2, {2,1})
rays X03

primitiveCollections(satoBase)
primitiveCollections(X03)

-- So there are certainly weak Fanos in this list that are not "weakened Fanos"

-- [1] Sato, The classification of smooth toric weakened Fano 3-folds	
	
	
	
	
	
	
	
	
-- Trying to define some new invariants now

X = VV48v3

Omega1D = (X, D) -> (D = -D; rank1 = OO D; omega = cotangentSheaf(X); prune sheafHom(rank1, omega))

HH^3(X, Omega1D(X, toricDivisor X))

randomInvariant1 = X -> HH^3(X, Omega1D(X, toricDivisor X))
randomInvariant2 = X -> HH^1(X, Omega1D(X, -toricDivisor X))

vectorFields = X -> HH^0(X, Omega1D(X, -toricDivisor X))

tensor2 = X -> prune sheafHom(prune sheafHom(cotangentSheaf(X), sheaf(X, ring X)), cotangentSheaf(X))
randomInvariantP = X -> HH^1(X, tensor2(X))


l2 = {VV48v1, VV48v2, VV48v3, VV56v1, VV56v2, VV64v1, VV64v2, VV64v3}

for X in l2 do (<<randomInvariant2(X)<<endl;)

-- All VV48s are different
-- Third VV64 is different from other 2
-- That's what I know so far

-- According to randomInvariant2 = X -> HH^1(X, Omega1D(X, -toricDivisor X))
-- all 48s are different, all 

	    
for X in l do (<<vectorFields(X)<<endl;)

vectorFields(V5)
vectorFields(V7)

vectorFields(Fg327)
vectorFields(Fg328)
vectorFields(VV1)
vectorFields(VV2)

(0, {2,2}), (1, {3,2})

XX = projectiveBundleConstructor(0, {2,2})
YY = projectiveBundleConstructor(1, {3,2})

for i from 0 to 3 do (
    for coef from -1 to 1 do (
O	a = HH^i(XX, Omega1D(XX, coef*toricDivisor XX));
	b = HH^i(YY, Omega1D(YY, coef*toricDivisor YY));
	if not a == b then << i << ", " << coef << " differentiates them!" << endl;
    )
)

degreeOfThreefold(Fg327)
degreeOfThreefold(Fg328)

-- 2,2 takes a *long* time!!


-- V2 vs V9
-- V5 vs V6 vs V7

vectorFields(V2)
vectorFields(V9)

vectorFields(V5)
vectorFields(V6)
vectorFields(V7)

------


-- Klein-schmidt... on the way to producing 4-folds...

-- Can we classify Fanos and weak Fanos (toric) of Picard rank 2?
-- Picard rank 1 is ez, it is always projective space
-- Picard rank 2, for curves there is none, for surfaces it is Hirzebruch surfaces? Fano for a <= 1, weak Fano for a <= 2
-- For three folds, you need the Klein-Schmidt algorithm. And we need it bc we are doing 3-folds. 

for a1 from 0 to 20 do (
    for a2 from 0 to 20 do (
    	X = kleinschmidt(3, {a1, a2});
	if (a1 <= a2) then (
	    if (isFano X) then 
	    << a1 << ", " << a2 << " is Fano"  << endl
	    else if (not isFano X and isNef(-toricDivisor X)) then 
	    << a1 << ", " << a2 << " is weak Fano but not Fano" << endl;
	)
    )
)

-- Only 6! This is entirely doable!
-- Output
-- 0, 0 is Fano
-- 0, 1 is Fano
-- 0, 2 is weak Fano but not Fano
-- 1, 1 is weak Fano but not Fano
X = kleinschmidt(3, {0,0})
Y = kleinschmidt(3, {0,1})

for a1 from 0 to 20 do (
    X = kleinschmidt(3, {a1});
    if (isFano X) then 
    << a1 << " is Fano"  << endl
    else if (not isFano X and isNef(-toricDivisor X)) then 
    << a1 << " is weak Fano but not Fano" << endl;
)

-- 0 is Fano
-- 1 is Fano
-- 2 is Fano
-- 3 is weak Fano but not Fano
Z = kleinschmidt(3, {0})
genusOfThreefold(X)
genusOfThreefold(Y)
genusOfThreefold(Z)

vectorFields(X)
vectorFields(Y)
vectorFields(Z)
vectorFields(W)

randomInvariantP(X)
randomInvariantP(Y)
randomInvariantP(Z)

W = cartesianProduct(toricProjectiveSpace 1, toricProjectiveSpace 2)
rays W
rays X
rays Y
rays Z

-- Interestingly... the rays of W match those of Z exactly!
-- Meanwhile X seems to be a simple symmetry away from W and Z
-- Y seems genuinely different.
-- same story again! sameToricVariety does the job! As it is with batyrev, it is with kleinschmidt

So those are our three Fanos. What are our weak Fanos??
-- 0, 2 is weak Fano but not Fano
-- 1, 1 is weak Fano but not Fano
-- 3 is weak Fano but not Fano


wX = kleinschmidt(3, {0,2})
wY = kleinschmidt(3, {1,1})
wZ = kleinschmidt(3, {3})

genusOfThreefold(wX) -- 28
genusOfThreefold(wY) -- 28
genusOfThreefold(wZ) -- 37

vectorFields(wX) -- 38
vectorFields(wY) -- 35
vectorFields(wZ) -- 53

randomInvariantP(wX) -- 3
randomInvariantP(wY) -- 3
randomInvariantP(wZ) -- 10

-- So they are all different
-- So our complete list of Fanos+weak Fano 3-folds of Picard rank 2 is
-- kleinschmidt(3, {0,1}) -- Fano deg 28 (2-33)
-- kleinschmidt(3, {0}) -- Fano deg 28 (2-34), P^1 x P^2
-- kleinschmidt(3, {1}) -- Fano deg 29 (2-35)
-- kleinschmidt(3, {2}) -- Fano deg 32 (2-36)
-- kleinschmidt(3, {3}) -- weak Fano
-- kleinschmidt(3, {0,2}) -- weak Fano
-- kleinschmidt(3, {1,1}) -- weak Fano

-- These are the 7 weak Fano 3-folds of Picard rank 3. We may be interested in doing projectivations of O+O(D) on these, the same way as before, to produce fourfolds



------
-- Fourfolds of Picard rank 2 (for later)
for a1 from 0 to 10 do (
    for a2 from a1 to 10 do (
	for a3 from a2 to 10 do (
    	    X = kleinschmidt(4, {a1, a2, a3});
	    if (isFano X) then 
	    << a1 << ", " << a2 << ", " << a3 << " is Fano" << endl
	    else if (not isFano X and isNef(-toricDivisor X)) then 
	    << a1 << ", " << a2 << ", " << a3 << " is weak Fano but not Fano" << endl;
	)
    )
)





---
for a from 0 to 123 do (
    X = smoothFanoToricVariety(4, a);
    if rank picardGroup X == 3 then << a << " is Picard rank 3!" << endl;
)


--
