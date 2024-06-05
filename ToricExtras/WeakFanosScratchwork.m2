--- just some stuff I was playing around with


	
for n from 2 to 5 list (
    for b from 0 to 5 list (
	V = batyrevConstructor({1,n,1,1,1}, {b}, {});
	<< isFano(V);
	<< isNef(-toricDivisor V); -- but it is weak Fano
	<< n << " " << b << endl;))

for c from 0 to 5 list (
    for b from 0 to 5 list (
	V = batyrevConstructor({1,1,2,1,1}, {b}, {c});
	<< isFano(V);
	<< isNef(-toricDivisor V); -- but it is weak Fano
	<< "b is " << b << ", c is " << c << endl;))
	
for n from 2 to 5 list (
    for b from 0 to 5 list (
	V = batyrevConstructor({1,1,1,1,n}, {b}, {});
	<< isFano(V);
	<< isNef(-toricDivisor V); -- but it is weak Fano
	<< n << " " << b << endl;))

for b0 from 0 to 5 list
	;
	
	
	
for b0 from 0 to 5 list (
    for b1 from 0 to 5 list (
	V = batyrevConstructor({1,1,1,2,1}, {b0, b1}, {});	
	<< isFano(V);
	<< isNef(-toricDivisor V); -- but it is weak Fano
	<< b0 << " " << b1 << endl;))

-- Complete list of smooth toric Fano 3-folds of Picard rank 3 (w/ 5 primitive collections)

V26 = batyrevConstructor({2,1,1,1,1}, {0}, {});
V29 = batyrevConstructor({1,1,2,1,1}, {1}, {});
V26p = batyrevConstructor({1,1,2,1,1}, {0}, {}); -- isomorphic to V26 above


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


-- Question: is this the complete list?

V = V5
isFano(V)
isNef(-toricDivisor V)

degreeOfThreefold = X -> (
    K := toricDivisor X;
    c := chern(1, OO (-K));
    integral (c*c*c)
)

degreeOfThreefold(V)

degreeOfThreefold(V1)
degreeOfThreefold(V2)
degreeOfThreefold(V3)
degreeOfThreefold(V4)

-- Moreover, are V2 and V4 isomorphic?

for b from 0 to 2 list (
    V = batyrevConstructor({2,1,1,1,1}, {b}, {});
    << isFano(V);
    << isNef(-toricDivisor V); -- but it is weak Fano
    << "b=  " << b << endl;)

for b from 0 to 2 list (
    V = batyrevConstructor({1,2,1,1,1}, {b}, {});
    << isFano(V);
    << isNef(-toricDivisor V); -- but it is weak Fano
    << "b=  " << b << endl;)

V = batyrevConstructor({1,1,2,1,1}, {0}, {0});
<< isFano(V);
<< isNef(-toricDivisor V); -- but it is weak Fano
<< "b=  " << b << endl;

V = batyrevConstructor({1,1,2,1,1}, {1}, {0});
<< isFano(V);
<< isNef(-toricDivisor V); -- but it is weak Fano
<< "b=  " << b << endl;

V = batyrevConstructor({1,1,2,1,1}, {0}, {1});
<< isFano(V);
<< isNef(-toricDivisor V); -- but it is weak Fano
<< "b=  " << b << endl;

for b from 0 to 1 list (
    V = batyrevConstructor({1,1,1,1,2}, {b}, {});
    << isFano(V);
    << isNef(-toricDivisor V); -- but it is weak Fano
    << "b=  " << b << endl;)


PP1 = toricProjectiveSpace 1
FF2 = hirzebruchSurface 2
prod = cartesianProduct (PP1, FF2)

degreeOfThreefold(prod)

-- write some code to compute Hodge number h^{1,2} for a 3-fold
-- write some code to compute primitive collections given a normalToricVariety (just loop over all subsets is a viable brute force approach)
