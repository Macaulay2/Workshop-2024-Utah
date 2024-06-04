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


-- Complete smooth toric weak Fano ( \ Fano) 3-folds of Picard rank 3 (w/ 5 primitive collections)

V1 = batyrevConstructor({2,1,1,1,1}, {2}, {});
V2 = batyrevConstructor({1,1,2,1,1}, {0}, {1});
V3 = batyrevConstructor({1,2,1,1,1}, {2}, {});
V4 = batyrevConstructor({1,1,1,2,1}, {0,0}, {});


isNef(-toricDivisor V)

degreeOfThreefold = X -> (
    K := toricDivisor X;
    c := chern(1, OO (-K));
    integral (c*c*c)
)


degreeOfThreefold(V1)
degreeOfThreefold(V2)
degreeOfThreefold(V3)
degreeOfThreefold(V4)


-- This one is actually Fano
V5 = batyrevConstructor({1,1,2,1,1}, {0}, {0});
degreeOfThreefold(V4)


