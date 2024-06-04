restart





correctSymmetricAlgebraHelper = method();
correctSymmetricAlgebraHelper(Matrix) := Y -> (
	-- take in a matrix Y
	-- S is the ring that Y is defined over
	-- assuming that S = k[Y]
	-- we will construct an isomorphic ring R = k[X]
	-- with a specific grading where deg(x_(i, j)) = (e_i, f_j) \in Z^n x Z^m = Z^{n+m}
	-- define phi : S -> R, where the entries of Y are mapped to the entries of X (correspondingly...)
	-- return (R, X, phi)
	
	x := symbol x;
	n := numRows Y;
	m := numColumns Y;
	S := ring Y;
	k := baseRing S;

	indices := toList((1, 1)..(n, m));
    myVars := apply(indices, ind -> x_(ind));
    edegree := i -> toList(insert(i-1, 1, ((n-1) : 0)));    -- e(i) = (0, ..., 0, 1, 0, ..., 0) \in Z^n
    fdegree := j -> toList(insert(j-1, 1, ((m-1) : 0)));    -- f(j) = (0, ..., 0, 1, 0, ..., 0) \in Z^m
    myDegrees := apply(indices, ind -> edegree(ind#0) | fdegree(ind#1));
    R := k[myVars, Degrees => myDegrees];
    X := transpose genericMatrix(R, m, n);

	phi := map(S, R, flatten entries Y);

	return (R, X, phi);
)


naiveClosure = method();---
naiveClosure (Matrix, Ideal) := (Y,I) ->( 
    (R,X,phi):=correctSymmetricAlgebraHelper(Y);
    kk=baseRing ring Y;
    n:=rank target X;
    m:=rank source X;
    II:=(inverse(phi))(I);
    A:=random(kk^n,kk^n, MaximalRank=>true);
    B:=random(kk^m,kk^m, MaximalRank=>true);
    act:=map(R,R,flatten entries(A*X*B));
    JJ:=II+act(II);   
    while not(II==JJ) do(
    II=JJ;
    A=random(kk^n,kk^n, MaximalRank=>true);
    B=random(kk^m,kk^m, MaximalRank=>true);
    act=map(R,R,flatten entries(A*X*B));
    JJ=II+act(II);    
    );
return (ideal mingens phi(II));
    )

idealILambda2 = method()
idealILambda2(Matrix,List) := (X,lam) -> (
     n:=rank target X;
     m:=rank source X;
     kk:=baseRing ring X; 
     if char kk !=0 then(
		error "Base ring is not characteristic 0";
	);
     d:=numgensILambda(X,lam);
     (R,Y,phi):=correctSymmetricAlgebraHelper(X);
     gen:=detLam(Y,lam);
     lis := for i from 0 to d-1 list
     (
    A:=random(kk^n,kk^n);
    B:=random(kk^m,kk^m);
    act:=map(R,R,flatten entries(A*Y*B));
    act(gen));
    J := ideal lis;
    minJ:=mingens J;
    	while rank source minJ != d do(
		A = random(kk^n,kk^n);
		B = random(kk^m,kk^m);
		act=map(R,R,flatten entries(A*Y*B));
		lis = append(lis,act(gen));
		J = ideal lis;
		minJ = mingens J;		
	);
	return phi(ideal minJ);
);

needsPackage"GLIdeals"

S=QQ[x_1..x_9];
Y=genericMatrix(S,3,3);


numgensILambda(Y,{2})

for i from 1 to 10 do(
    p=randomLam(3,2);
    print (p,numgensILambda(Y,p));
    elapsedTime(idealILambda(Y,p));
    elapsedTime(idealILambda2(Y,p));
    )

