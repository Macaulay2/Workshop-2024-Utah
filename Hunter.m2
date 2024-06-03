idealILambda = method()
idealILambda(Matrix,List) := (X,lam) -> (
     n:=rank target X;
     m:=rank source X;
     kk:=baseRing ring X; 
     if char kk !=0 then(
		error "Base ring is not characteristic 0";
	);
     conjlam := toList conjugate( new Partition from lam);
     d:=numgensILambda(X,lam);
     lis := for i from 0 to d-1 list
     (
    A := random(kk^n,kk^n);
    B := random(kk^m,kk^m);
    N := A * X * B;
    product for j from 0 to #conjlam-1 list det(N_{0..conjlam_j-1}^{0..conjlam_j-1}));
    J := ideal lis;
    minJ:=mingens J;
    if rank source minJ != d then(
	error "Did not compute full ideal";
	);
    ideal mingens J
    )


restart
loadPackage"GLIdeals"


R=QQ[x_1..x_12];
X=genericMatrix(R,3,4)
n=rank target X
m=rank source X

degs=flatten apply(m,i->entries id_(ZZ^n))

--degs=join(entries id_(ZZ^n),apply(n*(m-1),i->apply(n,i->0)))
S=QQ[xx_(1,1)..xx_(m,n), Degrees=>degs, MonomialOrder=>Lex]
vars S

phi1=map(R,S,flatten entries transpose X) 



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


S=QQ[x_1..x_6];
Y=genericMatrix(S,2,3);

baseRing ring Y

(a,b,c)=correctSymmetricAlgebraHelper(Y)
(inverse c)(x_1) 
