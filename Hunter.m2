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

Y=genericMatrix(S,n,m)

phi1(Y)
phi1'=inverse phi1

degree xx_(2,2)


numgensILambda(X,{2,1,0})
I1=idealILambda(X,{1,1,1});
I2=idealILambda(X,{2,1,0});
J
J=phi1'(I1+I2);

mingens J
unique degrees J


