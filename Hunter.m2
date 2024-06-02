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
    A := random(kk^m,kk^m);
    B := random(kk^n,kk^n);
    N := A * X * B;
    product for j from 0 to #conjlam-1 list det(N_{0..conjlam_j-1}^{0..conjlam_j-1}));
    J := ideal lis;
    minJ:=mingens J;
    if rank source minJ != d then(
	error "Did not compute full ideal";
	);
    ideal mingens J
    )
