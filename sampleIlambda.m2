restart
needsPackage"SchurRings"
m = n = 4
kk = ZZ/101
S = kk[x_(1,1)..x_(m,n)]

generateIlambda = method()
generateIlambda(ZZ,ZZ,List,PolynomialRing) := (n,m,lam,S) -> (
     r := local r;
     s := local s;
     x := local x;
     R := schurRing(r,n);
     T := schurRing(s,m);
     conjlam := toList conjugate( new Partition from lam);
     d := dim r_lam;
     e := dim s_lam;
     M := genericMatrix(S,m,n);
     lis := for i from 0 to d*e-1 list
     (
    A := random(kk^m,kk^m);
    B := random(kk^n,kk^n);
    N := A * M * B;
    product for j from 0 to #conjlam-1 list det(N_{0..conjlam_j-1}^{0..conjlam_j-1}));
    J := ideal lis;
    ideal mingens J
    )

I1 = generateIlambda(4,4,{1},S);

gens I1
I2 = generateIlambda(4,4,{2,2,2,1},S);
