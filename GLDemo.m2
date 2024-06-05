restart
loadPackage"GLIdeals"


n=3
m=3
R=QQ[x_1..x_(n*m)]
X=genericMatrix(R,n,m)



idealToChi(X,ideal(1_R))
idealToChi(X,ideal(0_R))

lam1={2,1,1}
numgensGLIdeal(X,lam1)
I1=GLIdeal(X,lam1,MaximalRank=>false);
--transpose mingens I1


lam2={2,2}
numgensGLIdeal(X,lam2)
I2=GLIdeal(X,lam2,MaximalRank=>false);
--transpose mingens I2

J=I1+I2;

I=minors(2,X);
idealToChi(X,I^2)
J==I1+I2

T=Tor_1(R^1/I1,R^1/I2);
AT=ann T;
idealToChi(X, AT)

needsPackage"SymbolicPowers"
K=symbolicPower(I,2);
idealToChi(X,K)
K==AT





