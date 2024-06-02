newPackage(
	"GLIdeals",
    	Version => "0.1", 
    	Date => "June 1, 2024",
    	Authors => {
	     {Name => "Hunter Simper", Email => "hunter.simper@utah.edu", HomePage => "https://www.huntersimper.com"},
	     {}
	     },
    	Headline => "Template for gl ideals package",
	PackageExports => {"SchurRings","SimpleDoc"},
    	DebuggingMode => false,
	Reload=>true
    	)
    
export{"generateIlambda"}


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

mingensILambda = method()
mingensILambda(Matrix, List) := (X, lam) - (
	myHookLength := P -> (
		-- compute prod (P_i - P_j + j - i) / (j - i) for all i < j
		num := 1;
		den := 1;
		for i in 0..#P-1 do(
			for j in i+1..#P-1 do(
				num = num * (P#i - P#j + j - i);
				den = den * (j - i);
			);
		);
		return num / den;
	);

	lam := local lam;
	size := #lam;

	r := numRows X;
	c := numColumns X;
	min_size := min(r, c);
	max_size := max(r, c);

	if min_size < size then(
		error "Partition is too large for the matrix";
	);

	lam = lam | apply(min_size - size, i -> 0); -- make it size min_size
	dimension := myHookLength(lam);
	lam = lam | apply(max_size - min_size, i -> 0); -- make it size max_size
	return dimension * myHookLength(lam);
) 
-----

beginDocumentation()


doc ///

///


end




