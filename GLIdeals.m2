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
    
export{"idealILambda"}

numgensILambda = method() -- TODO: Make it accept Partitions as well
numgensILambda(Matrix, List) := (X, lam) -> (
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

	size := #lam;

	r := numRows X;
	c := numColumns X;
	minSize := min(r, c);
	maxSize := max(r, c);

	if minSize < size then(
		error "Partition is too large for the matrix";
	);

	P := lam | apply(minSize - size, i -> 0); -- make it size minSize
	dimension := myHookLength(P);
	P = P | apply(maxSize - minSize, i -> 0); -- make it size maxSize
	return sub(dimension * myHookLength(P),ZZ);
);

minimizeChi = method();
minimizeChi(List) := (chi) -> (
	-- chi is a list of partitions (possibly of Type List)
	-- return a list of partitions (of Type List) 
	-- these are the minimal elements of chi

    minimals = new List from {};
	c := apply(chi, L -> new Partition from L);
    i := 0;
	while i < #c do(
        isMinimalPartition := true;
        lambda := c#i;
        j := i + 1;
        while j < #c do(
            compare = lambda ? c#j;
            j = j + 1;
            if compare == incomparable then continue;
            if compare == (symbol >) then (isMinimalPartition = false; break;);
            
            j = j - 1;
            -- we have lambda <= c#j at this point
            c = drop(c, {j, j});    -- remove c#j
        );
        if isMinimalPartition then minimals = append(minimals, lambda);
        i = i + 1;
    );
    return apply(minimals, P -> toList(P));
)

numgensIChi = method();
numgensIChi(Matrix, List) := (X, chi) -> (
 	-- chi is a list of partitions
 	-- return sum numgensILambda(chi_i) for all i
	c := minimizeChi(chi);
	s := 0;
	for lam in c do(
		s = s + numgensILambda(X, lam);
	);
	return s;
);

detLam = method()
detLam = (X,lam) -> (
    (m,n) := (numColumns(X), numRows(X));
    base := ring X;
    C := 1;
    for i in lam do(
	C = C*determinant(submatrix(X,{0..i},{0..i}));
    );
    return C 
    )

randomLam = method();
randomLam(ZZ, ZZ) := (n,k) -> (
    L=new MutableList from {};
    sumsofar = 0;
    for i from 0 to (n-2) do(
        x = random(0, k-sumsofar);
        L#i=x;
        sumsofar = sumsofar + x;
        );
    L#(n-1) = k - sumsofar;
    L=toList L;
    L=rsort L;
    return delete(0,L)
    )

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
    
	while rank source minJ != d do(
		A := random(kk^m,kk^m);
		B := random(kk^n,kk^n);
		N := A * X * B;
		lis = append(lis, product for j from 0 to #conjlam-1 list det(N_{0..conjlam_j-1}^{0..conjlam_j-1}));
		J = ideal lis;
		minJ = mingens J;
	);
	return ideal minJ;
);



partitionsLeq = method();
partitionsLeq(Partition, Partition) := (A, B) -> (
    -- A and B are partitions
    -- assuming weakly decreasing
    -- return if A <= B
    n := #B;
    for i in 0..(#A-1) do(
        if A#i == 0 then break;
        if i >= n then return false;
        if A#i > B#i then return false;
    );
    return true;
)

Partition ? Partition := (A, B) -> (
    AleB := partitionsLeq(A, B);
    BleA := partitionsLeq(B, A);
    if AleB and BleA then return symbol ==;
    if AleB then return symbol <;
    if BleA then return symbol >;
    return symbol incomparable;
)

Partition == Partition := (A, B) -> (
    return (A ? B) == (symbol ==);
)
-----

beginDocumentation()


doc ///

///


end



