newPackage(
	"GLIdeals",
    	Version => "0.1", 
    	Date => "June 1, 2024",
    	Authors => {
	     {Name => "Hunter Simper", Email => "hunter.simper@utah.edu", HomePage => "https://www.huntersimper.com"},
	     {Name => "Van Vo"}, 
		 {Name => "Aryaman Maithani", Email => "maithani@math.utah.edu", HomePage => "https://aryamanmaithani.github.io/math/"}
	     },
    	Headline => "Template for gl ideals package",
	PackageExports => {"SchurRings","SimpleDoc"},
    	DebuggingMode => false,
	Reload=>true
    	)
    
export{"idealILambda","numgensILambda", "idealToChi"}

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

    minimals := new List from {};
	c := apply(chi, L -> new Partition from L);
    i := 0;
	while i < #c do(
        isMinimalPartition := true;
        lambda := c#i;
        j := i + 1;
        while j < #c do(
            compare := lambda ? c#j;
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
detLam (Matrix, List) := (X,lam) -> (
	conjlam := toList conjugate( new Partition from lam);
    C := 1;
    for i from 0 to #conjlam-1 do(
	C = C*determinant(submatrix(X,{0..conjlam_i - 1},{0..conjlam_i - 1}));
    );
    return C 
)

randomLam = method();
randomLam(ZZ, ZZ) := (n,k) -> (
    L:=new MutableList from {};
    sumsofar := 0;
    for i from 0 to (n-2) do(
        x := random(0, k-sumsofar);
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
    A := random(kk^n,kk^n);
    B := random(kk^m,kk^m);
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

idealIChi = method()
idealIChi(Matrix,List) := (X,chi) -> (
	-- chi is a list of partitions
	-- return sum of ideals idealILambda(X, chi_i) for all i
	R := ring X;
	if #chi == 0 then return ideal(0_R);
	c := minimizeChi(chi);
	I := idealILambda(X, c#0);
	for i in 1..#c-1 do(
		I = I + idealILambda(X, c#i);
	);
	return I;
);

goodDegree = method();
goodDegree(ZZ, ZZ, List) := (n, m, L) -> (
    -- L is a length n + m list of integers
    -- return if L is a good degree for n, m
    -- i.e., if L = (lam, lam) for a partition lam (suitably padded)
    edegree := take(L, n);      -- the first n elements of L
    fdegree := drop(L, n);      -- the remaining elements
    
    isPartition := P -> (rsort(P) == P);
    
    if not (isPartition(edegree) and isPartition(fdegree)) then return false;
    -- edegree and fdegree are weakly decreasing

    edegree = delete(0, edegree);    -- remove the trailing 0s
    fdegree = delete(0, fdegree);    -- remove the trailing 0s
    return edegree == fdegree;
);

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

idealToChi = method();
idealToChi(Matrix, Ideal) := (Y, J) -> (
	-- J is an ideal
	-- return the list of partitions chi such that J = idealIChi(Y, chi)

	n := numRows Y;
	m := numColumns Y;
	(R, X, phi) := correctSymmetricAlgebraHelper(Y);
	I := phi^(-1)(J);

	-- now we have the ideal I in R
	L := flatten entries mingens I;
	-- print L;
	degrees := apply(L, degree);
	goodDegrees := select(degrees, deg -> goodDegree(n, m, deg));

	possiblePartitions := apply(goodDegrees, deg -> take(deg, n));
	possiblePartitions = unique(possiblePartitions);
	possiblePartitions = apply(possiblePartitions, P -> delete(0, P));

	return select(possiblePartitions, P -> (detLam(X, P) % I == 0));
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



