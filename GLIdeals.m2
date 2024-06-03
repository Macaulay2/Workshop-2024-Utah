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

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **DOCUMENTATION** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
beginDocumentation()


doc ///
    Key
        numgensILambda
	(numgensILambda, Matrix, List)
    Headline
        computes the number of generators of the I_Lambda ideal
    Usage
	numgensILambda(X, lam)
    Inputs
	X: Matrix
	    a matrix with dimensions r x c
	lam: List
	    a list of integers representing a partition
    Outputs
	: ZZ
	    the number of generators of the I_Lambda ideal
    Description
        Text
	    The `numgensILambda` function computes the number of generators of the I_Lambda ideal
	    associated with the input matrix `X` and partition `lam`. It uses the hook length formula
	    to compute the dimension of the space. The partition `lam` should be compatible with the 
	    dimensions of the matrix `X`.
	Example

    Caveat 

///


doc ///
    Key
        minimizeChi
	(minimizeChi, List)
    Headline
        computes the minimal elements of a list of partitions
    Usage
	minimizeChi(chi)
    Inputs
	chi: List
	    a list of partitions, where each partition is represented as a list of integers
    Outputs
	: List
	    a list of partitions, which are the minimal elements of the input list
    Description
        Text
	    The `minimizeChi` function takes a list of partitions and returns the minimal elements
	    of that list. A partition is minimal if it is not greater than any other partition
	    in the list according to the dominance order.
	Example
	    chi = {{3,2,1}, {3,1,1}, {2,2,2}, {1,1,1}}
	    minimizeChi(chi)
    Caveat 
///


doc ///
    Key
        numgensIChi
	(numgensIChi, Matrix, List)
    Headline
        computes the total number of generators for I_Chi ideals
    Usage
	numgensIChi(X, chi)
    Inputs
	X: Matrix
	    a matrix with dimensions r x c
	chi: List
	    a list of partitions, where each partition is represented as a list of integers
    Outputs
	: ZZ
	    the total number of generators for the I_Chi ideals
    Description
        Text
	    The `numgensIChi` function takes a matrix `X` and a list of partitions `chi` and computes the
	    total number of generators for the I_Chi ideals. It first minimizes the list of partitions using
	    `minimizeChi`, and then computes the sum of the number of generators for each minimal partition
	    using the `numgensILambda` function.
	Example

    Caveat 
///


doc ///
    Key
        minimizeChi
	(minimizeChi, List)
    Headline
        computes the minimal elements of a list of partitions
    Usage
	minimizeChi(chi)
    Inputs
	chi: List
	    a list of partitions, where each partition is represented as a list of integers
    Outputs
	: List
	    a list of partitions, which are the minimal elements of the input list
    Description
        Text
	    The `minimizeChi` function takes a list of partitions and returns the minimal elements
	    of that list. A partition is minimal if it is not greater than any other partition
	    in the list according to the dominance order. The function uses the hook length formula
	    to determine the dominance order among partitions.
	Example
	    chi = {{3,2,1}, {3,1,1}, {2,2,2}, {1,1,1}}
	    minimizeChi(chi)
    Caveat 
  		  
///


doc ///
    Key
        numgensIChi
        (numgensIChi, Matrix, List)
    Headline
        computes the total number of generators for I_Chi ideals
    Usage
	numgensIChi(X, chi)
    Inputs
	 X: Matrix
	    a matrix with dimensions r x c
	 chi: List
	    a list of partitions, where each partition is represented as a list of integers
    Outputs
	 : ZZ
	    the total number of generators for the I_Chi ideals
    Description
        Text
	    The `numgensIChi` function takes a matrix `X` and a list of partitions `chi` and computes the
	    total number of generators for the I_Chi ideals. It first minimizes the list of partitions using
	    `minimizeChi`, and then computes the sum of the number of generators for each minimal partition
	    using the `numgensILambda` function.
	Example 

    Caveat 
///


doc ///
    Key
        detLam
        (detLam, Matrix, List)
    Headline
        computes the product of determinants of specific submatrices
    Usage
	detLam(X, lam)
    Inputs
	 X: Matrix
	    a matrix with dimensions r x c
	 lam: List
	    a list of integers representing the sizes of submatrices
    Outputs
	 : RingElement
	    the product of the determinants of the specified submatrices
    Description
        Text
	    The `detLam` function takes a matrix `X` and a list of integers `lam`. For each integer in `lam`,
	    it computes the determinant of the submatrix of `X` formed by the first i rows and columns, where i
	    is an element of `lam`. It then returns the product of these determinants.
	Example 

    Caveat 
///


doc ///
    Key
        randomLam
        (randomLam, ZZ, ZZ)
    Headline
        generates a random partition of an integer
    Usage
	randomLam(n, k)
    Inputs
	 n: ZZ
	    the number of parts in the partition
	 k: ZZ
	    the integer to be partitioned
    Outputs
	 : List
	    a list representing a partition of the integer k into n parts
    Description
        Text
	    The `randomLam` function generates a random partition of the integer `k` into `n` parts. 
	    It constructs a list of `n` non-negative integers that sum to `k`, then sorts the list 
	    in non-increasing order and removes any zero entries.
	Example
	    randomLam(4, 10)
    Caveat 

///
	
doc ///
    Key
        idealILambda
        (idealILambda, Matrix, List)
    Headline
        constructs the I_Lambda ideal for a given matrix and partition
    Usage
	idealILambda(X, lam)
    Inputs
	 X: Matrix
	    a matrix with dimensions m x n
	 lam: List
	    a list of integers representing a partition
    Outputs
	 : Ideal
	    the I_Lambda ideal generated from the matrix and partition
    Description
        Text
	    The `idealILambda` function constructs the I_Lambda ideal for a given matrix `X` and a partition `lam`.
	    It first checks that the base ring of `X` has characteristic 0. Then, it generates random matrices `A` and `B`
	    to compute the matrix `N = A * X * B`. It calculates the determinants of submatrices of `N` based on the
	    conjugate of the partition `lam` and constructs the ideal from these determinants. The process ensures that
	    the rank of the source of the minimal generators of the ideal matches the required number of generators.
	Example

    Caveat 

///

doc ///
    Key
        idealIChi
        (idealIChi, Matrix, List)
    Headline
        constructs the sum of I_Lambda ideals for a list of partitions
    Usage
	idealIChi(X, chi)
    Inputs
	 X: Matrix
	    a matrix with dimensions m x n
	 chi: List
	    a list of partitions, where each partition is represented as a list of integers
    Outputs
	 : Ideal
	    the sum of the I_Lambda ideals for the given partitions
    Description
        Text
	    The `idealIChi` function constructs the sum of I_Lambda ideals for a given matrix `X` and a list of partitions `chi`.
	    It first minimizes the list of partitions using `minimizeChi`. If the list of partitions is empty, it returns the zero ideal.
	    For each minimal partition, it computes the corresponding I_Lambda ideal using `idealILambda` and sums these ideals to form the final result.
	Example

    Caveat 
///


doc ///
    Key
        partitionsLeq
        (partitionsLeq, Partition, Partition)
    Headline
        checks if one partition is less than or equal to another
    Usage
	partitionsLeq(A, B)
    Inputs
	 A: Partition
	    a partition represented as a list of integers in weakly decreasing order
	 B: Partition
	    another partition represented as a list of integers in weakly decreasing order
    Outputs
	 : Boolean
	    true if partition A is less than or equal to partition B, false otherwise
    Description
        Text
	    The `partitionsLeq` function checks if one partition `A` is less than or equal to another partition `B`
	    in the dominance order. It assumes that both partitions are represented as lists of integers in weakly
	    decreasing order. The function iterates through the parts of `A` and `B` and returns false if any part
	    of `A` is greater than the corresponding part of `B`.
	Example
	    A = {3, 2, 1}
	    B = {4, 2, 1}
	    partitionsLeq(A, B)
	    partitionsLeq(B, A)
    Caveat 
  
///

doc ///
    Key
        ?
        (?, Partition, Partition)
    Headline
        compares two partitions and returns their relationship
    Usage
	A ? B
    Inputs
	 A: Partition
	    a partition represented as a list of integers in weakly decreasing order
	 B: Partition
	    another partition represented as a list of integers in weakly decreasing order
    Outputs
	 : Symbol
	    the symbol representing the relationship between partitions A and B:
	    `<` if A is less than B, `>` if A is greater than B, `==` if A equals B, and `incomparable` if they are incomparable
    Description
        Text
	    The `Partition ? Partition` function compares two partitions `A` and `B` and returns a symbol
	    representing their relationship. It first checks if `A` is less than or equal to `B` and if `B` is
	    less than or equal to `A` using the `partitionsLeq` function. Depending on the results of these checks,
	    it returns `<` if `A` is less than `B`, `>` if `A` is greater than `B`, `==` if `A` equals `B`, and 
	    `incomparable` if neither partition is less than or equal to the other.
        Example
	    A = {3, 2, 1}
	    B = {4, 2, 1}
	    A ? B
	    B ? A
	    C = {3, 3, 1}
	    D = {3, 2, 2}
	    C ? D
    Caveat 
///

doc ///
    Key
        ==
        (==, Partition, Partition)
    Headline
        checks if two partitions are equal
    Usage
	A == B
    Inputs
	 A: Partition
	    a partition represented as a list of integers in weakly decreasing order
	 B: Partition
	    another partition represented as a list of integers in weakly decreasing order
    Outputs
	 : Boolean
	    true if partition A is equal to partition B, false otherwise
    Description
        Text
	    The `Partition == Partition` function checks if two partitions `A` and `B` are equal. 
	    It uses the comparison function `A ? B` and returns true if the result is the symbol `==`,
	    indicating that the partitions are equal.
        Example
            A = {3, 2, 1}
	    B = {3, 2, 1}
	    C = {4, 2, 1}
	    A == B
	    A == C
    Caveat 
 
///



------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **TESTS** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------

TEST ///

///



end



