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
	PackageExports => {"SimpleDoc"},
    	DebuggingMode => false,
	Reload=>true
    	)
    
export{"idealILambda","numgensILambda", "idealToChi", "naiveClosure","detLam","randomLam","IsMinimal", "GLIdeal", "numgensGLIdeal"}

numgensILambda = method()
numgensILambda(ZZ, ZZ, List) := (n, m, lam) -> (
	-- lam is a parition
	-- n, m are integers specifying the dimensions for GL_n x GL_m
	-- return the dimension of I_lam in the |lam|-th degree

	lam = rsort lam; -- make it weakly decreasing
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

	minSize := min(n, m);
	maxSize := max(n, m);

	if minSize < size then(
		error "Partition is too large for the matrix";
	);

	P := lam | apply(minSize - size, i -> 0); -- make it size minSize
	dimension := myHookLength(P);
	P = P | apply(maxSize - minSize, i -> 0); -- make it size maxSize
	return sub(dimension * myHookLength(P),ZZ);
);
numgensILambda(ZZ,ZZ,Partition) := (n,m,P)->(
    return numgensILambda(n,m,toList P);
);
    
numgensILambda(Matrix, List) := (X, lam) -> (
	return numgensILambda(numRows X, numColumns X, lam);
);
numgensILambda(Matrix,Partition):=(X,P)->(
    return numgensILambda(X,toList P);
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


numgensIChi = method(Options => {IsMinimal => false});
numgensIChi(ZZ,ZZ, List) := opts->(n,m, chi) -> (
 	-- chi is a list of partitions
 	-- return sum numgensILambda(chi_i) for all i
	c := apply(chi, rsort); -- make each element weakly decreasing
	if not (opts#IsMinimal) then c = minimizeChi(chi);
	s := 0;
	for lam in c do(
		s = s + numgensILambda(n,m, lam);
	);
	return s;
);
numgensIChi(Matrix, List) := opts->(X, chi) -> (
	return numgensIChi(numRows X, numColumns X, chi, opts);
);


detLam = method()
detLam (Matrix, List) := (X,lam) -> (
	conjlam := toList conjugate( new Partition from (rsort lam));
    C := 1;
    for i from 0 to #conjlam-1 do(
	C = C*determinant(submatrix(X,{0..conjlam_i - 1},{0..conjlam_i - 1}));
    );
    return C;
)
detLam (Matrix, Partition) := (X, P) -> (
	return detLam(X, toList P);
);

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
    return delete(0,L);
);

idealILambda = method(Options => {MaximalRank => true});
idealILambda(Matrix,List) := opts -> (X,lam) -> (
    kk:=baseRing ring X; 
    if char kk !=0 then(
		I := ideal(detLam(X,lam));
		return naiveClosure(X,I,MaximalRank => opts#MaximalRank);
	);

	lam = rsort lam; -- make it weakly decreasing
    n:=rank target X;
    m:=rank source X;
    conjlam := toList conjugate( new Partition from lam);
    d:=numgensILambda(X,lam);
    
	myRand := n -> random(kk^n,kk^n,MaximalRank=>opts#MaximalRank);

	lis := for i from 0 to d-1 list(
    	A := myRand(n);
    	B := myRand(m);
    	N := A * X * B;
    	product for j from 0 to #conjlam-1 list det(N_{0..conjlam_j-1}^{0..conjlam_j-1})
	);
    J := ideal lis;
    minJ:=mingens J;
    
	while rank source minJ != d do(
		A := myRand(n);
    	B := myRand(m);
		N := A * X * B;
		lis = append(lis, product for j from 0 to #conjlam-1 list det(N_{0..conjlam_j-1}^{0..conjlam_j-1}));
		J = ideal lis;
		minJ = mingens J;
	);
	return ideal minJ;
);
idealILambda(ZZ, ZZ, List) := opts -> (n, m, lam) -> (
	X := symbol X;
	R := QQ[X_(1,1)..X_(n,m)];
	return idealILambda(transpose genericMatrix(R, m, n), lam, opts);
);

idealIChi = method(Options => {IsMinimal => false, MaximalRank => true});
idealIChi(Matrix,List) := opts -> (X,chi) -> (
	-- chi is a list of partitions
	-- return sum of ideals idealILambda(X, chi_i) for all i
	if #chi == 0 then return ideal(0_(ring X));
	c := apply(chi, rsort); -- make each element weakly decreasing
	if not (opts#IsMinimal) then c = minimizeChi(chi);
	I := idealILambda(X, c#0, MaximalRank => opts#MaximalRank);
	for i in 1..#c-1 do(
		I = I + idealILambda(X, c#i, MaximalRank => opts#MaximalRank);
	);
	return I;
);
idealIChi(ZZ, ZZ, List) := opts -> (n, m, chi) -> (
	X := symbol X;
	R := QQ[X_(1,1)..X_(n,m)];
	return idealIChi(transpose genericMatrix(R, m, n), chi, opts);
);

naiveClosure = method(Options => {MaximalRank => true, Limit => false});
naiveClosure (Matrix, Ideal) := opts -> (Y,I) ->( 
    (R,X,phi):=correctSymmetricAlgebraHelper(Y);
    kk:=baseRing ring Y;
	myRand := n -> random(kk^n,kk^n,MaximalRank=>opts#MaximalRank);

	keepGoing := i -> true;
	if not (opts#Limit === false) then keepGoing = (i -> i < opts#Limit);

    n:=rank target X;
    m:=rank source X;
    II:=(inverse(phi))(I);
    A:=myRand(n);
    B:=myRand(m);
    act:=map(R,R,flatten entries(A*X*B));
    JJ:=II+act(II);
	counter := 0;
	success := (II == JJ);
    while keepGoing(counter) and not(success) do(
		II=JJ;
		A=myRand(n);
		B=myRand(m);
		act=map(R,R,flatten entries(A*X*B));
		JJ=II+act(II);
		counter = counter + 1;
		success = (II == JJ);
    );
	if opts#Limit === false then
		return (ideal mingens phi(II));
	return (ideal mingens phi(II), success);
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
	-- define phi : R -> S, where the entries of Y are mapped to the entries of X (correspondingly...)
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
	degrees := apply(L, degree);
	-- goodDegrees := select(degrees, deg -> goodDegree(n, m, deg));

	goodDegrees := {};
	for deg in degrees do(
		if ((not member(deg, goodDegrees)) and goodDegree(n, m, deg)) then goodDegrees = append(goodDegrees, deg);
	);

	possiblePartitions := apply(goodDegrees, deg -> take(deg, n));
	possiblePartitions = apply(possiblePartitions, P -> delete(0, P));

	return select(possiblePartitions, P -> (detLam(X, P) % I == 0));
);

IsPartition = method();
IsPartition(List) := (L) -> (
	-- L is a list
	-- return if L is a partition (either as a List or a Partition) or a list of partitions
	-- {} returns false
	-- (new Partition from {}) returns true
	-- DOESN'T CHECK IF THE PARTITION IS WEAKLY DECREASING
	if class(L) === Partition then return true;
	if #L == 0 then return false;
	return class(L#0) === ZZ;
);

GLIdeal = method(Options => {IsMinimal => false, MaximalRank => true}); -- call the appropriate idealILambda or idealIChi
GLIdeal(ZZ, ZZ, List) := opts -> (n, m, L) -> (
	if IsPartition(L) then return idealILambda(n, m, L, MaximalRank => opts#MaximalRank);
	return idealIChi(n, m, L, opts);
);
GLIdeal(ZZ, ZZ, Partition) := opts -> (n, m, P) -> (
	return idealILambda(n, m, P, MaximalRank => opts#MaximalRank);
);
GLIdeal(Matrix, List) := opts -> (X, L) -> (
	if IsPartition(L) then return idealILambda(X, L, MaximalRank => opts#MaximalRank);
	return idealIChi(X, L, opts);
);
GLIdeal(Matrix, Partition) := opts -> (X, P) -> (
	return idealILambda(X, P, MaximalRank => opts#MaximalRank);
);

numgensGLIdeal = method(Options => {IsMinimal => false});
numgensGLIdeal(ZZ, ZZ, List) := opts -> (n, m, L) -> (
	if IsPartition(L) then return numgensILambda(n, m, L);
	return numgensIChi(n, m, L, opts);
);
numgensGLIdeal(ZZ, ZZ, Partition) := opts -> (n, m, P) -> (
	return numgensILambda(n, m, P);
);
numgensGLIdeal(Matrix, List) := opts -> (X, L) -> (
	if IsPartition(L) then return numgensILambda(X, L);
	return numgensIChi(X, L, opts);
);
numgensGLIdeal(Matrix, Partition) := opts -> (X, P) -> (
	return numgensILambda(X, P);
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










------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **DOCUMENTATION** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
beginDocumentation()

doc ///
    Key
        numgensILambda
	(numgensILambda, ZZ, ZZ, List)
	(numgensILambda, Matrix, List)
    Headline
        computes the number of generators of the I_Lambda ideal.
    Usage
	numgensILambda(n,m, lam)
    Inputs
	n: ZZ
	    an integer 
	m: ZZ
	    an integer 
	lam: List
	    a list of integers representing a partition
    Outputs
	: ZZ
	    the number of generators of the I_Lambda ideal, 
    Description
        Text
	    This function computes the number of generators of the I_Lambda ideal in the ring Sym(CC^n,CC^m). Inputting a matrix X  will set n and m to be the number of rows and columns of X.
    	Example
	    S=QQ[x_1..x_5,y_1..y_5];
	    X=transpose genericMatrix(S,5,2)
	    L={2,1}
	    numgensILambda(X,L)
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
            This function takes a list of partitions and returns the minimal elements of that list.
        Example
            Chi = {{3,2,1}, {3,1,1}, {2,2,2}, {1,1,1}}
            minimizeChi(Chi)
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
            This function takes a matrix X and a list of partitions chi and computes the total number of generators for the I_Chi ideals.
        Example
            S=QQ[x_(1,1)..x_(3,5)];
            X=transpose genericMatrix(S,5,3)
            L=apply(partitions(3), toList)
            numgensIChi(X,L)
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
            This function takes a matrix X and a list of integers lam.
        Example
            S=QQ[x_(1,1)..x_(3,5)];
            X=transpose genericMatrix(S,5,3)
            L={3,2}
            detLam(X,L)
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
	    a list representing a partition of the integer k into at most n parts
    Description
        Text
            This function generates a random partition of the integer k into at most n parts.
        Example
            randomLam(4, 10)

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
            This function constructs the I_Lambda ideal for a given matrix X and a partition lam.

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
            This function constructs the sum of I_Lambda ideals for a given matrix X and a list of partitions chi.
///

doc ///
    Key
        naiveClosure
        (naiveClosure, Matrix, Ideal)
    Headline
        computes the naive closure of an ideal under a matrix action
    Usage
        naiveClosure(Y, I)
    Inputs
         Y: Matrix
            a matrix representing the mapping
         I: Ideal
            the ideal to be closed
    Outputs
         : Ideal
            the ideal representing the naive closure
    Description
        Text
            This function computes the naive closure of the given ideal I under the action of the matrix Y.
            The process involves repeatedly applying random maximal rank matrices to transform the ideal until it stabilizes.
        Example
            naiveClosure(Y, I)

///

doc ///
    Key
        goodDegree
        (goodDegree, ZZ, ZZ, List)
    Headline
        checks if a list of integers is a good degree for given dimensions
    Usage
        goodDegree(n, m, L)
    Inputs
         n: ZZ
            the length of the first part of the list
         m: ZZ
            the length of the second part of the list
         L: List
            a list of integers of length n + m
    Outputs
         : Boolean
            true if the list represents a good degree, false otherwise
    Description
        Text
            This function checks if a given list of integers L is a good degree for dimensions n and m.
            A list is considered a good degree if it can be split into two parts, each forming a weakly decreasing sequence, and these sequences are equal after removing trailing zeros.
        Example
            goodDegree(3, 3, {3, 2, 1, 3, 2, 1}) 
            goodDegree(3, 3, {3, 2, 1, 3, 2, 0})
    Caveat
        The input list L must have a length equal to n + m.

///

doc ///
    Key
        correctSymmetricAlgebraHelper
        (correctSymmetricAlgebraHelper, Matrix)
    Headline
        constructs an isomorphic ring with a specific grading
    Usage
        correctSymmetricAlgebraHelper(Y)
    Inputs
         Y: Matrix
            a matrix over a ring S = k[Y]
    Outputs
         : Sequence
            a sequence (R, X, phi) where R is the isomorphic ring, X is a matrix over R, and phi is the map from S to R
    Description
        Text
            This function constructs an isomorphic ring R with a specific grading for a given matrix Y over a ring S.
        Example
            correctSymmetricAlgebraHelper(Y)

///

doc ///
    Key
        idealToChi
        (idealToChi, Matrix, Ideal)
    Headline
        finds partitions corresponding to an ideal
    Usage
        idealToChi(Y, J)
    Inputs
         Y: Matrix
            a matrix over a ring
         J: Ideal
            the ideal to be analyzed
    Outputs
         : List
            a list of partitions chi such that J = idealIChi(Y, chi)
    Description
        Text
            This function finds the list of partitions chi such that the given ideal J is equal to idealIChi(Y, chi).
        Example
            idealToChi(Y, J)
    Caveat
        The function assumes that the input matrix Y and ideal J are defined over compatible rings.

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
            This function checks if one partition A is less than or equal to another partition B in the dominance order.
        Example
            A = {3, 2, 1}
            B = {4, 2, 1}
            partitionsLeq(A, B)
            partitionsLeq(B, A)
  
///



------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **TESTS** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------

TEST ///

///



end



