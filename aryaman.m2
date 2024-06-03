-- for my (aryaman's) personal use

partitionsLeq = method();
partitionsLeq(Partition, Partition) := (A, B) -> (
    -- A and B are partitions
    -- assuming weakly decreasing
    -- return if A <= B
    n := #B;
    for i in 0..#A-1 do(
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

c = apply(5, i -> randomLam(3, 2));
c = minimizeChi(c);
R = QQ[X_(1, 1)..X_(3, 4)];
XM = transpose genericMatrix(R, 4, 3);
I = idealIChi(XM, c);

n = 3;
m = 4;
K = QQ;
x = symbol x;
c = apply(5, i -> randomLam(3, 2));
(R, X) = myConstructor(n, m, K, x);
I = idealIChi(X, c);
J = I^2;
L = flatten entries mingens J;
gL = select(L, f -> goodDegree(n, m, degree f));
unique(apply(gL, degree))

myConstructor = (n, m, K, x) -> (
    -- n, m are ZZ, K is a ring, x is a variable/symbol?

    if not isField K then error("K must be a field");
    -- n is number of rows, m is number of columns
    -- want to degree to x_ij as (e_i, f_j) \in Z^n x Z^m = Z^(n + m)
    
    indices := toList((1, 1)..(n, m));
    myVars := apply(indices, ind -> x_(ind));
    edegree := i -> toList(insert(i-1, 1, ((n-1) : 0)));    -- e(i) = (0, ..., 0, 1, 0, ..., 0) \in Z^n
    fdegree := j -> toList(insert(j-1, 1, ((m-1) : 0)));    -- f(j) = (0, ..., 0, 1, 0, ..., 0) \in Z^m
    myDegrees := apply(indices, ind -> edegree(ind#0) | fdegree(ind#1));
    R := K[myVars, Degrees => myDegrees];
    X := transpose genericMatrix(R, m, n);

    return (R, X);
)

goodDegree = method();
goodDegree(ZZ, ZZ, List) := (n, m, L) -> (
    -- L is a length n + m list of integers
    -- return if L is a good degree for n, m
    -- i.e., if L = (lam, lam) for a partition lam (suitably padded)
    edegree := take(L, n);      -- the first n elements of L
    fdegree := drop(L, n);      -- the remaining elements
    
    isPartition = P -> (rsort(P) == P);
    
    if not (isPartition(edegree) and isPartition(fdegree)) then return false;
    -- edegree and fdegree are weakly decreasing

    edegree = delete(0, edegree);    -- remove the trailing 0s
    fdegree = delete(0, fdegree);    -- remove the trailing 0s
    return edegree == fdegree;
)



restart
load "GLIdeals.m2";
R = QQ[x_(1,1)..x_(2,3)];X = transpose genericMatrix (R, 3, 2);I = (minors(2, X))^3;
-- detLam(X, {3, 3})
idealToChi(X, I)
J = annihilator Ext^3(R^1/I, R)
idealToChi(X, J)