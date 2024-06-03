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