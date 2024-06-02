-- for my (aryaman's) personal use

partitionsLeq = method();
partitionsLeq(Partition, Partition) := (A, B) -> (
    -- A and B are partitions
    -- assuming weakly decreasing
    -- return if A <= B
    n := #B;
    for i in 0..#A-1 do(
        if A#i == 0 then break;
        if i > n then return false;
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