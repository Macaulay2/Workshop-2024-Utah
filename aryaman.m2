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
     print ("expected: ",d*e);
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