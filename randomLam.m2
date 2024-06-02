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
