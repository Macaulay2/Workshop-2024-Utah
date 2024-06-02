

makeBatyrevMatrix = method();
makeBatyrevMatrix(List, List, List) := (P, B, C) -> (
    
    listRows := {};
    r10 := toList flatten(P#0 : {1});
    r11 := toList flatten(P#1 : {1});
    r12 := apply(C, i -> -i);
    r13 := apply(B, i -> -i - 1);
    r14 := toList flatten(P#4 : {0});

    R1 := r10 | r11 | r12 | r13 | r14;

    r20 := toList flatten(P#0 : {0});
    r21 := toList flatten(P#1 : {1});
    r22 := toList flatten(P#2 : {1});
    r23 := toList flatten(P#3 : {0});
    r24 := toList flatten(P#4 : {-1});

    R2 := r20 | r21 | r22 | r23 | r24;

    r30 := toList flatten(P#0 : {0});
    r31 := toList flatten(P#1 : {0});
    r32 := toList flatten(P#2 : {1});
    r33 := toList flatten(P#3 : {1});
    r34 := toList flatten(P#4 : {0});

    R3 := r30 | r31 | r32 | r33 | r34;

    r40 := toList flatten(P#0 : {0});
    r41 := toList flatten(P#1 : {-1});
    r42 := toList flatten(P#2 : {0});
    r43 := toList flatten(P#3 : {1});
    r44 := toList flatten(P#4 : {1});

    R4 := r40 | r41 | r42 | r43 | r44;

    r50 := toList flatten(P#0 : {1});
    r51 := toList flatten(P#1 : {0});
    r52 := apply(C, i -> -i);
    r53 := apply(B, i -> -i);
    r54 := toList flatten(P#4 : {1});
    
    R5 := r50 | r51 | r52 | r53 | r54;

    return matrix{R1, R2, R3, R4, R5}
)

batyrevRays = method();
batyrevRays(Matrix, List) := (A, P) -> (

    K := ker A;
    bRays := entries gens K;
    
    primitiveCollection1 := for i from 0 to P#0 + P#1 - 1 list i;
    primitiveCollection2 := for i from P#0 to P#0 + P#1 + P#2 -1 list i;
    primitiveCollection3 := for i from P#0 + P#1 to P#0 + P#1 + P#2 + P#3 - 1 list i;
    primitiveCollection4 := for i from P#0 + P#1 + P#2 to P#0 + P#1 + P#2 + P#3 + P#4 - 1 list i;
    primitiveCollection5a := for i from 0 to P#0 - 1 list i;
    primitiveCollection5b := for i from P#0 + P#1 + P#2 + P#3 to P#0 + P#1 + P#2 + P#3 + P#4 - 1 list i;
    primitiveCollection5 := primitiveCollection5a | primitiveCollection5b;

    return (bRays, {primitiveCollection1, primitiveCollection2, primitiveCollection3, primitiveCollection4, primitiveCollection5});
)

batyrevMaxCones = method();
batyrevMaxCones(List, List) := (bRays, primitiveCollections) -> (
    n := #bRays;
    x := symbol x;
    R := ZZ/2[x_0..x_(n-1)];
    Xs := first entries vars R;
    Ps := apply(primitiveCollections, i -> ideal apply(i, j -> x_j));
    I := intersect Ps;
    complementsOfCones := apply(first entries gens I, i -> support i);
    variableCones := apply(complementsOfCones, i -> toList(set Xs - set i));
    return apply(variableCones, i -> apply(i, j -> index j));

) 

batyrevConstructor = method();
batyrevConstructor(List, List, List) := (P, B, C) -> (
    A := makeBatyrevMatrix(P, B, C);
    (bRays, primitiveCollections) := batyrevRays(A, P);
    maxCones := batyrevMaxCones(bRays, primitiveCollections);
    return normalToricVariety(bRays, maxCones);
)


