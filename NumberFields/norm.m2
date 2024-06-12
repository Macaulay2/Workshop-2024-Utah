loadPackage ("NumberFields", Reload=>true)
norm(RingElement) := (elt) ->(
    S = ring elt;
    return det pushFwd(map(S^1, S^1, {{elt}}));
);
R = QQ[x]/(x^3-2)
norm(1_R)

