loadPackage ("NumberFields", Reload=>true)
norm(RingElement) := (elt) ->(
    S = ring elt;
    return det pushFwd(map(S^1, S^1, {{elt}}));
);

