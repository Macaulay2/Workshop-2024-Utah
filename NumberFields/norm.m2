loadPackage ("NumberFields", Reload=>true)
norm(Ring, RingElement) := (S, elt) ->(
    det pushFwd(map(S^1, S^1, {{elt}}))
);
