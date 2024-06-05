loadPackage ("NumberFields", Reload=>true)
trace(RingElement) := (elt) -> (
    S = ring elt;
    return trace pushFwd(map(S^1, S^1, {{elt}}));
);
R = QQ[x]/(x^3-2)
trace(x)

