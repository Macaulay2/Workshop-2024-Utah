loadPackage ("NumberFields", Reload=>true)
trace(Ring, RingElement) := (S, elt) -> (
    trace pushFwd(map(S^1, S^1, {{elt}}))
);
