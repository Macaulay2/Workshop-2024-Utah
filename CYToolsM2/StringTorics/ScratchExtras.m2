-- Scratch code from Extras.m2
-- Unused decompose function for simplicial complexes

-*
decompose SimplicialComplex := C -> (
    -- return a list of sub-complexes corresponding to the
    -- connected components of C
    S := ring C;
    V := faces(0,C);
    G := graph(V, for e in faces(1,C) list support e);
    comps := connectedComponents G;
    for comp in comps list subcomplex(C,comp)
    )
*-

///
  S = QQ[x_0..x_6]
  C = simplicialComplex {x_0*x_1, x_2*x_3*x_4, x_0*x_2, x_5*x_6}
  decompose C
///

--    fac := flatten entries facets C;
--    fac = fac/support/(v -> v/index//sort);
--    faces := unique for f in fac list sort toList (set D * set f);
--    -- Next we take these faces, and make monomials out of them.
--    simplicialComplex flatten for f in faces list (1_S * product (f/(i -> S_i)))
--    )
