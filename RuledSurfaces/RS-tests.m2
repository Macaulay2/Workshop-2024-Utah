TEST ///
    X = Proj(QQ[x,y])
    F = (ring X)^2
    projectiveBundle(X,sheaf F)
    --TODO: add checks for variety, sheaf, rank of this bundle
    T = QQ[x]
    E = T^1
    assert instance(try projectiveBundle(X, sheaf E), Nothing)
///

TEST ///
    loadPackage "Divisor"
    loadPackage "RuledSurfaces"

    S = (ZZ/2)[x,y,z]
    I = ideal(x^3+y^3+z^3)
    R = S/I
    X = Proj R
    D0 = divisor x
    E = prune source first yonedaSheafExtension matrix( (Ext^1(OO_X,OO_X))_0)
    PE = projectiveBundle E
    q = homomorphism((Hom(E,OO_X^1))_0)
    sectionFromLineBundleQuotient(q)

    imageOfLinearSeries(PE, D0, 1)

    P = divisor ideal(y-11*z,x-33*z)
    D = 2*P
    E = OO_X^2
    PE = projectiveBundle E
    imageOfLinearSeries(PE, D0, 1)
    

///

TEST ///
    loadPackage "Divisor"
    loadPackage "RuledSurfaces"
    P2 = Proj ZZ/17[x,y,z]
    P1 = Proj ZZ/17[x,y]
    E = OO_P1^2
    E' = E(-1)
    PE = projectiveBundle E
    PE' = projectiveBundle E'
    assert instance(try imageOfLinearSeries(PE, OO_P2(1), 1), Nothing)
    imageOfLinearSeries(PE,OO_P1(1),1)
    imageOfLinearSeries(PE',OO_P1(2),1)
    q0=E^{0}
    q1=E^{1}
    
///
