isTango = method()


planeTangoCurve = method()
planeTangoCurve(ZZ) := p -> (
    if not isPrime p_ZZ then error "expected a prime number p";
    X := symbol X;
    Y := symbol Y;
    Z := symbol Z;
    S:= (ZZ/p)[X,Y,Z];
    Proj quotient ideal(Y^p - Y*X^(p-1) - Z^(p-1)*X)
)

frobeniusSheafMap = method()
frobeniusSheafMap(ProjectiveVariety) := X ->(
    H := Hom(OO_X^1, frobeniusPushforward(1, OO_X));
    homomorphism(H_0)
)

--TODO make the below work for arbitrary e
planeTangoCurve(ZZ, ZZ) := (p,e) -> planeTangoCurve(p)

