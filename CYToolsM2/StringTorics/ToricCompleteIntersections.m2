-- Need: better Hodge numbers, and also cohomologies (for complete intersections)
-- Also: can we be do better at getting the entire Picard group?

------------------------------------------------------------------------------------
-- link to Schubert2, as well as the ability to deal with complete intersections ---
------------------------------------------------------------------------------------

--------------------------------------------------------
-- Code for complete intersections in toric varieties --
--------------------------------------------------------
completeIntersection = method(Options => {
        Equations => true,
        Basis => null, -- a list of integer indices of rays that form a basis
        Variables => null -- variable names for each basis element
        })

completeIntersection(NormalToricVariety, List) := opts -> (Y,CIeqns) -> (
    if not all(CIeqns, d -> instance(d, ToricDivisor))
    then error "expected a list of toric divisors";
    if not all(CIeqns, d -> variety d === Y)
    then error "expected a list of toric divisors on the given toric variety";
    eqns := if opts.Equations then (
            S := ring Y;
            for D in CIeqns list random(degree D, S)
        ) else 
            null;
    B := if opts.Basis =!= null then (
        symbs := for i from 0 to #opts.Basis - 1 list opts.Variables_i;
        base toSequence symbs
        ); -- set to null if not being set
    X := new CompleteIntersectionInToric from {
        symbol Ambient => Y,
        symbol CI => CIeqns, -- these are the degrees
        symbol Equations => eqns,
        symbol Basis => opts.Basis,
        symbol Base => B,
        symbol cache => new CacheTable
        };
    X
    )


dim CompleteIntersectionInToric := (X) -> dim X.Ambient - #X.CI
ambient CompleteIntersectionInToric := (X) -> X.Ambient

equations CompleteIntersectionInToric := List => X -> (
    X.Equations
    )

lineBundle(CompleteIntersectionInToric, List) := (X, deg) -> (
    if not all(deg, x -> instance(x, ZZ)) or #deg =!= degreeLength ring ambient X
    then error("expected multidegree of length "|degreeLength ring ambient X);
    new LineBundle from {
        symbol cache => new CacheTable,
        symbol variety => X,
        symbol degree => deg
        }
    )

degree LineBundle := L -> L.degree
variety LineBundle := L -> L.variety

installMethod(symbol _, OO, CompleteIntersectionInToric, LineBundle => 
     (OO,X) -> lineBundle(X, (degree 1_(ring ambient X)))
     )

LineBundle Sequence := (L, deg) -> (
    lineBundle(variety L, degree L + toList deg)
    )


abstractVariety(CompleteIntersectionInToric, AbstractVariety) := opts -> (X,B) -> (
    if not X.cache#?(abstractVariety, B) then X.cache#(abstractVariety, B) = (
        aY := abstractVariety(ambient X, B);
        -- Question: how best to define F??
        bundles := X.CI/(d -> OO d);
        F := bundles#0;
        for i from 1 to #bundles-1 do F = F ++ bundles#i;
        aF := abstractSheaf(ambient X, B, F);
        sectionZeroLocus aF
        );
    X.cache#(abstractVariety, B)
    )

abstractVariety(CompleteIntersectionInToric) := opts -> (X) -> (
    if not X.cache#?(abstractVariety) then X.cache#(abstractVariety) = (
        aY := abstractVariety(ambient X, X.Base);
        -- Question: how best to define F??
        bundles := X.CI/(d -> OO d);
        F := bundles#0;
        for i from 1 to #bundles-1 do F = F ++ bundles#i;
        aF := abstractSheaf(ambient X, X.Base, F);
        Xa := sectionZeroLocus aF;
        X.cache.LinearForm = if X.Basis =!= null then (
            I := intersectionRing Xa;
            coeffsI := coefficientRing I;  -- this should be thevariables for the basis
            sum for i from 0 to #X.Basis - 1 list coeffsI_i * I_(X.Basis#i)
            );
        Xa
        );
    X.cache#(abstractVariety)
    )

linearForm = method()
linearForm CompleteIntersectionInToric := RingElement => X -> (
    Xa := abstractVariety X;
    X.cache.LinearForm
    )

intersectionRing CompleteIntersectionInToric := X -> (
    Xa := abstractVariety X;
    intersectionRing Xa
    )

intersectionForm = method()
intersectionForm CompleteIntersectionInToric := RingElement => X -> (
    if X.Basis === null then error "expected a basis to have been given";
    h := linearForm X;
    integral(h^(dim X))
    )

-- todo: this makes most sense for 3-folds...?
c2Form CompleteIntersectionInToric := RingElement => X -> (
    Xa := abstractVariety X;
    c2element := chern_2 tangentBundle Xa; -- I want a curve class here...
    h := linearForm X;
    integral(h^(dim X - first degree c2element) * c2element)
    )


-----------------------------
-- Place elsewhere ----------
-----------------------------
variety(CalabiYauInToric, Ring) := CompleteIntersectionInToric => (X, kk) -> (
    if not X.cache#?(variety, kk) then X.cache#(variety, kk) = (
        V := normalToricVariety X;
        X1 := completeIntersection(V, { - toricDivisor V});
        X1.cache.CalabiYauInToric = X;
        X1);
    X.cache#(variety, kk)
    )
variety CalabiYauInToric  := X -> variety(X, QQ)

-- CI/CalabiYauInToric compatibility exploration moved to ScratchToricCIs.m2
