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
        Basis => null, -- a list of integer indices of rays taht form a basis
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

TEST ///
-*
  restart
*-
  debug needsPackage "StringTorics" -- remove 'debug'
  -- Let's consider smooth toric surfaces.
  for i from 0 to 4 list (
      V = smoothFanoToricVariety(2, 1);
      picardGroup V
      )
  V = smoothFanoToricVariety(2, 2)
  transpose matrix degrees ring V
  X = completeIntersection(V, {-toricDivisor V}, Basis => {0,3}, Variables => {symbol a, symbol b})
  linearForm X
  intersectionForm X
  linearForm X
  hh^* OO_X(0,0) -- X is an elliptic curve
  hh^*(OO_X(-1,1))

  L = OO_X(-1,1)
  assert(hh^0(L) == 1)
  assert(hh^1(L) == 0)
  assert(hh^* L == {1,0})
  intersectionForm X
  chern_1 tangentBundle abstractVariety X
  c2Form X
///

TEST ///
-*
restart
*-
  debug needsPackage "StringTorics" -- remove 'debug'

  -*    -- code to generate this example
    topes = kreuzerSkarke(3, Limit => 50);    
    A = matrix topes_30
    P = convexHull A
    (V, basisElems) = reflexiveToSimplicialToricVarietyCleanDegrees(P, CoefficientRing => ZZ/32003)
  *-
  
  verts = {{-1, -1, 0, 0}, {-1, -1, 0, 1}, {-1, -1, 2, 0}, {-1, 0, 0, 0}, {1, -1, -1, 1}, {1, 2, -1, -1}, {-1, -1, 1, 0}}
  maxcones = {{0, 1, 3, 4}, {0, 1, 3, 6}, {0, 1, 4, 6}, {0, 3, 4, 5}, {0, 3, 5, 6}, {0, 4, 5, 6}, {1, 2, 3, 5}, {1, 2, 3, 6}, {1, 2, 4, 5}, {1, 2, 4, 6}, {1, 3, 4, 5}, {2, 3, 5, 6}, {2, 4, 5, 6}}
  V = normalToricVariety(verts, maxcones, CoefficientRing => ZZ/32003)
  glsm = transpose matrix degrees ring V
  basiselems = {0, 5, 6}
  X = completeIntersection(V, {-toricDivisor V}, Basis => {0, 5, 6}, Variables => {symbol a, symbol b, symbol c})
  assert(dim X == 3)
  intersectionRing X
  h = linearForm X
  integral(h^3)
  assert(intersectionForm X == a^3-6*a^2*b+6*a*b^2+12*b^3-9*a^2*c+24*a*b*c+21*a*c^2-24*b*c^2-16*c^3)
  assert(c2Form X == 10*a + 60*b + 8*c)

  X = completeIntersection(V, {-toricDivisor V}, Basis => {0, 5, 6}, Variables => symbol a)

  D = completeIntersection(V, {-toricDivisor V, V_0}, Basis => {0, 5, 6}, Variables => {symbol a, symbol b, symbol c})
  dim D == 2
  saturate(ideal equations D, ideal V)
  intersectionRing D
  intersectionForm D
  c2Form D
  chern_2 tangentBundle abstractVariety D

  L = OO_X(1,1,2)
  assert(hh^* L == {8, 17, 0, 0}) -- TODO: recheck these numbers!
  assert(variety L === X)
  assert(degree L == {1,1,2})

  -- TODO: ADD BACK IN once DanilovKhovanskii more functional
  -- needsPackage "DanilovKhovanskii"
  -- computeHodgeDeligne(-toricDivisor V)
  -- oo#1
  -- matrix for i from 0 to dim X list for j from 0 to dim X list (-1)^(i+j) * oo#(i,j)

  -- want to be able to turn this into a CalabiYauInToric...
  -- then we can check computations against each other too.
  
///

TEST ///
-*
restart
*-
  debug needsPackage "StringTorics" -- remove 'debug'
  V = kleinschmidt(3, {2,1}, CoefficientRing => ZZ/101)
  rays V
  max V
  picardGroup V
  isSmooth V
  transpose matrix degrees ring V
  X = completeIntersection(V, {-toricDivisor V}, Basis => {3,0}, Variables => {a,b})
  dim X == 2
  hh^*(OO_X(0,0)) == {1, 0, 1}
  F = first equations X
  saturate(ideal F + ideal jacobian F, ideal V) -- X is smooth
  -- hh^*(OO_V(-1,3)) -- ouch!  needs to work...
///


"TEST"
///
-- DanilovKhovanskii
-*
restart
*-
  debug needsPackage "StringTorics" -- remove 'debug'
  needsPackage "DanilovKhovanskii"
  V = kleinschmidt(3, {2,1}, CoefficientRing => ZZ/101)
  rays V
  max V
  picardGroup V
  isSmooth V
  transpose matrix degrees ring V
  X = completeIntersection(V, {-toricDivisor V}, Basis => {3,0}, Variables => {a,b})
  dim X == 2
  hh^*(OO_X(0,0)) == {1, 0, 1}
  F = first equations X
  saturate(ideal F + ideal jacobian F, ideal V) -- X is smooth
  -- hh^*(OO_V(-1,3)) -- ouch!  needs to work...

  computeHodgeDeligne(-toricDivisor V)
  oo#1
  matrix for i from 0 to dim X list for j from 0 to dim X list oo#(i,j)

  
///

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

///
 -- how compatible with CalabiYauInToric is this?
  -- I guess we need to know if the triangulation comes from a triangulation of the polytope...
  
-*
restart
*-
  debug needsPackage "StringTorics" -- remove 'debug'

  -*    -- code to generate this example
    topes = kreuzerSkarke(3, Limit => 50);    
    A = matrix topes_30
    P = convexHull A
    (V, basisElems) = reflexiveToSimplicialToricVarietyCleanDegrees(P, CoefficientRing => ZZ/32003)
  *-
  
  verts = {{-1, -1, 0, 0}, {-1, -1, 0, 1}, {-1, -1, 2, 0}, {-1, 0, 0, 0}, {1, -1, -1, 1}, {1, 2, -1, -1}, {-1, -1, 1, 0}}
  maxcones = {{0, 1, 3, 4}, {0, 1, 3, 6}, {0, 1, 4, 6}, {0, 3, 4, 5}, {0, 3, 5, 6}, {0, 4, 5, 6}, {1, 2, 3, 5}, {1, 2, 3, 6}, {1, 2, 4, 5}, {1, 2, 4, 6}, {1, 3, 4, 5}, {2, 3, 5, 6}, {2, 4, 5, 6}}
  V = normalToricVariety(verts, maxcones, CoefficientRing => ZZ/32003)
  glsm = transpose matrix degrees ring V
  basiselems = {0, 5, 6}
  X = completeIntersection(V, {-toricDivisor V}, Basis => {0, 5, 6}, Variables => {symbol a, symbol b, symbol c})

  -- we need a function that determines if this is a Batryev CY3.
  -- and to return the corresponding CalabiYauInToric...
  transpose matrix rays V
  P = convexHull oo
  vertices P
  Q = cyPolytope P
  rays Q
  X1 = first findAllCYs Q
  hodgeDiamond X1
  hodgeDiamond X
  hh^(1,2) X1 == 69

  -- TODO: add back in once DanilovKhovanskii is in the system
  -- needsPackage "DanilovKhovanskii"
  -- computeHodgeDeligne (-toricDivisor V) -- this is not the correct answer I think!
  -- oo#1
  -- matrix for i from 0 to dim X list for j from 0 to dim X list (-1)^(i+j) * oo#(i,j)
  -- assert(hh^(1,2) Q == 69)
  -- assert(hh^(1,1) Q == 3)
///
