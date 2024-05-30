------------------------------------------------------------
-- This file contains various invariants coming from the  --
-- cubic intersection form, and c2 form. -------------------
------------------------------------------------------------
-- Invariants include:
--  invariantsHubsch -- Hubsch invariants of L, F.

--  1. simple ones: content(L),content(F) (h11, h12 are either passed in, or not considered).
--  2. Hubsch invariants
--  3. Point counts
--  4. Shape of (a ZZ-factorization of) the Hessian (and its content).
--  5. Singularity info the cubic
--  6. Singularity info of the intersection of the linear form and cubic.
--  7. values of invariants (at leat for h11=3, maybe h11=4, what about higher?
-- Each invariants type can be called with (L, F, h11, h12), (L, F), or a CalabiYauInToric, or CYToolsCY3
-- Point counting though requires a PointCount class in addition (otherwise, too slow).

---------------------------------------
-- Hubsch invariants ------------------
---------------------------------------
hubsch1 = X -> (
    -- X is a CalabYauInToric
    N := hh^(1,1) X;
    i := hashTable intersectionNumbers X;
    gcdAll := gcd toSequence values i;
    gcdDistincts := gcd toSequence(
      flatten for A from 0 to N-3 list 
      flatten for B from A+1 to N-2 list 
      for C from B+1 to N-1 list (
        triple := sort {A,B,C};
        if i#?triple then i#triple else 0
        ));
    gcdPairs := gcd toSequence(
        -- this includes singles
        flatten for A from 0 to N-1 list for B from 0 to N-1 list (
            triple := sort {A,A,B};
            if i#?triple then i#triple else 0
            ));
    gcdPairs2 := gcd toSequence(
        -- this includes singles
        flatten for A from 0 to N-1 list flatten for B from 0 to N-1 list (
            triple1 := sort {A,A,B};
            triple2 := sort {A,B,B};
            val1 := if i#?triple1 then i#triple1 else 0;
            val2 := if i#?triple2 then i#triple2 else 0;
            {val1 + val2, val1 - val2}
            ));
    gcdSingles := gcd toSequence(
        for A from 0 to N-1 list (
            triple := sort {A,A,A};
            if i#?triple then i#triple else 0
        ));
    --return (gcdAll, gcdDistincts, gcdPairs, gcdSingles, gcdPairs);
    d1 := gcdAll;
    d2 := gcd(gcdPairs, 2*gcdAll);
    d3 := gcd(gcdSingles, 3*gcdPairs2, 6*gcdAll);
    (d1, d2, d3)
    )

quadlinear = X -> (
    N := hh^(1,1) X;
    i := hashTable intersectionNumbers X;
    L := -2 * c2 X; -- p1
    ifcn := triple -> (triple = sort triple; if i#?triple then i#triple else 0);
    unique for abcd in (0,0,0,0)..(N-1,N-1,N-1,N-1) list (
        (a,b,c,d) := abcd;
        val := ifcn {a,b,c} * L#d + ifcn {b,c,d} * L#a + ifcn {c,d,a} * L#b + ifcn {d,a,b} * L#c;
        if val != 0 then sort{a,b,c,d} => val else continue
        )
    )

hubsch2 = X -> (
    N := hh^(1,1) X;
    H := hashTable quadlinear X;
    hfcn := ind -> (ind = sort ind; if H#?ind then H#ind else 0);
    d4 := gcd(values H);
    onedouble := gcd for k in keys H list if #unique k <= 3 then H#k else continue;
    d5 := gcd(onedouble, 2*d4);
    triples := for k in keys H list if any(values tally k, val -> val >= 3) then k else continue;
    onetriple := if #triples > 0 then gcd for k in triples list H#k else 0;
    d6part2 := toSequence for acd in (0,0,0)..(N-1,N-1,N-1) list (
        (a,c,d) := acd;
        gcd(hfcn {a,a,c,d} + hfcn {a,c,c,d}, hfcn {a,a,c,d} - hfcn {a,c,c,d})
        );
    d6part2 = gcd d6part2;
    d6 := gcd(onetriple, 3*d6part2, 6*d4);
    -- Now for d7...
    d7part1 := gcd for a from 0 to N-1 list hfcn {a,a,a,a};
    --  2(2<a,a,a,d> \pm 3 <a,a,d,d> \pm <a,d,d,d>).
    d7part2s := for ad in (0,0)..(N-1,N-1) list (
        (a,d) := ad;
        (a,d) => 2 * {2 * hfcn {a,a,a,d} + 3 * hfcn {a,a,d,d} + hfcn {a,d,d,d},
            2 * hfcn {a,a,a,d} + 3 * hfcn {a,a,d,d} - hfcn {a,d,d,d},
            2 * hfcn {a,a,a,d} - 3 * hfcn {a,a,d,d} + hfcn {a,d,d,d},
            2 * hfcn {a,a,a,d} - 3 * hfcn {a,a,d,d} - hfcn {a,d,d,d}}
        );
    d7part2 := gcd for ad in (0,0)..(N-1,N-1) list (
        (a,d) := ad;
        2 * gcd(2 * hfcn {a,a,a,d} + 3 * hfcn {a,a,d,d} + 2 * hfcn {a,d,d,d},
            2 * hfcn {a,a,a,d} + 3 * hfcn {a,a,d,d} - 2 * hfcn {a,d,d,d},
            2 * hfcn {a,a,a,d} - 3 * hfcn {a,a,d,d} + 2 * hfcn {a,d,d,d},
            2 * hfcn {a,a,a,d} - 3 * hfcn {a,a,d,d} - 2 * hfcn {a,d,d,d})
        );
    d7part3 := gcd for acd in (0,0,0)..(N-1,N-1,N-1) list (
        (a,c,d) := acd;
        12 * gcd(
            hfcn {a,a,c,d} + hfcn {a,c,c,d} + hfcn {a,c,d,d},
            hfcn {a,a,c,d} + hfcn {a,c,c,d} - hfcn {a,c,d,d},
            hfcn {a,a,c,d} - hfcn {a,c,c,d} + hfcn {a,c,d,d},
            hfcn {a,a,c,d} - hfcn {a,c,c,d} - hfcn {a,c,d,d})
        );
    --<< (d7part1, d7part2, d7part3, 24*d4) << endl;
    d7 := gcd(d7part1, d7part2, d7part3, 24*d4);
    (d4, d5, d6, d7)
    )

hubschInvariants = method()
hubschInvariants CalabiYauInToric := X -> join(hubsch1 X, hubsch2 X, {gcd c2 X})

---------------------------------------
-- Point counts -----------------------
---------------------------------------
allPrimitivePoints = (ht, n) -> (
    -- all nonzero primitive integer points in n-space, whose first
    -- nonzero value is > 0
    -- and whose gcd is 1, with entries a satisfying |a| <= ht
    )

allPointsToHeight = (ht, n) -> (
    -- all points in ZZ^n, with each entry a s.t. |a| <= ht.
    pts := for a from -ht to ht list {a};
    if n === 1 then return pts;
    if n === 0 then return {{}};
    if n < 0 then error "internal logic error";
    b := allPointsToHeight(ht, n-1);
    flatten for a from -ht to ht list (b/(b1 -> prepend(a, b1)))
    )

allPrimitivePointsToHeight = (ht, n) -> (
    pts := allPointsToHeight(ht, n);
    for p in pts list (
        nonzero := select(1, p, a ->  a != 0);
        if #nonzero == 0 then continue;
        --if nonzero#0 < 0 then continue;
        if gcd p != 1 then continue;
        p
        )
    )

allPoints = (p, n) -> (
    -- all points in kk = ZZ//p in kk^n
    pts := for a from 0 to p-1 list {a};
    if n === 1 then return pts;
    if n === 0 then return {{}};
    if n < 0 then error "internal logic error";
    b := allPoints(p, n-1);
    flatten for a from 0 to p-1 list (b/(b1 -> prepend(a, b1)))
    )

allProjectivePoints = (p, n) -> (
    -- all points in n-space (with p elements in the field) with first non-zero value = 1, except 0.
    flatten for i from 1 to n list (
        -- collect all the points with first non-zero value at location i+1.
        firstPart := splice{(i-1): 0, 1};
        --print firstPart;
        for pt in allPoints(p, n-i) list join(firstPart, pt)
        ))

createPointMaps = method(Options => {Projective => true})
createPointMaps(ZZ, Ring) := opts -> (p, R) -> (
    pts := if opts.Projective then allProjectivePoints(p, numgens R) else allPoints(p, numgens R);
    kk := ZZ/p;
    for pt in pts list map(kk, R, pt)
    )
createPointMaps(Sequence, Ring) := opts -> (pr, R) -> (
    (p,r) := pr;
    q := p^r;
    pts0 := if opts.Projective then allProjectivePoints(q, numgens R) else allPoints(q, numgens R);
    -- now for each one we translate to GF(p,nr);
    t := local t;
    kk := GF(p^r, Variable => t);
    print kk_0;
    H := new MutableList;
    H#0 = 0_kk;
    H#1 = kk_0; -- the variable;
    for i from 2 to q-1 do H#i = kk_0 * H#(i-1);
    for pt0 in pts0 list (
        pt := for a in pt0 list H#a;
        map(kk, R, pt)
        )
    )

PointCounter = new Type of HashTable
pointCounter = method(Options=> {
--        "Primes" => {2,3,5,7,11,13,(2,2),(3,2),(2,3),(2,4)}, -- h11=3 gives 166 diff.
--        "Primes" => {2,3,5,7,11,13,(2,2),(3,2),(2,3),17}, -- h11=3 gives 167 diff
          "Primes" => {2,3,5,7,11,13,(2,2),(3,2),(2,3),17,19}, -- h11=3 gives 168 diff
--          "Primes" => {2,3,5,7,11,13,(2,2),(3,2),(2,3),17,19,23,29,31,37}, -- h11=3 still gives 168 diff
--        "Primes" => {2,3,5,7,11,13,(2,2),(3,2),(2,3),(2,4),17,19,23,(5,2),(3,3),29},
        Projective => true
        })
pointCounter Ring := PointCounter => opts -> RZ -> (
    -- we expect that RZ is a polynomial ring over ZZ, in n variables.
    -- each element of LprimePowers is a prime or a sequence (p,r).
    -- In the latter case we find points over GF(p^r), in the former, over ZZ/p.
    n := numgens RZ;
    LprimePowers := opts#"Primes";
    Lmaps := for pr in LprimePowers list createPointMaps(pr, RZ, Projective => opts.Projective);
    PC := new PointCounter from {
        symbol Ring => RZ,
        "Primes" => LprimePowers,
        "Elements" => Lmaps
        };
    PC
    )

pointCounts = method()
pointCounts(PointCounter, RingElement, RingElement) := List => (PC, L, F) -> (
    -- We return a list with 3 lists in it:
    cL := polynomialContent L;
    cF := polynomialContent F;
    if cL != 1 then (
        L = L // cL;
        --<< "content L = " << cL << " and L/cL = " << L << endl;
        );
    if cF != 1 then (
        F = F // cF;
        --<< "content F = " << cF << " and F/cF = " << F << endl;
        );
    R := PC.Ring;
    if R =!= ring L or R =!= ring F then error "expected elements to be over the same ring";
    elems := PC#"Elements";
    for eachset in elems list (
        nzerosL := 0;
        nzerosF := 0;
        nzerosBothLandF := 0;
        for phi in eachset list (
            a := phi L;
            b := phi F;
            if a == 0 then nzerosL = nzerosL + 1;
            if b == 0 then nzerosF = nzerosF + 1;
            if a == 0 and b ==0 then nzerosBothLandF = nzerosBothLandF + 1;
            );
        {nzerosL, nzerosF, nzerosBothLandF}
        )
    )
pointCounts(PointCounter, CalabiYauInToric) := List => (PC, X) -> (
    elapsedTime pointCounts(PC, c2Form X, cubicForm X)
    )
TEST ///
-*
 restart
 debug needsPackage "StringTorics"
*-
  debug StringTorics -- for allPoints, createPointMaps
  assert(
      allPoints(3, 2) 
      === 
      {{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}}
      )

  assert(allPoints(3,1) === {{0}, {1}, {2}})

  allPoints(4, 3) 

  R = ZZ[a,b]
  createPointMaps(3, R)
  createPointMaps(3, R, Projective => false)
  createPointMaps((2,2), R)
  createPointMaps((2,2), R, Projective => false)
  createPointMaps((2,3), ZZ[a,b,c])

  R = ZZ[a,b,c];
  elapsedTime PC = pointCounter(R, "Primes" => {2,3,5,7,11,13,(2,2),(3,2),(2,3)});
  transpose matrix pointCounts(PC, a+b, a^3+b^3+c^3-3*a*b*c)
  transpose matrix pointCounts(PC, a+b, a^3+b^3+c^3-4*a*b*c)
///

invariantsH11H12 = method()
invariantsH11H12 CalabiYauInToric := X -> List => {hh^(1,1) X, hh^(1,2) X}

invariantsContents = method()
invariantsContents CalabiYauInToric := List => X -> {
    polynomialContent c2Form X,
    polynomialContent cubicForm X
    }

-----------------------------
-- Hessian invariants -------
-----------------------------
hessianInvariants = method()
hessianInvariants CalabiYauInToric := X -> (
    factorShape det hessian cubicForm X
    )

-----------------------------
-- Singular set invariants --
-----------------------------

linearcontent := (I) -> (
    if I == 0 then return 0;
    lins := select(I_*, f -> f != 0 and first degree f <= 1);
    if #lins == 0 then return 0;
    gcd for ell in lins list (trim content ell)_0
    );

contentToDegree := (d, I) -> (
    if I == 0 then return 0;
    lins := select(I_*, f -> f != 0 and first degree f <= d);
    if #lins == 0 then return 0;
    gcd for ell in lins list (trim content ell)_0
    );

cubicConductorInvariants = method()
cubicConductorInvariants CalabiYauInToric := X -> (
    F := cubicForm X;
    jac := ideal F + ideal jacobian F;
    jacsat := saturate jac;
    {integerPart jacsat, linearcontent jacsat}
    )

cubicLinearConductorInvariants = method()
cubicLinearConductorInvariants CalabiYauInToric := X -> (
    L := c2Form X;
    F := cubicForm X;
    I := ideal(L,F);
    jac := I + minors(2, jacobian I);
    jacsat := saturate jac;
    integerPart jacsat
    )

singularContents = method()
singularContents CalabiYauInToric := (X) -> (
    F := cubicForm X;
    F = F // polynomialContent F;
    jac := ideal F + ideal jacobian F;
    jacsat := saturate jac;
    H := partition(f -> first degree f, jacsat_*);
    degs := sort keys H;
    prevgcd := 0;
    done := false;
    for i from 0 to max degs list (
        if prevgcd == 1 then break;
        if not H#?i then prevgcd
        else (
            gcd1 := gcd((H#i)/polynomialContent);
            prevgcd = gcd(gcd1, prevgcd);
            prevgcd
            ))
    )

singularContents Ideal := (J) -> (
    if not isHomogeneous J then error "expected homogeneous ideal";
    -- Also expect: coefficients are ZZ.
    -- Grading is singly graded.
    H := partition(f -> first degree f, J_*);
    degs := sort keys H;
    prevgcd := 0;
    done := false;
    for i from 0 to max degs list (
        if prevgcd == 1 then break;
        if not H#?i then prevgcd
        else (
            gcd1 := gcd((H#i)/polynomialContent);
            prevgcd = gcd(gcd1, prevgcd);
            prevgcd
            ))
    )

singularContentsQuartic = method()
singularContentsQuartic CalabiYauInToric := (X) -> (
   << "singularContentsQuartic: doing " << label X << endl;
   F := cubicForm X;
   F = F // polynomialContent F;
   L := c2Form X;
   L = L // polynomialContent L;
   F = F*L;
   jac := ideal F + ideal jacobian F;
   jacsat := saturate jac;
   singularContents jacsat
   )
 
----------------------------------------------------------------------------------

--Ternary cubic form for
-- a*x^3 + b*y^3 + c*z^3 + 3*d*x^2*y + 3*e*y^2*z + 3*f*z^2*x + 3*g*x*y^2 + 
--   3*h*y*z^2 + 3*i*z*x^2 + 6*j*x*y*z

aronholdS = method()
aronholdS RingElement := F -> (
    -- F should be a cubic in 3 variables.
    S := ring F;
    KK := coefficientRing S;
    if numgens S =!= 3 then error "expected polynomial in 3 variables";
    if first degree F != 3 then error "expected cubic polynomial";
    t := local t;
    S1 := KK[t_(1,0)..t_(4,2)];
    DM := genericMatrix(S1, 3, 4);
    Fproduct := product for i from 1 to 4 list sub(F, {S_0 => t_(i,0), S_1 => t_(i,1), S_2 => t_(i,2)});
    f1 := diff((det DM_{3,0,1}), Fproduct);
    f2 := diff((det DM_{2,3,0}), f1);
    f3 := diff((det DM_{1,2,3}), f2);
    f4 := diff((det DM_{0,1,2}), f3);
    ans := 1/31104 * f4;
    lift(ans, KK)
    )

aronholdT = method()
aronholdT RingElement := F -> (
    -- F should be a cubic in 3 variables.
    S := ring F;
    KK := coefficientRing S;
    if numgens S =!= 3 then error "expected polynomial in 3 variables";
    if first degree F != 3 then error "expected cubic polynomial";
    t := local t;
    S1 := KK[t_(1,0)..t_(6,2)];
    DM := genericMatrix(S1, 3, 6);
    G := product for i from 1 to 6 list sub(F, {S_0 => t_(i,0), S_1 => t_(i,1), S_2 => t_(i,2)});
    g1 := elapsedTime diff((det DM_{3,4,5}), G);
    g2 := elapsedTime diff((det DM_{3,4,5}), g1);
    g3 := elapsedTime diff((det DM_{2,0,5}), g2);
    g4 := elapsedTime diff((det DM_{1,2,4}), g3);
    g5 := elapsedTime diff((det DM_{0,1,3}), g4);
    g6:= elapsedTime diff((det DM_{0,1,2}), g5);
    lift(1/279936 * g6, KK)
    )

-- j-invariant is: j = 1728*(4S)^3/((4S)^3 - T^2) (from https://www.notzeb.com/aronhold.html)
-- still need to check that.

aronhold = method()
aronhold RingElement := List => F -> (
    -- formulas taken from https://www.notzeb.com/aronhold.html, checked by a test in this file.
    H := hashTable toCOO F;
    hf := a -> if H#?a then H#a else 0;
    -- F = a*x^3 + b*y^3 + c*z^3 + 3*d*x^2*y + 3*e*y^2*z + 3*f*z^2*x + 
    --     3*g*x*y^2 + 3*h*y*z^2 + 3*i*z*x^2 + 6*j*x*y*z
    a := hf {0,0,0}; -- x^3
    b := hf {1,1,1}; -- y^3
    c := hf {2,2,2};
    d := hf {0,0,1}; -- x^2*y
    e := hf {1,1,2}; -- y^2*z
    f := hf {0,2,2}; -- z^2*x
    g := hf {0,1,1};
    h := hf {1,2,2};
    i := hf {0,0,2};
    j := hf {0,1,2};
    S := a*g*e*c - a*g*h^2 - a*j*b*c + a*j*e*h + a*f*b*h - a*f*e^2 - 
      d^2*e*c + d^2*h^2 + d*i*b*c - d*i*e*h + d*g*j*c - d*g*f*h - 
      2*d*j^2*h + 3*d*j*f*e - d*f^2*b - i^2*b*h + i^2*e^2 - 
      i*g^2*c + 3*i*g*j*h - i*g*f*e - 2*i*j^2*e + i*j*f*b + 
      g^2*f^2 - 2*g*j^2*f + j^4;
    T := a^2*b^2*c^2 - 3*a^2*e^2*h^2 - 6*a^2*b*e*h*c + 4*a^2*b*h^3 + 4*a^2*e^3*c - 
    6*a*d*g*b*c^2 + 18*a*d*g*e*h*c - 12*a*d*g*h^3 + 12*a*d*j*b*h*c - 24*a*d*j*e^2*c + 
    12*a*d*j*e*h^2 - 12*a*d*f*b*h^2 + 6*a*d*f*b*e*c + 6*a*d*f*e^2*h + 
    6*a*i*g*b*h*c - 12*a*i*g*e^2*c + 6*a*i*g*e*h^2 + 12*a*i*j*b*e*c + 
    12*a*i*j*e^2*h - 6*a*i*f*b^2*c + 18*a*i*f*b*e*h - 24*a*g^2*j*h*c - 
    24*a*i*j*b*h^2 - 12*a*i*f*e^3 + 4*a*g^3*c^2 - 12*a*g^2*f*e*c + 
    24*a*g^2*f*h^2 + 36*a*g*j^2*e*c + 12*a*g*j^2*h^2 + 12*a*g*j*f*b*c - 
    60*a*g*j*f*e*h - 12*a*g*f^2*b*h + 24*a*g*f^2*e^2 - 20*a*j^3*b*c - 
    12*a*j^3*e*h + 36*a*j^2*f*b*h + 12*a*j^2*f*e^2 - 24*a*j*f^2*b*e + 
    4*a*f^3*b^2 + 4*d^3*b*c^2 - 12*d^3*e*h*c + 8*d^3*h^3 + 24*d^2*i*e^2*c - 
    12*d^2*i*e*h^2 + 12*d^2*g*j*h*c + 6*d^2*g*f*e*c - 24*d^2*j^2*h^2 - 
    12*d^2*i*b*h*c - 3*d^2*g^2*c^2 - 24*g^2*j^2*f^2 + 24*g*j^4*f - 
    12*d^2*g*f*h^2 + 12*d^2*j^2*e*c - 24*d^2*j*f*b*c - 27*d^2*f^2*e^2 + 
    36*d^2*j*f*e*h + 24*d^2*f^2*b*h + 24*d*i^2*b*h^2 - 12*d*i^2*b*e*c - 
    12*d*i^2*e^2*h + 6*d*i*g^2*h*c - 60*d*i*g*j*e*c + 36*d*i*g*j*h^2 + 
    18*d*i*g*f*b*c - 6*d*i*g*f*e*h + 36*d*i*j^2*b*c - 12*d*i*j^2*e*h - 
    60*d*i*j*f*b*h + 36*d*i*j*f*e^2 + 6*d*i*f^2*b*e + 12*d*g^2*j*f*c - 
    12*d*g*j^3*c - 12*d*g*j^2*f*h + 36*d*g*j*f^2*e - 12*d*g*f^3*b + 
    24*d*j^4*h + 12*d*j^2*f^2*b + 4*i^3*b^2*c + 24*i^2*g^2*e*c - 
    27*i^2*g^2*h^2 - 36*d*j^3*f*e - 12*i^3*b*e*h + 8*i^3*e^3 - 24*i^2*g*j*b*c + 
    36*i^2*g*j*e*h + 6*i^2*g*f*b*h + 12*i^2*j^2*b*h - 3*i^2*f^2*b^2 - 
    12*d*g^2*f^2*h - 12*i^2*g*f*e^2 - 24*i^2*j^2*e^2 + 12*i^2*j*f*b*e - 
    12*i*g^3*f*c + 12*i*g^2*j^2*c + 36*i*g^2*j*f*h - 12*i*g^2*f^2*e - 
    36*i*g*j^3*h - 12*i*g*j^2*f*e + 12*i*g*j*f^2*b + 24*i*j^4*e - 
    12*i*j^3*f*b + 8*g^3*f^3 - 8*j^6;
    {S, T}
    )
    
/// -- TEST: takes too long
  debug StringTorics -- for aronholdS, aronholdT.
  KK = QQ[a,b,c,d,e,f,g,h,i,j]

  S = a*g*e*c - a*g*h^2 - a*j*b*c + a*j*e*h + a*f*b*h - a*f*e^2 - 
      d^2*e*c + d^2*h^2 + d*i*b*c - d*i*e*h + d*g*j*c - d*g*f*h - 
      2*d*j^2*h + 3*d*j*f*e - d*f^2*b - i^2*b*h + i^2*e^2 - 
      i*g^2*c + 3*i*g*j*h - i*g*f*e - 2*i*j^2*e + i*j*f*b + 
      g^2*f^2 - 2*g*j^2*f + j^4
      
  T = a^2*b^2*c^2 - 3*a^2*e^2*h^2 - 6*a^2*b*e*h*c + 4*a^2*b*h^3 + 4*a^2*e^3*c - 
    6*a*d*g*b*c^2 + 18*a*d*g*e*h*c - 12*a*d*g*h^3 + 12*a*d*j*b*h*c - 24*a*d*j*e^2*c + 
    12*a*d*j*e*h^2 - 12*a*d*f*b*h^2 + 6*a*d*f*b*e*c + 6*a*d*f*e^2*h + 
    6*a*i*g*b*h*c - 12*a*i*g*e^2*c + 6*a*i*g*e*h^2 + 12*a*i*j*b*e*c + 
    12*a*i*j*e^2*h - 6*a*i*f*b^2*c + 18*a*i*f*b*e*h - 24*a*g^2*j*h*c - 
    24*a*i*j*b*h^2 - 12*a*i*f*e^3 + 4*a*g^3*c^2 - 12*a*g^2*f*e*c + 
    24*a*g^2*f*h^2 + 36*a*g*j^2*e*c + 12*a*g*j^2*h^2 + 12*a*g*j*f*b*c - 
    60*a*g*j*f*e*h - 12*a*g*f^2*b*h + 24*a*g*f^2*e^2 - 20*a*j^3*b*c - 
    12*a*j^3*e*h + 36*a*j^2*f*b*h + 12*a*j^2*f*e^2 - 24*a*j*f^2*b*e + 
    4*a*f^3*b^2 + 4*d^3*b*c^2 - 12*d^3*e*h*c + 8*d^3*h^3 + 24*d^2*i*e^2*c - 
    12*d^2*i*e*h^2 + 12*d^2*g*j*h*c + 6*d^2*g*f*e*c - 24*d^2*j^2*h^2 - 
    12*d^2*i*b*h*c - 3*d^2*g^2*c^2 - 24*g^2*j^2*f^2 + 24*g*j^4*f - 
    12*d^2*g*f*h^2 + 12*d^2*j^2*e*c - 24*d^2*j*f*b*c - 27*d^2*f^2*e^2 + 
    36*d^2*j*f*e*h + 24*d^2*f^2*b*h + 24*d*i^2*b*h^2 - 12*d*i^2*b*e*c - 
    12*d*i^2*e^2*h + 6*d*i*g^2*h*c - 60*d*i*g*j*e*c + 36*d*i*g*j*h^2 + 
    18*d*i*g*f*b*c - 6*d*i*g*f*e*h + 36*d*i*j^2*b*c - 12*d*i*j^2*e*h - 
    60*d*i*j*f*b*h + 36*d*i*j*f*e^2 + 6*d*i*f^2*b*e + 12*d*g^2*j*f*c - 
    12*d*g*j^3*c - 12*d*g*j^2*f*h + 36*d*g*j*f^2*e - 12*d*g*f^3*b + 
    24*d*j^4*h + 12*d*j^2*f^2*b + 4*i^3*b^2*c + 24*i^2*g^2*e*c - 
    27*i^2*g^2*h^2 - 36*d*j^3*f*e - 12*i^3*b*e*h + 8*i^3*e^3 - 24*i^2*g*j*b*c + 
    36*i^2*g*j*e*h + 6*i^2*g*f*b*h + 12*i^2*j^2*b*h - 3*i^2*f^2*b^2 - 
    12*d*g^2*f^2*h - 12*i^2*g*f*e^2 - 24*i^2*j^2*e^2 + 12*i^2*j*f*b*e - 
    12*i*g^3*f*c + 12*i*g^2*j^2*c + 36*i*g^2*j*f*h - 12*i*g^2*f^2*e - 
    36*i*g*j^3*h - 12*i*g*j^2*f*e + 12*i*g*j*f^2*b + 24*i*j^4*e - 
    12*i*j^3*f*b + 8*g^3*f^3 - 8*j^6

  R = KK[x,y,z]
  F = a*x^3 + b*y^3 + c*z^3 + 3*d*x^2*y + 3*e*y^2*z + 3*f*z^2*x + 3*g*x*y^2 + 3*h*y*z^2 + 3*i*z*x^2 + 6*j*x*y*z

  assert(aronholdS F == S)
  elapsedTime assert(aronholdT F == T) -- 21 seconds, but comes out correct.
  
  (S,T) = toSequence aronhold F;
  assert(aronholdS F == S)
  assert(aronholdT F == T) -- aronholdT is really too intensive for impatient people.
///

hesseForm = method()
hesseForm RingElement := F -> (
    R := ring F;
    if numgens R != 3 or first degree F =!= 3 or not isHomogeneous F
    then error "expected homogeneous cubic in a ring with 3 variables";
    -- TODO: this doesn't check that the degrees are all 1 yet...
    
    )

TEST ///
-- from cubicForm Xs#(12,0)
  RQ = QQ[a,b,c]
  F = -a^3+9*a^2*b-9*a*b^2+3*b^3+3*a^2*c-3*a*c^2+c^3
  factor det hessian F
  -- 3 factors:
  a --> a, a-c -> b, a-b -> c
  phi = map(RQ, RQ, {a, a-c, a-b})
  phi^-1 F -- same!
  phi F
///

end--

-- notes on Elsenhans-Jahnel: Computing invariants of cubic surfaces (2020).

-- Aronhold invariants S, T of cubic plane curves.
S = ZZ/32003[x_1..z_4]
DM = transpose genericMatrix(S, 4, 3)

C = x_1^3 + y_1^3 + z_1^3 - 7 * x_1 * y_1 * z_1
F4 = product for i from 1 to 4 list sub(F, {x_1 => x_i, y_1 => y_i, z_1 => z_i});
size  F4
diff((det DM_{3,0,1}), F4);
diff((det DM_{2,3,0}), oo);
diff((det DM_{1,2,3}), oo);
diff((det DM_{0,1,2}), oo);


KK = QQ[a,b,c,d,e,f,g,h,i,j]
S4 = KK[x_1..z_4]
DM = transpose genericMatrix(S4, 4, 3)

F = a*x_1^3 + b*y_1^3 + c*z_1^3 + 3*d*x_1^2*y_1 + 3*e*y_1^2*z_1 + 3*f*z_1^2*x_1 + 3*g*x_1*y_1^2 + 3*h*y_1*z_1^2 + 3*i*z_1*x_1^2 + 6*j*x_1*y_1*z_1

-- Compute S
F4 = product for i from 1 to 4 list sub(F, {x_1 => x_i, y_1 => y_i, z_1 => z_i});
size  F4
diff((det DM_{3,0,1}), F4);
diff((det DM_{2,3,0}), oo);
diff((det DM_{1,2,3}), oo);
diff((det DM_{0,1,2}), oo);
aronholdS = 1/31104 * oo

F = a*x^3 + b*y^3 + c*z^3 + 3*d*x^2*y + 3*e*y^2*z + 3*f*z^2*x + 3*g*x*y^2 + 3*h*y*z^2 + 3*i*z*x^2 + 6*j*x*y*z


aronholdS == a*g*e*c - a*g*h^2 - a*j*b*c + a*j*e*h + a*f*b*h - a*f*e^2 - d^2*e*c + d^2*h^2 + d*i*b*c - d*i*e*h + d*g*j*c - d*g*f*h - 2*d*j^2*h + 3*d*j*f*e - d*f^2*b - i^2*b*h + i^2*e^2 - i*g^2*c + 3*i*g*j*h - i*g*f*e - 2*i*j^2*e + i*j*f*b + g^2*f^2 - 2*g*j^2*f + j^4

restart
debug needsPackage "StringTorics"

KK = QQ[a,b,c,d,e,f,g,h,i,j]
S4 = KK[x,y,z]
F = a*x^3 + b*y^3 + c*z^3 + 3*d*x^2*y + 3*e*y^2*z + 3*f*z^2*x + 3*g*x*y^2 + 3*h*y*z^2 + 3*i*z*x^2 + 6*j*x*y*z
aronholdS F
aronholdS(13*6*x*y*z)

G = elapsedTime aronholdT F;
DM = genericMatrix(ring G, 3, 6)
elapsedTime diff((det DM_{3,4,5}), G);
elapsedTime diff((det DM_{3,4,5}), oo);
elapsedTime diff((det DM_{2,0,5}), oo);
elapsedTime diff((det DM_{1,2,4}), oo);
elapsedTime diff((det DM_{0,1,3}), oo);
elapsedTime diff((det DM_{0,1,2}), oo);


