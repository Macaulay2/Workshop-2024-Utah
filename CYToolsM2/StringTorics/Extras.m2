------------------------------------------------
-- TODO: what to do with the code below this? -- for now, placed it here...
-- i.e. it is older code, but potentially useful
------------------------------------------------

protect MYRING
---------------------------------------------------
-- Simplicial complex-like code
subcomplex = method()
subcomplex(NormalToricVariety, List) := (V, D) -> (
    if not V.?MYRING then (
        x := getSymbol "x";
        V.MYRING = (ZZ/32003) (monoid([x_0..x_(#rays V-1)]));
        );
    S := V.MYRING;
    faces := unique for f in max V list sort toList (set D * set f);
    -- Next we take these faces, and make monomials out of them.
    monoms := flatten for f in faces list (1_S * product (f/(i -> S_i)));
    simplicialComplex monoms
    )
subcomplex(NormalToricVariety, List) := (V, D) -> (
    if not V.?MYRING then (
        x := getSymbol "x";
        V.MYRING = (ZZ/32003) (monoid([x_0..x_(#rays V-1)]));
        );
    S := V.MYRING;
    R := (coefficientRing S)(monoid[for d in D list S_d]);
    trD := hashTable for i from 0 to #D-1 list D#i => i;
    faces := unique for f in max V list sort toList (set D * set f);
    facesTr := for f in faces list (f/(i -> trD#i));
    -- Next we take these faces, and make monomials out of them.
    monoms := flatten for f in facesTr list (1_R * product (f/(i -> R_i)));
    simplicialComplex monoms
    )
skeleton(ZZ, SimplicialComplex) := (d,C) -> (
    S := ring C;
    fac := flatten entries facets C;
    fac1 := for f in fac list (support f)/index//sort;
    print fac1;
    subfacs := unique flatten for f in fac1 list subsets(f, d+1);
    simplicialComplex for f in subfacs list (1_S * product(f/(i -> S_i)))
    )
simplicialComplex NormalToricVariety := (V) -> (
    if not isSimplicial V then error "expected a simplicial polytope";
    S := ring V;
    simplicialComplex for f in max V list (1_S * product(f/(i -> S_i)))
    )
subcomplex(SimplicialComplex, List) := (C,D) -> (
    S := ring C;
    if not all(D, v -> instance(v,S)) then
      D = D/(i -> S_i);
    elems := toList (set gens S - D);
    if #elems == 0 then return C;
    J := monomialIdeal elems;
    simplicialComplex(J + monomialIdeal C)
    )

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
---------------------------------------------------
hodgeOfCYToricDivisor = method();

-------------------------------------------
-- Newer code -----------------------------
-------------------------------------------
hodgeCY = (P, pt) -> (
    -- pt should be a lattice point: list of integers, a point in P2 = polar P.
    --if dim P != 4 then error "expected a 4d reflexive polytope";
    pt = for a in pt list if instance(a,QQ) then lift(a,ZZ) else a;
    if not all(pt, i -> instance(i,ZZ))
    then error "expected a list of integers";
    P2 := polar P;
    L := latticePointHash P2;
    if not L#?pt then error "point provided is not a lattice point of the dual polytope";
    f := minimalFace(P2, pt);
    i := dim(P2,f);
    dualf := dualFace(P2, f);
    g := # interiorLatticePointList(P, dualf);
    result := new MutableList from splice{1,(dim P-2):0};
--    if i == 0 then {1,0,g}
--    else if i == 1 then {1,g,0}
--    else if i == 2 then {1+g,0,0};
    if i >= 0 and i <= dim P-2 then result#(dim P - 2 - i) = result#(dim P - 2 - i) + g
    else return null;
    toList result
    )

hodgeOfCYToricDivisors = method()
hodgeOfCYToricDivisors Polyhedron := (P) -> (
    P2 := polar P;
    L := latticePointList P2;
    -- only go to #L-2, since the origin is L#(#L-1)
    hashTable for i from 0 to #L-2 list i => hodgeCY(P, L#i)
    )

hodgeOfCYToricDivisor(Polyhedron,List) := (P,pt) -> hodgeCY(P, pt)

h11OfCY Polyhedron := (P1) -> (
    elapsedTime np := 1 + #(latticePointList polar P1) - 1;  -- -1 for the origin
    elapsedTime A0 := annotatedFaces(0,P1);
    elapsedTime A1 := annotatedFaces(1,P1);
    t := A0/last//sum;
    t1 := A1/(v -> v#2 * v#3)//sum;
    np - dim P1 - 1 - t + t1
    )

h21OfCY Polyhedron := (P1) -> (
    np := 1 + #(latticePointList P1) - 1; -- -1 for the origin
    A2 := annotatedFaces(2,P1);
    A3 := annotatedFaces(3,P1);
    t := A3/(v -> v#2)//sum;
    t1 := A2/(v -> v#2 * v#3)//sum;
    np - 5 - t + t1
    )

isFavorable Polyhedron := (P1) -> (
    -- This is from the Batyrev formula for h^11, 
    -- The term ell^*(theta) * ell^*(theta^*) gives new divisors.
    --  here theta is an edge of P, theta^* is a (dim P)-2 face of (polar P)
    --  and ell^* is the number of interior lattice points.
    A1 := annotatedFaces(1,P1);
    t1 := A1/(v -> v#2 * v#3)//sum;
    t1 == 0
    )

complexWithInteriorEdgesRemoved = method()
complexWithInteriorEdgesRemoved(SimplicialComplex,Polyhedron,ZZ) := (C,P2,facedim) -> (
    LPs := latticePointList P2;
    edgeList := faces(1,C);
    badEdges := select(edgeList, e -> (
            f := (support e)/index; 
            fpts := f/(i -> LPs#i);
            dim(P2,minimalFace(P2,fpts)) >= facedim)
        );
    if #badEdges == 0 then C else simplicialComplex(monomialIdeal C + monomialIdeal badEdges)
    )

complexWithInteriorFacesRemoved = method()
complexWithInteriorFacesRemoved(SimplicialComplex,Polyhedron,ZZ) := (C,P2,facedim) -> (
    LPs := latticePointList P2;
    facelist := flatten for i from 0 to dim P2-1 list faces(i,C);
    badfaces := select(facelist, e -> (
            f := (support e)/index;
            fpts := f/(i -> LPs#i);
            dim(P2,minimalFace(P2,fpts)) >= facedim)
        );
    if #badfaces == 0 then C else simplicialComplex(monomialIdeal C + monomialIdeal badfaces)
    )

tentativeHodgeVector = method()
    
-- V: a normal toric variety corresponding to a fine star
-- triangulation of P2 using all of the non-zero lattice points of P2,
-- except possibly those that are interior to 3d faces.
-- P2: reflexive polytope of dimension 4.
-- D is a list of indices into 'latticePointList P2'
-- OUTPUT:
-- the vector {h^0(D, OO_D), h^1(D, OO_D), h^2(D, OO_D)}, where
-- D is considered to be the sum of the effective toric divisors 
-- in the input D.

hodgeVectorViaTheorem = method()
hodgeVectorViaTheorem(NormalToricVariety, Polyhedron, List) := (V,P2,D) -> (
    -- we assume that the list of rays of V is identical to the first 
    -- elements of latticePointsList P2:
    if debugLevel >= 1 then << "starting tentativeHodgeVector" << endl;
    if rays V != take(latticePointList P2, # rays V) then 
      error "logic error: we need to translate from rays V to lattice points of P2";
    (h0,h1,h2) := (0,0,0);
    Delta' := subcomplex(V,D);
    Delta := complexWithInteriorFacesRemoved(Delta',P2,3);
    deltaIndices := hashTable for i from 0 to #D-1 list D#i => i; -- need to use these indices when dealing with Delta',Delta.
    h0 = 1 + rank HH_0(Delta);
    h1 = rank HH_1(Delta);
    h2 = rank HH_2(Delta);
    if debugLevel >= 1 then
      << "orig h: " << {h0,h1,h2} << endl;
    -- Contribution #1: loop over all vertices of P2
    --   for each element v of the set D also a verte    x of P2:
    --     add in genus(P2,{v}) to h2
    vertsP2 := (faceList(0,P2))/first; -- indices of vertices in P2
    vertsKD := sort toList (set vertsP2 * set D);
    if debugLevel >= 1 then << "about to do vertex contributions" << endl;
    for v in vertsKD do (
          g := genus(P2, {v});
          if debugLevel >= 1 and g > 0 then
            << "vertex: " << v << " adding " << g << " to h2" << endl;
          h2 = h2 + g;
        );
    -- Contribution #2: loop over all edges of P2 that meet D in some way
      -- note: {edge, lattice pts on edge, genus}
    edgeGenera := for e in faceList(1,P2) list {e, latticePointList(P2,e), genus(P2,e)};
    if debugLevel >= 2 then print netList edgeGenera;
    for e in edgeGenera do (
          g := e#2;
          if g == 0 then continue;
          if isSubset(e#1, D) then (
              if debugLevel >= 1 and g > 0 then
                << "edge: " << e << " adding " << g << " to h2" << endl;
              h2 = h2 + g
              )
          else (
              nelems := #((set e#1) * (set D));
              nvertices := #((set e#0) * (set D));
              if nelems > 0 then (
                  eindices := for i in e#1 list if deltaIndices#?i then deltaIndices#i else continue;
                  Delta1 := subcomplex(Delta,eindices);
                  if debugLevel >= 1 then
                    << "edge: " << e#1 << " nvertices: " << nvertices << " Delta1: " << Delta1 << endl;
                  amt := (1 + rank HH_0(Delta1) - nvertices) * g;
                  if debugLevel >= 1 then 
                    << "    : " << e << " adding " << amt << " to h1" << endl;
                  h1 = h1 + amt
                  );
              );
          );
      faceGenera := for t in faceList(2,P2) list {t, latticePointList(P2,t), interiorLatticePointList(P2,t), genus(P2,t)};
      if debugLevel >= 2 then print netList faceGenera;
      if debugLevel >= 1 then << "in facGenera" << endl;
      for t in faceGenera do (
          g := t#3;
          if g > 0 then (
              nverts := #((set t#0) * (set D));
              nlattice := #((set t#1) * (set D));
              ninterior := #((set t#2) * (set D));
              -- 3 cases: 
              if nlattice == 0 then (
                  -- do nothing
                  )
              else if nlattice > 0 and nlattice == ninterior then (
                  if debugLevel >= 1 and g > 0 then 
                    << "2face: " << t#0 << " adding " << g << " to h0" << endl;
                  h0 = h0 + g;
                  )
              else if nlattice == #t#1 then ( -- the whole face t is in D
                  if debugLevel >= 1 and g > 0 then                   
                    << "2face: " << t#0 << " adding " << g << " to h2" << endl;
                  h2 = h2 + g;
                  )
              else (
                  tindices := for i in t#1 list if deltaIndices#?i then deltaIndices#i else continue;
                  KDf := subcomplex(Delta,tindices); -- K_D \cap f
                  df := sort toList(set t#1 - set t#2); -- boundary(f), as a list of lattice points
                  comps := decompose KDf; -- each connected component y of K_D \cap f
                  comps1 := for c in comps list (
                      dfindices := for i in df list if deltaIndices#?i then deltaIndices#i else continue;
                      c1 := subcomplex(c,dfindices);
                      -- now remove all edges that go in the interior of the 2-face
                      -- y \cap df, for each y
                      c2 := complexWithInteriorEdgesRemoved(c1,P2,2);
                      c2);
                  n1 := sum for c in comps1 list rank HH_0(c);
                  if debugLevel >= 1 and g > 0 and n1 > 0 then                   
                    << "2face: " << t#0 << " adding " << n1*g << " to h1" << endl;                  
                  h1 = h1 + n1*g;
                  )
              );
          );
    {h0,h1,h2}
    )

tentativeHodgeVector(NormalToricVariety, Polyhedron, List) := (V, P2, D) -> hodgeVectorViaTheorem(V,P2,D)
