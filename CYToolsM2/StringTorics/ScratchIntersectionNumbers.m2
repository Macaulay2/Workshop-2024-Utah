-- Scratch code from IntersectionNumbers.m2

-- Old version of toBasisIntersectionNumbers (superseded by version taking nonfavsHash)
-*
toBasisIntersectionNumbers = (toricIntersectionNumbers, basIndices) -> (
    H := hashTable for i from 0 to #basIndices-1 list basIndices#i => i;
    for t in toricIntersectionNumbers list (
        if isSubset(t#0, basIndices) then t#0/(a -> H#a)//sort => t#1 else continue
        )
    )
*-

------------------------------
-- REMOVE: tripleProductsCY --
------------------------------
-- This code is no longer simpler than current code.
-- Simpler code, used to debug the algorithm/implementation above.
tripleProductsCY = method()
tripleProductsCY NormalToricVariety := (V) -> (
    elapsedTime AV := abstractVariety(V, point);
    IV := intersectionRing AV;
    h := sum gens IV; -- Calabi-Yau hyperplane class in V.
    J := ideal select((ideal IV)_*, f -> size f == 1);
    forceGB gens J;
    gens gb J;
    A := (ring J)/J;
    monoms := ideal basis(3, A);
    elapsedTime JV := sub(monoms,IV);
    elapsedTime (JVh := h ** (gens JV));
    flatJVh := flatten entries JVh;
    hashTable for i from 0 to numgens monoms - 1 list (
        m := monoms_i;
        d := integral(flatJVh#i);
        if d > 0 then m => d else continue
        )
    )

------------------------------
-- REMOVE: possibleNonZeros --
------------------------------
possibleNonZeros = (V) -> (
    -- assumption currently: V has dim 4, is reflexive, and X is the anti-canonical CY3 divisor.
    -- returns a list of lists of 3 integers (0 <= i1 <= i2 <= i3 <= N-1)
    --  where N = #rays V.
    -- and all triples other than those on this list must have triple intersection
    -- on X being zero.
    P2 := convexHull transpose matrix rays V;
    F := annotatedFaces P2;
    faces2 := select(F, f -> f#0 == 2);
    faces2 = faces2/(x -> x#2); -- this is a list of all 2-faces in the polytope,
    -- with which rays are on each face.
    -- any triple not supported on a 2-face will have triple intersection zero.
    triangles := (max V)/(t -> subsets(t,3))//flatten//unique;
    edges := (max V)/(t -> subsets(t,2))//flatten//unique//sort;
    triples := sort flatten for f in faces2 list select(triangles, t -> isSubset(t,f));
    singles := for i from 0 to # rays V - 1 list {i,i,i};
    doubles := sort flatten for f in faces2 list select(edges, t -> isSubset(t,f));
    doubles = unique flatten for x in doubles list {{x#0,x#0,x#1},{x#0,x#1,x#1}};
    {singles,doubles,triples}
    )

--------------------------------------
-- REMOVE: CY3NonzeroMultiplicities --
--------------------------------------
  CY3NonzeroMultiplicities = method()
  CY3NonzeroMultiplicities NormalToricVariety := (V) -> (
      RAYS := transpose matrix rays V;
      P2 := convexHull RAYS;
      (singles,doubles,triples) := toSequence possibleNonZeros V;
      doubles = doubles/unique/sort//unique;
      singles = singles/unique/sort//unique;
      mult3 := new MutableHashTable from for x in triples list (
           x => 1 + genus(P2, minimalFace(P2, (rays V)_x))
           );
      multvec := ij -> (
          for ell from 0 to #rays V-1 list (
              if member(ell,ij) then 0
              else (
                  s := sort append(ij,ell);
                  if mult3#?s then mult3#s else 0
                  ))
          );
      for d in doubles do (
          RHS := - RAYS *  transpose (matrix{multvec d});
          vals := flatten entries solve(RAYS_d, RHS);
          d1 := prepend(d#0,d);
          d2 := append(d,d#1);
          if vals#0 != 0 then mult3#d1 = vals#0;
          if vals#1 != 0 then mult3#d2 = vals#1;
          );
      for d in singles do (
          s := {d#0,d#0};
          RHS := - RAYS *  transpose (matrix{multvec s});
          vals := flatten entries solve(RAYS_d, RHS);
          if vals#0 != 0 then mult3#{d#0,d#0,d#0} = vals#0;
          );
      new HashTable from mult3
      )

------------------------------
-- REMOVE: CY3Intersections --
------------------------------
CY3Intersections = method()
CY3Intersections(NormalToricVariety, List) := (V, indexOfDs) -> (
    H := CY3NonzeroMultiplicities V;
    loc := new HashTable from for i from 0 to #indexOfDs-1 list indexOfDs#i => i;
    tr := k -> sort for k1 in k list loc#k1;
    new Array from for kv in pairs H list (
      if not isSubset(kv#0, indexOfDs) then continue;
      new Array from append(tr kv#0, kv#1)
      )
  )
