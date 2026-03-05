-- Scratch code from CalabiYauInToric.m2
-- Testing restrictTriangulation, automorphisms, triangulation equivalence

///
  -- Let's test restrictTriangulation, automorphisms and getting only triangulations
  -- unique up to linear automorphism of the lattice polytope.
  topes = kreuzerSkarke(3, Limit => 20);
  topes_15
  tope = KSEntry "4 7  M:74 7 N:8 7 H:3,63 [-120] id:15
   1   0   0   0  -3  -3   3
   0   1   1   1   1  -2  -2
   0   0   3   0   3   0  -6
   0   0   0   3  -3  -6   6
  "
  Q = cyPolytope tope
  Xs = findAllCYs Q
  #Xs
  auts = for e in automorphisms Q list matrix e
  rayvecs = for v in rays Q list transpose matrix {v}
  rayHash = hashTable for i from 0 to #rayvecs - 1 list rayvecs#i => i
  g = auts#0
  perms = for g in auts list
    for v in rayvecs list rayHash#(g*v)
  sort unique flatten restrictTriangulation(3, Xs#0)
  sort for x in oo list sort perms_0_x
///
