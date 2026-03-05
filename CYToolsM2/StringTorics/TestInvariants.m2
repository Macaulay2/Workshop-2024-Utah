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
