TEST ///
  debug StringTorics
  R = QQ[a..d]
  (A, phi) = genericLinearMap R
  TR = target phi
  assert(source phi === TR)
  assert(ring A === coefficientRing TR)
  for i from 0 to 3 do
    assert(phi TR_i == (A^{i} * (transpose vars TR))_(0,0))

  linearEquationConstraints(A, phi, {}, {
          {{1,0,0,0}, {1,1,0,0}},
          {{0,1,0,0}, {1,1,3,7}},
          {{0,0,1,0}, {5,6,-2,8}},
          {{0,0,0,1}, {0,1,0,0}}}
      )

  F1 = -2*a^3-6*a^2*b+6*b^2*c-12*a^2*d+12*a*b*d+12*b^2*d+36*b*c*d+30*a*d^2+60*b*d^2+54*c*d^2+76*d^3
  F7 = -2*a^3+6*a^2*b-6*a*b^2+2*b^3-6*a*c^2-4*c^3+6*a^2*d-6*a*d^2+2*d^3

  (A0, phi0, I1) = linearEquationConstraints(A, phi, {
          {b+2*d, a-d},
          {b+3*d, a+c},
          {a+b+2*d, a-b},
          {F1, F7}
          }, {
          }
      )
  phi0 F1 == F7

  -- this one isn't correct yet.
  (A0, phi0, I1) = linearEquationConstraints(A, phi, {
          {b+2*d, a-d},
          {b+3*d, a+c},
          {a+b+2*d, a-b}
          }, {
          {{0,0,1,0}, {1,1,-1,1}}
          }
      )
--

  A0 % sub(trim ideal last coefficients(phi0 F1 - F7), coefficientRing TR)
///
