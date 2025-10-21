TEST ///
  R = ZZ/101[x_0] / ideal(x_0^3)
  M = coker vars R
  assert try ( singularityToDerived(1, M, LengthLimit => 1); false ) else true
///

TEST ///
///
