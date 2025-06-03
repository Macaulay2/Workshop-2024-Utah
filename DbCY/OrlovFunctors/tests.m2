TEST ///
  R = ZZ/101[x_0] / ideal(x_0^3)
  M = coker vars R
  assert try ( singularityToDerived(M, 1, 1); false ) else true
///

TEST ///
///
