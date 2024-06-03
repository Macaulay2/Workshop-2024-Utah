TEST ///
    X = Proj(QQ[x,y])
    F = (ring X)^2
    projectiveBundle(X,sheaf F)
    --TODO: add checks for variety, sheaf, rank of this bundle
    T = QQ[x]
    E = T^1
    assert instance( try projectiveBundle (X,sheaf E), Nothing)
///

TEST ///

///
