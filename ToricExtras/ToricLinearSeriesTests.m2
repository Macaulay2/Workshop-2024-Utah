-- test 0: test linear series type
TEST ///
    series = new ToricLinearSeries from hashTable {};
    assert (series === series)
///

-- test 1: test the monomial constructor for ToricLinearSeries
TEST ///
    P2 = toricProjectiveSpace 2;
    S = ring P2;

    m = {x_0^2, x_0*x_1, x_1^2};
    s = toricLinearSeries m;

    assert(monomials(s) == m)
///

-- test 2: test the isComplete function 
TEST ///
    P2 = projectiveToricVariety 2;
    S = ring P2;
    M1 = toricLinearSeries flatten entries basis(1,S);
    M2 = toricLinearSeries {x_0};
    M3 = toricLinearSeries {x_0, x_0, x_0};
    M4 = toricLinearSeries {x_0, x_1, x_2, x_2};
    M5 = toricLinearSeries {x_2, x_1, x_0};
    assert isComplete M1;
    assert not isComplete M2;
    assert not isComplete M3;
    assert not isComplete M4;
    assert isComplete M5;
///