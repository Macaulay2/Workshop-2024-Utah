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
    P2 = toricProjectiveSpace 2;
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

-- test 3: test the baseLocusIdeal function
TEST ///
    P2 = toricProjectiveSpace 2;
    S = ring P2;
    M1 = toricLinearSeries flatten entries basis(1,S);
    M2 = toricLinearSeries {x_0};
    M3 = toricLinearSeries {x_0, x_0, x_0};
    M4 = toricLinearSeries {x_0, x_1, x_2, x_2};
    M5 = toricLinearSeries {x_2, x_1, x_0};
    assert (baseLocusIdeal M1 == ideal(x_0, x_1, x_2));
    assert (baseLocusIdeal M2 == ideal(x_0));
    assert (baseLocusIdeal M3 == ideal(x_0));
    assert (baseLocusIdeal M4 == ideal(x_0, x_1, x_2));
    assert (baseLocusIdeal M5 == ideal(x_0, x_1, x_2));
///

-- test 4: test the isBasePointFree function
TEST ///
    P2 = toricProjectiveSpace 2;
    S = ring P2;
    M1 = toricLinearSeries flatten entries basis(1,S);
    M2 = toricLinearSeries {x_0};
    M3 = toricLinearSeries {x_0, x_0, x_0};
    M4 = toricLinearSeries {x_0, x_1, x_2, x_2};
    M5 = toricLinearSeries {x_2, x_1, x_0};
    assert (isBasepointFree M1)
    assert (not isBasepointFree M2)
    assert (not isBasepointFree M3)
    assert (isBasepointFree M4)
    assert (isBasepointFree M5)


    P1 = toricProjectiveSpace 1;
    X = P1 ** P1;
    T = ring X;
    M1 = toricLinearSeries {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2},T_{2,0,2,0}};
    M2 = toricLinearSeries {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2}};
    M3 = toricLinearSeries {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2},T_{1,1,2,0}};
    M4 = toricLinearSeries {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2},T_{1,1,1,1}};
    assert(isBasepointFree M1)
    assert(not isBasepointFree M2)
    assert(not isBasepointFree M3)
    assert(not isBasepointFree M4)
///