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

-- test 5: test the toricMap function
TEST ///
    P1 = toricProjectiveSpace 1;
    X = P1 ** P1;
    T = ring X;
    M1 = toricLinearSeries {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2},T_{2,0,2,0}};
    M2 = toricLinearSeries {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2}};
    M3 = toricLinearSeries {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2},T_{1,1,2,0}};
    --M4 = toricLinearSeries {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2},T_{1,1,1,1}};
    M4 = toricLinearSeries {T_{2,0,0,1}, T_{2,0,1,0}, T_{0,2,1,0}, T_{0, 2, 0, 1}};

    map1 = toricMap(M1);
    map2 = toricMap(M2);
    map3 = toricMap(M3);
    map4 = toricMap(M4);
    
    assert(isWellDefined(map1));
    assert(not isWellDefined(map2));
    assert(not isWellDefined(map3));
    assert(isWellDefined(map4));

    assert(source(map1) === X);
    assert( (#(rays target(map1)) == 4) and isSmooth(target(map1)) and (dim(target(map1)) == 3)); --checks if target is P3

    assert(source(map2) === X);
    assert( (#(rays target(map2)) == 3) and isSmooth(target(map2)) and (dim(target(map2)) == 2)); --checks if target is P2
    
    assert(source(map3) === X);
    assert( (#(rays target(map3)) == 4) and isSmooth(target(map3)) and (dim(target(map3)) == 3)); --checks if target is P3
    
    assert(source(map4) === X);
    assert( (#(rays target(map4)) == 4) and isSmooth(target(map4)) and (dim(target(map4)) == 3)); --checks if target is P3

///

-- test 6: test ideal function gives ideal of image of the induced toric map by a toric linear series
TEST ///
    P1 = toricProjectiveSpace 1;
    X = P1 ** P1;
    T = ring X;
    M1 = toricLinearSeries {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2},T_{2,0,2,0}};
    I = ideal toricMap M1
    S = ring I
    assert(idealOfImage(M1, TargetRing => S) == I)

    S = ring P1;
    v3 = toricLinearSeries first entries basis(3,S);
    
    needsPackage "Resultants";
    idealOfVeronese3 = kernel veronese(1, 3, QQ);
    assert(idealOfImage(v3, TargetRing => ring idealOfVeronese3) == idealOfVeronese3);
///