doc ///
    Key
        toricLinearSeries
        (toricLinearSeries, List)
    Headline
        constructor for toric linear series
    Usage
        toricLinearSeries L
    Inputs
        L : List
            of monomials of the Cox ring    
    Description
    	Text
            Constructs a toric linear series from a list of monomials from the Cox ring of a toric variety.
            It checks that all monomials are of the same degree.
            It expects the monomials to be in ring(X) where X is a normal toric variety.
        Example
            P2 = toricProjectiveSpace 2;
            coxRing = ring P2;
            mons = {x_0^2, x_0*x_1, x_1^2};
            toricLinearSeries mons
///

doc ///
    Key
        ToricLinearSeries
    Headline
        Linear series on a toric variety
    Description
    	Text
	        Work in progress implementation of linear series on a toric variety
    SeeAlso 
        ToricDivisor
///

doc ///
    Key
        (monomials,ToricLinearSeries)
    Headline
        method to extract the list of monomials encoding a (generating set for) linear series
    Usage
        monomials s
    Inputs
        s : ToricLinearSeries
    Description
        Text
            Returns the list of monomials encoding a linear series
        Example
            P2 = toricProjectiveSpace 2;
            coxRing = ring P2;
            mons = {x_0^2, x_0*x_1, x_1^2}
            monomials toricLinearSeries mons
///
