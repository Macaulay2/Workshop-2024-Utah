ToricLinearSeries = new Type of HashTable
ToricLinearSeries.synonym = "toric linear series"

toricLinearSeries = method(TypicalValue => ToricLinearSeries)
toricLinearSeries List := ToricLinearSeries => m -> (
    s := new ToricLinearSeries from {
        "monomials" => m  
    };
    s
)

monomials ToricLinearSeries := List => o -> s -> (
    s#"monomials"
)
-- the monomial map that defines a toric map
monomials ToricMap := List => o -> f -> first entries matrix inducedMap f
-- helper for listing monomials of given degree in the ring
-- TODO: move to Core
monomials(ZZ,   Ring) :=
monomials(List, Ring) := List => o -> (d, S) -> first entries basis(d, S)


-- getting map of tori from a divisor or linear series
map(NormalToricVariety, NormalToricVariety, ToricLinearSeries) :=
map(NormalToricVariety, NormalToricVariety, ToricDivisor) := ToricMap => opts -> (Y, X, D) -> map(Y, X, monomials D)
map(NormalToricVariety, NormalToricVariety, List)         := ToricMap => opts -> (Y, X, L) -> map(Y, X, transpose(
	divmap := matrix transpose(first \ exponents \ L); -- map CDiv Y -> CDiv X
	(divmap * matrix rays Y) // matrix rays X))        -- map    M_Y -> M_X

-- if the source is not given, get the variety from the linear series or divisor
map(NormalToricVariety, Nothing, ToricLinearSeries) :=
map(NormalToricVariety, Nothing, ToricDivisor) := ToricMap => opts -> (Y, X, D) -> map(Y, variety D, monomials D)
map(NormalToricVariety, Nothing, List)         := ToricMap => opts -> (Y, X, L) -> map(Y, variety ring L#0, L)

-- if target is not given, use the appropriate projective space
map(Nothing, NormalToricVariety, ToricLinearSeries) :=
map(Nothing, NormalToricVariety, ToricDivisor) := ToricMap => opts -> (Y, X, D) -> map(, X, monomials D)
map(Nothing, NormalToricVariety, List)         := ToricMap => opts -> (Y, X, L) -> map(toricProjectiveSpace(#L - 1), X, L)

-- if neither source or target is given, deduce both of them!
map ToricLinearSeries :=
map ToricDivisor      := ToricMap => opts -> D -> map(, variety D, monomials D)
map List              := ToricMap => opts -> L -> map(, variety ring L#0, L)

-- whether a linear series over a Cox ring is basepoint free
-- TODO: make Ideal == Ideal work when one is ideal () with ring ZZ
isBasepointFree = (X, L) -> set(intersect decompose ideal L)_* == set(ideal X)_*

-- lists all toric linear series of Proj S in degree d, including the complete one
allToricLinearSeries = (d, S) -> select(subsets monomials(d, S), mons -> isBasepointFree(variety S, mons))
