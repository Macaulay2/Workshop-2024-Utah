ToricLinearSeries = new Type of HashTable
ToricLinearSeries.synonym = "toric linear series"

toricLinearSeries = method(TypicalValue => ToricLinearSeries)
toricLinearSeries List := ToricLinearSeries => m -> (
    if #m == 0 then error "toricLinearSeries expects a nonempty list";
    d := degree m#0;
    if not all(m, x -> d == degree x) then error "toricLinearSeries expects a list of monomials of the same degree";
    s := new ToricLinearSeries from {
        "monomials" => m,
        "degree" => d
    };
    s
)

monomials ToricLinearSeries := List => o -> s -> (
    s#"monomials"
)

degree ToricLinearSeries := s-> (
    s#"degree"
)
-- the monomial map that defines a toric map
-- monomials ToricMap := List => o -> f -> first entries matrix inducedMap f

--isComplete = method(TypicalValue => Boolean)

isComplete ToricLinearSeries := linSeries -> (
    m := monomials linSeries;
    setM := set m;
    d := degree linSeries;
    if #m != #setM then return false;
    degDMonomials := flatten entries basis(d, ring m#0);
    setM == set degDMonomials
)

baseLocusIdeal = method(TypicalValue => Ideal)

baseLocusIdeal ToricLinearSeries := linSeries -> (
    ideal monomials linSeries
)

-- -- helper for listing monomials of given degree in the ring
-- -- TODO: move to Core
-- monomials(ZZ,   Ring) :=
-- monomials(List, Ring) := List => o -> (d, S) -> first entries basis(d, S)


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
-- isBasepointFree = (X, L) -> set(intersect decompose ideal L)_* == set(ideal X)_*

-- -- lists all toric linear series of Proj S in degree d, including the complete one
-- allToricLinearSeries = (d, S) -> select(subsets monomials(d, S), mons -> isBasepointFree(variety S, mons))

-- TEST ///
--   -- embedding of P1xP1 in P8 via |O(2,2)|
--   D = smallAmpleToricDivisor(2, 0)
--   X = variety D -- P1xP1
--   S = ring X
--   phi = map D
--   assert isWellDefined phi
--   assert(first entries matrix inducedMap phi === monomials D)

--   -- embedding from ample divisor
--   X = (toricProjectiveSpace 1)^**2
--   D = toricDivisor({1,0,1,0}, X)
--   phi = map D
--   assert isWellDefined phi
--   assert(first entries matrix inducedMap phi === monomials D)

--   -- embedding from a complete linear series
--   X = (toricProjectiveSpace 1)^**2
--   -- FIXME: make this work with a ToricLinearSeries object
--   L = monomials ({1,1}, ring X)
--   phi = map L
--   assert isWellDefined phi
--   assert(first entries matrix inducedMap phi === L)

--   -- rational map from an incomplete linear series
--   X = toricProjectiveSpace 1
--   Y = toricProjectiveSpace 2
--   R = ring Y
--   S = ring X
--   -- cuspidal cubic, not projectively normal
--   phi = map(Y, X, {S_{1,2}, S_{3,0}, S_{0,3}})
--   assert isWellDefined phi
--   assert(first entries matrix inducedMap phi === monomials phi)
--   -- smooth conic, projectively normal
--   phi = map(Y, X, {S_{2,0}, S_{1,1}, S_{0,2}})
--   assert isWellDefined phi
--   assert(first entries matrix inducedMap phi === monomials phi)

--   -- the twistic cubic as the embedding of P1 in P3
--   Y = toricProjectiveSpace 3
--   R = ring Y
--   phi = map(Y, X, {S_{3,0}, S_{2,1}, S_{1,2}, S_{0,3}})
--   assert isWellDefined phi
--   assert(first entries matrix inducedMap phi === monomials phi)

--   -- rational map from an incomplete series
--   X = toricProjectiveSpace 1
--   Y = toricProjectiveSpace 3
--   phi = map(Y, X, matrix vector {1,3,5})
--   assert isWellDefined phi
--   assert(first entries matrix inducedMap phi === monomials phi)
-- ///
