-- file with ideals to benchmark radical code

-----------------------------------------------------------------------------------------
Ideals associated to determinantal varieties
-----------------------------------------------------------------------------------------
I=ideal(x_2*y_1-x_1*y_2,x_2*y_2-x_1*y_3,y_2^2-y_1*y_3)
J=I^4
------computational time of radical of J is pretty fast.

------------------------------------------------------------------------------------------
-- random principal ideal
------------------------------------------------------------------------------------------

restart

randomMonomial = indets -> (
    upperExponentBound := 10;
    
    numVars := #indets;
    
    possibleExponents := random splice join apply(toList(0..upperExponentBound), n -> (numVars-2*n):n);
    (random ZZ) * (product apply(toList(0..numVars-1), i -> (indets_i)^(possibleExponents_i) ))
    )

randomPolynomial = indets -> (
    numTerms := 5;
    sum apply(toList(1..numTerms), j -> (random ZZ) * randomMonomial(indets))
    )
 
R = QQ[x_1..x_10]
f = product for i from 1 to 3 list ((randomPolynomial(gens R))^(first random toList (1..5)));
I = ideal f;
time radical(I, Strategy=>CompleteIntersection);


-- approx. avg. time for Strategy=>Unmixed: not finishing, quit after a minute or so
-- approx. avg. time for Strategy=>Decompose: about 0.4s, usually between 0.2s and 0.7s
-- approx. avg. time for Strategy=>CompleteIntersection: at most 0.001s


restart


------------------------------------------------------------------------------------------
-- power of a determinantal ideal
------------------------------------------------------------------------------------------

R=QQ[a..z];
ROW=3;
COLUMN=3;
POWER = 5;
I=minors(2,genericMatrix(R,a,ROW, COLUMN));
I=I^POWER;
time radical I

----Decompose is faster than Unmixed here

------------------------------------------------------------------------------------------
-- permanental ideals
------------------------------------------------------------------------------------------

R=QQ[a..z]
loadPackage "Permanents";
ROW=3;
COLUMN=3;
POWER=4;
I=pminors(2,genericMatrix(R,a,ROW, COLUMN));
I=I^POWER;
time radical I

---Decompose is faster than Unmixed


------------------------------------------------------------------------------------------
-- Skew-symmetric matrix Schubert varieties are exactly matrix Schubert varieties
-- with the underlying generic matrix replaced by a generic skew-symmetric matrix.
-- The ideal constructed in this way is in general not radical.
-- Reference: Eric Marberg and Brendan Pawlowski. Groebner geometry for skew-symmetric
--	       matrix Schubert varieties.
------------------------------------------------------------------------------------------

-- two remarks:
-- 1) skew-symmetric matrix Schuberts only make sense for permutations that are
--     fixed-point-free involutions, defined in the method isFPFInv
-- 2) there are no fixed-point-free involutions in S_n for odd n, so we restrict to even n

restart

needsPackage "MatrixSchubert"

isFPFInv = w -> (
    -- returns true iff the permutation w (a list) is a Fixed-Point-Free Involution
    -- def: w is fixed-point-free iff w(i) never equals i
    -- def: w is an involution iff w^2 = id
    -- so w is a fixed-point-free involution if its permutation matrix is symmetric with 0's on the diagonal
    n := #w;
    wMatrix := map(ZZ^n, n, (i,j) -> if w#j == i+1 then 1 else 0);

    (wMatrix == transpose wMatrix) and not any(w, i -> i == w#(i-1))
    )


skewSymmEssentialSet = w -> (
    -- the usual essential set, but only the (i, j) for which i < j
    D := select(rotheDiagram w, L -> (L_0 < L_1));  -- the skew-symmetric diagram
    select(D, L -> ( not( isMember((L_0+1, L_1), D) or  isMember((L_0, L_1+1), D) )))
    )


skewSymmDetIdeal = (R, w) -> (
    -- returns the skew-symmetric matrix Schubert determinantal ideal
    -- the radical of this ideal defines the skew-symmetric matrix Schubert variety
   
    n := #w;
    M := genericSkewMatrix(R, n);

    -- get the essential set and re-index 
    essSet := apply(skewSymmEssentialSet w, L -> (L_0 - 1, L_1 - 1));

    -- find the determinantal ideal
    ranks := rankTable w;
    ideal apply(essSet, L -> (
	    row := L_0;
	    col := L_1;
	    r := ranks_(row, col);

	    minors(r+1, M^(toList(0..row))_(toList(0..col)))
	    )
        )
    )


R = QQ[x_1..x_100];
n = 8;  -- n must be even
Sn = permutations toList (1..n);
fpfInvs = select(Sn, isFPFInv);

for w in fpfInvs do (
    I := skewSymmDetIdeal(R, w);
    time radical(I, Strategy=>Decompose);
    )


-- avg. approx. time for Strategy=>Decompose: most are a fraction of a second (between 0.1s and 0.5s), some are quicker, a couple are 20-30s
-- assumptions not met for CompleteIntersection method



------------------------------------------------------------------------------------------
-- KERNEL of a map from a polynomial ring k[x_1..x_n] to a ring of the form
-- k[t^S, x^[p]], where S is a numerical semigroup and x^[p] is the ideal of pth-powers
-- of some variables
------------------------------------------------------------------------------------------




