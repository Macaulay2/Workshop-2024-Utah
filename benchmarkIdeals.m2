-- file with ideals to benchmark radical code


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
-- kernel of a map from a polynomial ring k[x_1..x_n] to a ring of the form
-- k[t^S, x^[p]], where S is a numerical semigroup and x^[p] is the ideal of pth-powers
-- of some variables
------------------------------------------------------------------------------------------




