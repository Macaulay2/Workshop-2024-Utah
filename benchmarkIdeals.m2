-- file with ideals to benchmark radical code


------------------------------------------------------------------------------------------
-- random PID
------------------------------------------------------------------------------------------




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




------------------------------------------------------------------------------------------
-- kernel of a map from a polynomial ring k[x_1..x_n] to a ring of the form
-- k[t^S, x^[p]], where S is a numerical semigroup and x^[p] is the ideal of pth-powers
-- of some variables
------------------------------------------------------------------------------------------




