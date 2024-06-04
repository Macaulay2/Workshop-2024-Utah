newPackage(
    "RadicalAlgorithms",
    AuxiliaryFiles => false,
    Version => "0.1", 
    Date => "",
    Keywords => {},
    Authors => {
        {Name => "Ayah Almousa", 
            Email => "aalmousa@mailbox.sc.edu", 
            HomePage => "http://sites.google.com/view/ayah-almousa"},
	{Name => "Trung Chau", 
            Email => "chauchitrung1996@gmail.com"},
	{Name => "Mike Cummings",
	    Email => "cummim5@mcmaster.ca",
	    HomePage => "https://math.mcmaster.ca/~cummim5/"},
	{Name => "David Eisenbud",
	    Email => "de@berkeley.edu"},
	{Name => "Justin Fong",
	    Email => "jafong1@gmail.com"},
	{Name => "Manohar Kumar",
	     Email => "manhar349@gmail.com",
             HomePage => "https://sites.google.com/view/manohar-kumar-pmrf-update/"},
	{Name => "Adam LaClair", 
            Email => "alaclair@purdue.edu", 
            HomePage => "https://sites.google.com/view/adamlaclair/home"}
    },
    Headline => "Implementing and benchmarking radical algorithms by Craig Huneke",
    PackageExports => {
    },
    DebuggingMode => true
)

export{
    -- UNTESTED EXPORTS
    "hunekeAlgorithm"
    
    -- TESTED EXPORTS

}


------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **CODE** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------

-------------------------------------------
-- ** Part 1 : Benchmarks **
-------------------------------------------


-----------------------------------
-- ** Part 2: Huneke's Algorithm **
-----------------------------------

genericSocle := I -> (
    d:=dim I;
    return (I:minors(dim ring I - d, jacobian I));
)

hunekeAlgorithm = method(TypicalValue=>Ideal)
hunekeAlgorithm(Ideal) := I -> (
	R := ring I;
	J:=ideal 0_R;
	previousI:=ideal 0_R;
	while J!=ideal 1_R do (
	    M:=presentation module I;
	    c:=codim I;
	    n:=rank target M;
	    k:=n-1;
	    if codim minors(n-c,M) >= c+1 then
		return genericSocle(I);
	    while codim minors(n-k,M)!=c do 
		k=k-1;
	    J=I:minors(n-k,M);
	    previousI=I;
	    I=J;
	    );
	return genericSocle(previousI);
)


-----------------------------------
-- ** Part 3: Experimenting **
----------------------------------


------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **DOCUMENTATION** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
beginDocumentation()

doc ///
    Key
        RadicalAlgorithms
    Headline
        Implement and benchmark algorithms by Huneke for computing radicals
    Description
        Text
            This package provides functions for benchmarking and implementing the algorithms
	    described in sections 3-10 of the short paper by Craig Huneke [Hun24], which
	    are derived from Levin's theorem (see [Vas94, Thm 10.3.16]).
        Text
            @UL {
	    {"[Hun24] Huneke, Craig (2024). Computing Radicals with Wolmer Vasconcelos"},
	    {"[Vas94] Vasconcelos, Wolmer V (1994). Arithmetic of blowup algebras. Vol. 195. Cambridge University Press."}
	    }@
    Subnodes
        hunekeAlgorithm
///

doc ///
    Key
        hunekeAlgorithm
	(hunekeAlgorithm, Ideal)
    Headline
        computes the radical using an algorithm of Huneke
    Usage
	hunekeAlgorithm(I)
    Inputs
	I: Ideal
	    an unmixed ideal
    Outputs
	: Ideal
	    the radical of the given ideal
    Description
        Text
	    Let $k$ be a field.
	    Given an ideal $I \subseteq k[x_1, \ldots, x_n]$, this method computes $\sqrt I$
	    using the result of [Hun24, Theorem 10.3].
    References
	[Hun24] Huneke, Craig (2024). Computing Radicals with Wolmer Vasconcelos        
///


------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **TESTS** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------

-- hunekeAlgorithm
-- can we write some tests?
TEST ///

///

end---------------------------------------------------------------------------     

------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **SCRATCH SPACE** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------


--------------------------------
--------------------------------
-- **EXAMPLES / CODE DEMO** --
--------------------------------
--------------------------------

EXAMPLE \\\
--
\\\


-----------------------------------
--Ayah's sandbox
----------------------------------

restart
--for this class, unmixed is faster
R = ZZ/101[x_(1,1)..x_(2,3)]
M = genericMatrix(R,2,3)
time radical (minors(2,M))^4
time radical((minors(2,M))^4, Strategy => Unmixed)
I = (minors(2,M))^4 -- 15 gens
I' = ideal (random(I_*))_0;
time I : ideal(jacobian(I'))
radical I'
--for this class, taking the colon with the jacobian
--of a single generator gives the correct radical

R = ZZ/101[a,b,c]
J = ideal(a^124-b^124,a^123*b-c^124);
time radical I
time radical(J, Strategy => Unmixed)
time J: ideal(jacobian(ideal J_1))

------------------------------------
--Development Section
------------------------------------


------------------------------------------------------------------------------------------------------------
-- test: computing radicals with colon of partial minors of the jacobian of a subset of generators
------------------------------------------------------------------------------------------------------------

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
    );


R = QQ[x_1..x_28];
--w = {3, 6, 1, 8, 7, 2, 5, 4};
--I = trim skewSymmDetIdeal(R, w);  -- 21 mingens, not radical

w = {3, 8, 1, 7, 6, 5, 4, 2};
I = trim skewSymmDetIdeal(R, w);  -- 34 mingens, not radical

radI = radical I;

needsPackage "FastMinors"

K = ideal( (random I_*)_{0..4} );

codim minors(4, jacobian K)

#select(
    apply(
    toList(1..100), i -> 
    radI == (I : chooseGoodMinors(9, 4, jacobian K))
    ),
    B -> B
)


restart

-- patch ideals of regular nilpotent Hessenberg varieties
needsPackage "MatrixSchubert";
needsPackage "FastMinors";

patchIdealFromMatrix = (X, h) -> (
    -- X, matrix; h, Hessenberg function
    n := #h;
    idealGenerators := flatten(for i from 1 to n list ( for j from 1 to n list (if i > h#(j-1) then X_(i-1,j-1) else 0 ) ));
    return ideal(delete(0, idealGenerators));
    )

patchIdeal = (w, h) -> (
    W := permToMatrix w;
    X := (inverse(W*M)) * N * (W*M);
    patchIdealFromMatrix(X, h)
    )

R = QQ[x_21,x_31,x_32,x_41..x_43];
M = matrix{
    {1, 0, 0, 0},
    {x_21, 1, 0, 0},
    {x_31, x_32, 1, 0},
    {x_41, x_42, x_43, 1}
    }

N = matrix{
    {0,1,0,0},
    {0,0,1,0},
    {0,0,0,1},
    {0,0,0,0}
    }

h = {2, 3, 4, 4}
S4 = permutations toList (1..4);

patchIdeals = apply(S4, w -> patchIdeal(w, h));
primesPatchIdeals = apply(patchIdeals, decompose);

I = trim (intersect((random patchIdeals)_{0..1}))^2;
radI = radical I;
I == radI  -- ensure that this is false

J = chooseGoodMinors(3, 3, jacobian I);
J == radI

------------------------------------------------------------------------



-- this is very rough and not tested at all


needsPackage "AnalyzeSheafOnP1";  -- don't really need this, just being lazy for now

radicalAlgTest = I -> (
    -- assume I is an ideal in a polynomial ring
    R := ring I;
    supp := support I;
        
    R' := (coefficientRing R) monoid([supp]);
    I' := sub(I, R');
    M' := comodule I';
    jacI := jacobian I';
    
    n := #supp;
    k := min(n, #(I'_*));
    c := codim I';
    
    while k > c do (
	J := minors(min(n, #(I'_*)), jacI);

	-- check that at least one minor is a nonzero divisor mod I', if not, decrement k and go to top of loop

	-- get first nonzero divisor
	nonzerodiv := select(1, J_*, j -> if isNZD(j, M'));

	if nonzerodiv == {} then (
	    k := k-1;
	    continue;
	    )
	break;
	)

    I1 := I' : (ideal first nonzerodiv) ;
    c1 := codim I1;

    if I1 == I1:(minors(n-c1, jacobian I1)) then return sub(I1, R);
    return radicalAlgTest(sub(I1, R));
    )



------------------------------------------------------------------------








restart
uninstallPackage "RadicalAlgorithms"
restart
installPackage "RadicalAlgorithms"
restart
needsPackage "RadicalAlgorithms"
elapsedTime check "RadicalAlgorithms"
viewHelp "RadicalAlgorithms"
