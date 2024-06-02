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
-- ++ means tested
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

 hunekeAlgorithm := (I) -> (
    J := ideal(0_R);
    previous := ideal(0_R);
    while (J != ideal(1_R)) do (
	M := presentation module I;
	c := codim I;
	n := rank target M;
	k := n-1;
	while codim minors(n-k,M) != c do (
	    k = k-1;
	    );
	if (k == c) then (
	    print("generically CI");
	    return null;
	    )
	else if (k < c) then (
	    print("k<c, what now?");
	    return null;
	    )
	else (
	    J = I : minors(n-k,M);
	    previous = I;
	    I = J;
	    );
    );
    return previousI
 )





-----------------------------------
-- ** Part 3: Experimenting **
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
///

doc ///
    Key
        hunekeAlgorithm
	hunekeAlgorithm(Ideal)
    Headline
        Implement the algorithm of Huneke using minors of the presentation matrix
    Usage
	hunekeAlgorithm(I)
    Inputs
	I: Ideal
	    an unmixed ideal
    Outputs
	:
    Description
        Text
            Given an ideal, it gives the outputs of Huneke's algorithm.
///


------------------------------------------------------------------------------
------------------------------------------------------------------------------
-- **TESTS** --
------------------------------------------------------------------------------
------------------------------------------------------------------------------
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

--trung test comment
------------------------------------
--Development Section
------------------------------------

restart
uninstallPackage "RadicalAlgorithms"
restart
installPackage "RadicalAlgorithms"
restart
needsPackage "RadicalAlgorithms"
elapsedTime check "RadicalAlgorithms"
viewHelp "RadicalAlgorithms"
