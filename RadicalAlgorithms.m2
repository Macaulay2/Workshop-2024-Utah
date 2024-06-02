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
	{Name => "Adam LaClair", 
            Email => "alaclair@purdue.edu", 
            HomePage => "https://sites.google.com/view/adamlaclair/home"},
	{Name => "Trung Chau", 
            Email => "chauchitrung1996@gmail.com"},
	{Name => "Mike Cummings",
	    Email => "cummim5@mcmaster.ca",
	    HomePage => "https://math.mcmaster.ca/~cummim5/"},
	{Name => "David Eisenbud",
	    Email => "de@berkeley.edu"},
	 {Name => "Manohar Kumar",
	     Email => "manhar349@gmail.com"}
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
