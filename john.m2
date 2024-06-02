-* This should serve as a generic "development" file for creating packages.

On the package file, you place the following: 

uninstallPackage "Varieties"
restart
loadPackage("Truncations", FileName => currentDirectory() | "Truncations.m2", Reload => true)
loadPackage("Complexes",   FileName => currentDirectory() | "Complexes.m2",   Reload => true)
loadPackage("Varieties",   FileName => currentDirectory() | "Varieties.m2",   Reload => true)
installPackage("Varieties",   FileName => currentDirectory() | "Varieties.m2")
viewHelp "Varieties"

restart
debug needsPackage "Varieties"
check "Varieties"

These are all the commands that should be run to register changes to the package. That is, open Macaulay2 FROM `Varieties.m2'. If you did it currently, running currentDirectory() will give you the directory where Varieties.m2 is located.

-- The following puts in a debugger and makes the error codes actually useful. The default mode I guess is "not useful".
-- Lower the error depth if it is still not useful -- lower the error depth, the more information you get.
*- 
errorDepth = 2

-*
Now we can test code here.
*-

R = QQ[x,y,z]

K = koszulComplex vars R
sheafK = sheaf K

complex {sheaf vars R}