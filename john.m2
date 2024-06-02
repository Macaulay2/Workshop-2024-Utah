-*
Run the following once to install the packages
*-

restart
uninstallPackage "Truncations"
uninstallPackage "Complexes"
uninstallPackage "Varieties"
debug installPackage("Truncations", FileName => "/Users/John/Documents/Github/Workshop-2024-Utah/packages/Truncations.m2")
debug installPackage("Complexes", FileName => "/Users/John/Documents/Github/Workshop-2024-Utah/packages/Complexes.m2")
debug installPackage("Varieties", FileName => "/Users/John/Documents/Github/Workshop-2024-Utah/packages/Varieties.m2")
check "Varieties"

-* 
Run the next thing to make sure quick changes in the packages are reflected here.
*-

needsPackage "Truncations"
needsPackage "Complexes"
needsPackage "Varieties"

-*
Now we can test code here.
*-

R = QQ[x,y,z]

K = koszulComplex vars R
sheafK = sheaf K