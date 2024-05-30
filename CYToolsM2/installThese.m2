restart
uninstallAllPackages()
restart
installPackage "IntegerEquivalences"
restart
installPackage "StringTorics"
restart
installPackage "DanilovKhovanskii"

-- Currently, GV invariants code not functional here, until we can get computeGV compiled and placed here...
-- I have commented out a few tests that use this.  All the following tests run.
restart
check "IntegerEquivalences"
check "StringTorics"
check "DanilovKhovanskii"
