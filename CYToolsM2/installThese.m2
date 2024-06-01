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

-- PALP interface
-- for now, need to install PALP on your computer to use this.
restart
uninstallPackage "PALPInterface"
restart
needsPackage "PALPInterface"
restart
installPackage "PALPInterface" -- no documentation or tests at all
check "PALPInterface"

-- We can make a similar interface to GV invariants?
-- maybe called GromovWitten or GVInvariants.
