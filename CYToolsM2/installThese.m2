-- To install StringTorics:
-- Do these lines in order.
restart
uninstallAllPackages()
restart
installPackage "IntegerEquivalences"
restart
installPackage "StringTorics"
check "IntegerEquivalences"
check "StringTorics"

-- At this point, GV invariants functions won't work yet, but everything else should.
-- To get GV invariants going, compile the code in ComputeGV
-- and place computeGV on your PATH.

-- to install DanilovKhovanskii
restart
uninstallPackage "DanilovKhovanskii"
restart
installPackage "DanilovKhovanskii"
restart
check "DanilovKhovanskii"

-- to install the PALP interface
-- for now, need to install PALP on your computer to use this
-- (place the palp executables, e.g. poly.x, cws.x, on your PATH)
-- then:
restart
uninstallPackage "PALPInterface"
restart
needsPackage "PALPInterface"
restart
installPackage "PALPInterface" -- no documentation or tests at all
check "PALPInterface"

