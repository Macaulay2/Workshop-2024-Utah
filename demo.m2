uninstallPackage "ToricExtras"
restart
installPackage "ToricExtras"
check "ToricExtras"
restart

-- demooooooooooooooo
needsPackage "ToricExtras"

-- linear series


-- projectivation of line bundles

rayListPP1 = {{1} , {-1}}
coneListPP1 = {{1}, {0}}
PP1 = normalToricVariety (rayListPP1, coneListPP1)
D0PP1 = toricDivisor ( { 0 , 0}, PP1)
D1PP1 = toricDivisor ( {0 , 7} , PP1)
H7 = projectivizationOfBundle({D0PP1, D1PP1})
rays H7
max H7

rayListPP2 = {{1 , 0}, {0 , 1}, {-1, -1}}
coneListPP2 = {{0, 1}, {1, 2}, {2 , 0}}
PP2 = normalToricVariety (rayListPP2, coneListPP2)
D0PP2 = toricDivisor ( { 9, 3 , 2}, PP2)
D1PP2 = toricDivisor ( {1 , 4, 7} , PP2)
Y = projectivizationOfBundle({D0PP2, D1PP2})
rays Y
max Y

-- Batyrev classification for smooth, projective toric varieties of Picard rank 3
V = batyrevConstructor({1,1,1,1,1}, {0}, {})
