uninstallPackage "ToricExtras"
restart
installPackage "ToricExtras"
check "ToricExtras"
restart

-- demooooooooooooooo
needsPackage "ToricExtras"

-- linear series


-- projectivation of line bundles


-- Batyrev classification for smooth, projective toric varieties of Picard rank 3
V = batyrevConstructor({1,1,1,1,1}, {0}, {})
isWellDefined V
dim V
isSmooth V
isProjective V
picardGroup V


V = batyrevConstructor({2,1,1,1,1}, {1}, {})
isWellDefined V
dim V
isSmooth V
isProjective V
picardGroup V
isFano(V)

degreeOfThreefold = X -> (
    K := toricDivisor X;
    c := chern(1, OO (-K));
    integral (c*c*c)
)

degreeOfThreefold V
