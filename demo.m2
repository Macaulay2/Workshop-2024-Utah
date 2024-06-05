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
