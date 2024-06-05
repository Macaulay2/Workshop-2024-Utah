uninstallPackage "ToricExtras"
restart
installPackage "ToricExtras"
check "ToricExtras"
restart

-- demooooooooooooooo
needsPackage "ToricExtras"

-- linear series
P2 = toricProjectiveSpace 2;
D = 3 * P2_0; -- 3th Veronese Divisor of P2
linSeries = toricLinearSeries D;
monomials linSeries
isComplete linSeries
isBasepointFree linSeries
f = toricMap linSeries;
matrix f
isWellDefined f
ideal f == idealOfImage(linSeries, TargetRing => ring target f)

P1 = toricProjectiveSpace 1;
X = P1 ** P1;
T = ring X;
mons1 = {T_{0,2,2,0},T_{0,2,0,2},T_{2,0,0,2},T_{2,0,2,0}}
linSeries1 = toricLinearSeries mons1;
isComplete linSeries1
isBasepointFree linSeries1
g = toricMap linSeries1
isWellDefined g


-- projectivation of line bundles


-- Batyrev classification for smooth, projective toric varieties of Picard rank 3
V = batyrevConstructor({1,1,1,1,1}, {0}, {})
