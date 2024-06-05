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
isComplete linSeries2
isBasepointFree linSeries2
g = toricMap linSeries2
isWellDefined g


-- projectivization of line bundles
rayListPP1 = {{1}, {-1}}
coneListPP1 = {{1}, {0}}
PP1 = normalToricVariety (rayListPP1, coneListPP1)
D0PP1 = toricDivisor ({0, 0}, PP1)
D1PP1 = toricDivisor ({0, 7} , PP1)
testH7 = projectivizationOfBundle({D0PP1, D1PP1})
isWellDefined testH7
H7 = hirzebruchSurface 7

rays testH7
rays H7
max testH7
max H7

isSmooth testH7
isProjective testH7

rayListPP2 = {{1, 0}, {0, 1}, {-1, -1}}
coneListPP2 = {{0, 1}, {1, 2}, {2, 0}}
PP2 = normalToricVariety (rayListPP2, coneListPP2)
D0PP2 = toricDivisor ({9, 3, 2}, PP2)
D1PP2 = toricDivisor ({1, 4, 7}, PP2)
Y = projectivizationOfBundle({D0PP2, D1PP2})
isWellDefined Y
rays Y
max Y

isSmooth Y
isProjective Y

-- Batyrev classification for smooth, projective toric varieties of Picard rank 3
-- This constructs the unique smooth toric del Pezzo surface of Picard rank 3, 
-- which is the Bl_1(P^1 x P^1) = Bl_2(P^2)
V = batyrevConstructor({1,1,1,1,1}, {0}, {})
isWellDefined V
dim V
isSmooth V
isProjective V
picardGroup V

-- This constructs a smooth toric Fano of Picard rank 3 which is 3-29 on Fanography. 
-- It can be described as the blowup of Bl_1 P^3 along a line in the exceptional divisor
V = batyrevConstructor({2,1,1,1,1}, {1}, {})
isWellDefined V
dim V
isSmooth V
isProjective V
picardGroup V
isFano V

degreeOfThreefold = X -> (
    K := toricDivisor X;
    c := chern(1, OO (-K));
    integral (c*c*c)
)

degreeOfThreefold V
