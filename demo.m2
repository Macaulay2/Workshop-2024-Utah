loadPackage "Divisor"
loadPackage "RuledSurfaces"

P1 = Proj ZZ/101[x,y]
E = OO_P1^2
PE = projectiveBundle E
peek PE
X11 = imageOfLinearSeries(PE,OO_P1(1),1)
ideal X11
X22 = imageOfLinearSeries(PE,OO_P1(2),2)
ideal X22
(dim X22, dim ambient X22, degree X22)
divisorToLineBundle(-canonicalDivisor(ring X22,IsGraded=>true))

S = (ZZ/101)[x,y,z]
C = Proj quotient ideal(x^3+y^3+z^3)
E = OO_C^1++OO_C(1)
PE = projectiveBundle E
findGlobalGeneratorsOfTwist PE 

ideal imageOfLinearSeries(PE,OO_C(0),1)
ideal imageOfLinearSeries(PE,OO_C(1),1)
ideal minimalEmbedding PE

