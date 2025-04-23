loadPackage ("NumberFields", Reload=>true)




-- This is an example of a number field
R = QQ[w,v]/ ideal(w^3-2,v^2+v+1)
R = (simpleExtension R)_0 
G = galoisGroup R
fixedFields R
getNormalSubgroups(G_2)

