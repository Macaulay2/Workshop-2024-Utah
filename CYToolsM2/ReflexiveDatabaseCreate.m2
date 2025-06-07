debug needsPackage "ReflexivePolytopesDB"

makeMongoDoc = method()
makeMongoDoc KSEntry := (tope) -> (
    -- example we want to extract from: "4 5  M:376 5 N:14 5 H:7,271 [-528]"
    re := "([0-9]+) +([0-9]+) +M:([0-9]+) +([0-9]+) +N:([0-9]+) +([0-9]+) +H:([0-9])+,([0-9]+) +";
    s := description tope;
    vals := regex(re, s);
    subs := for x from 0 to #vals-1 list substring(vals#x, s);
    verts := toString new Array from ((entries transpose matrix tope)/(x -> new Array from x));
    str := "{vertex_count:$VC, facet_count:$FC, point_count:$PC, dual_point_count:$DPC, h11:$H11, h12:$H12, vertices:$VERTS}";
    str = replace("\\$VC", subs#4, str);
    str = replace("\\$FC", subs#6, str);
    str = replace("\\$PC", subs#3, str);
    str = replace("\\$DPC", subs#5, str);
    str = replace("\\$H11", subs#7, str);
    str = replace("\\$H12", subs#8, str);
    str = replace("\\$VERTS", verts, str);
    str
    )

makeMongoDoc(String, List) := (filename, topes) -> (
    F := openOut filename;
    F << "db.reflexive4d.insertMany([" << endl;
    for i from 0 to #topes-1 do (
        t := topes#i;
        F << makeMongoDoc t << (if i === #topes-1 then "" else ",") << endl;
        );
    F << "])" << endl;
    close F
    )

createMongoInserts = method()
createMongoInserts(String, String) := (infilename, outfilename) -> (
    if fileExists outfilename then error "cannot write over existing file";
    contents := get infilename;
    topes := parseKSRawString contents;
    makeMongoDoc(outfilename, topes)
    )

-- conts = get "~/kreuzerSkarke/v05";


-- topes = parseKSRawString conts;
-- topes/makeMongoDoc;
-- elapsedTime makeMongoDoc("foov5", topes)




-- ///
-- tope = KSEntry "4 5  M:376 5 N:14 5 H:7,271 [-528]
--    1   0   0   9 -27
--    0   1   1  10 -26
--    0   0   3  12 -24
--    0   0   0  18 -18
-- "    
-- header tope
-- methods KSEntry
-- description tope
-- re = "([0-9]+) +([0-9]+) +M:([0-9]+) +([0-9]+) +N:([0-9]+) +([0-9]+) +H:([0-9])+,([0-9]+) +"
-- regex(re, description tope)
-- for x in drop(oo, 1) list value substring(x, description tope)
-- ///

end--

restart
load "~/utah/CYToolsM2/ReflexiveDatabaseCreate.m2"
createMongoInserts("~/kreuzerSkarke/v05", "foov05");
time mongosh mike --file foov05

elapsedTime createMongoInserts("~/kreuzerSkarke/v06", "foov06");
time mongosh mike --file foov06


restart
debug needsPackage "StringTorics"
topes = kreuzerSkarke(2, Limit => 1000);
topes_0
ks = topes_0
lab = 0
Q = cyPolytope(ks, ID => lab);
toJson Q
elapsedTime Qs = for t in topes list cyPolytope t;
"foo2" << elapsedTime toJson Qs << close;

concatenate oo
