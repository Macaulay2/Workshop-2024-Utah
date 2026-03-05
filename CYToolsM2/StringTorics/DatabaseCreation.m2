----------------------------------
-- Code for creating data bases --
----------------------------------
label KSEntry := ZZ => (ks) -> (
    str := toString ks;
    ans := regex("id:([0-9]+)", str);
    if ans === null or #ans != 2 then 
    null 
    else 
    value substring(str, ans#1#0, ans#1#1)
    )
  
hodgeNumbers = method()
hodgeNumbers KSEntry := (ks) -> (
    str := toString ks;
    ans := regex("H:([0-9]+),([0-9]+)", str);
    if #ans != 3 then error "expected 3 matches";
    (value substring(str, ans#1#0, ans#1#1),
        value substring(str, ans#2#0, ans#2#1))
    )
------------------------------------------------
-- Combining several Database files into one ---
------------------------------------------------
combineCYDatabases = method()
combineCYDatabases(Database, Database) := (db1, db2) -> (
    -- appends all keys of db2 to db1
    for k in keys db2 do db1#k = db2#k;
    )
combineCYDatabases(String, String) := (dbname1, dbname2) -> (
    -- appends all keys of db2 to db1
    db1 := openDatabaseOut dbname1;
    db2 := openDatabase dbname2;
    combineCYDatabases(db1, db2);
    close db1;
    close db2;
    )
combineCYDatabases List := (dbL) -> (
    -- either all elements are String filename's or are databases.
    for i from 1 to #dbL-1 do combineCYDatabases(dbL#0, dbL#i)
    )

addToCYDatabase = method(Options => {NTFE => true, "CYs" => true})

-- This function adds the ReflexivePolytope 'ks' to the database, if it is not there yet.
-- Actually, it only looks at the ID label in the 'ks' entry, not at the polytope itself.
-- Under default conditions, all NTFE triangulations are found, and all corresponding CY's
-- are placed into the data base.
-- This function returns the ReflexivePolytope found or created.
addToCYDatabase(String, KSEntry) := ReflexivePolytope => opts -> (dbfilename, ks) -> (
    lab := label ks;
    F := openDatabaseOut dbfilename;
    if not F#?(toString lab) then (
        << "computing for polytope " << lab << endl;
        Q := reflexivePolytope(ks, ID => lab); -- note that the polytope data is really that of the dual to topes#i.
        -- now fill it with data we want
        computeBasics Q;
        -- now write it
        F#(toString lab) = dump Q;
        )
    else (
        Q = reflexivePolytope F#(toString lab);
        );
    close F;
    if opts#"CYs" then addToCYDatabase(dbfilename, Q, NTFE => opts.NTFE);
    Q
    )

addToCYDatabase(String, String, Sequence) := String => opts -> (dbfilenamePrefix, topesFilename, lohi) -> (
    -- dbfilenamePrefix will include lo,hi in the name of the created database.
    -- creates (or appends to) a database, and returns the name of the database file.
    topes := value get topesFilename;
    (lo,hi) := lohi;
    if lo < 0 then error "expected range of non-negative integers";
    if hi >= #topes then hi = #topes-1; -- last one
    mytopes := take(topes, toList lohi);
    dbname := dbfilenamePrefix | "-range-" | lohi#0 | "-" | lohi#1 | ".dbm";
    t := elapsedTiming addToCYDatabase(dbname, mytopes, opts);
    << "filename " << dbname << " has been constructed in " << t#0 << " sec" << endl;
    dbname
    )

addToCYDatabase(String, String, ZZ, ZZ) := opts -> (dbfilenamePrefix, topesFilename, whichpart, numparts) -> (
    -- dbfilenamePrefix will include lo,hi in the name of the created database.
    -- creates (or appends to) a database, and returns the name of the database file.
    topes := value get topesFilename;
    nPerPart := ceiling(#topes / (numparts + 0.0)); -- all but the last...
    -- now recompute #parts...
    
    lo := whichpart * nPerPart;
    hi := (whichpart+1)  * nPerPart - 1;
    if lo < 0 then error "expected range of non-negative integers";
    if hi >= #topes then hi = #topes-1; -- last one
    << "doing part=" << whichpart << " #parts=" << numparts << " range:" << lo << " " << hi << endl;
    mytopes := take(topes, {lo,hi});
    dbname := dbfilenamePrefix | "-" | whichpart | "-of-" | numparts | ".dbm";
    t := elapsedTiming addToCYDatabase(dbname, mytopes, opts);
    << "filename " << dbname << " has been constructed in " << t#0 << "s" << endl;
    dbname
    )

createGroups = (numtotal, numgroups) -> (
    q := numtotal // numgroups;
    r := numtotal % numgroups;
    --print (q,r);
    set1 := for i from 0 to r-1 list (i*(q+1), i*(q+1) + q);
    set2 := for j from 0 to numgroups-r-1 list (r*(q+1) + j*q, r*(q+1) + j*q + q-1);
    join(set1, set2)
    )

cmdLine = ///M2 --silent --stop -e 'needsPackage "StringTorics"' -e 'lohi = (LO,HI)' -e 'addToCYDatabase("DBNAMEPREFIX", "TOPESFILE", lohi)' -e 'exit 0' &///
createM2Lines = (dbfilenamePrefix, topesFilename, numtotal, numgroups) -> (
    sets := createGroups(numtotal, numgroups);
    cmds := for s in sets list (
        replace("TOPESFILE", topesFilename,
        replace("DBNAMEPREFIX", dbfilenamePrefix,
        replace("HI", toString s#1, 
            replace("LO", toString s#0, cmdLine)))));
    concatenate between("\n", cmds)
    )

///
restart
debug needsPackage "StringTorics"
createM2Lines("AFile", "BFile", 40, 11)
createGroups(17101, 100)
///

addToCYDatabase(String, ReflexivePolytope) := opts -> (dbfilename, Q) -> (
    -- This version also finds "moriConeCap" which is a cone containing the actual mori cone: it is the
    -- intersection of all mori cones coming from triangulations equivalent to the given one.
    Xs := findAllCYs Q; -- TODO: check: is findALlCYs still correct.
    -- << "  " << #Xs << " triangulations total" << endl;
    -- if opts.NTFE then (
    --     elapsedTime H := partition(restrictTriangulation, Xs);
    --     << "  " << #(keys H) << " NTFE triangulations" << endl;
    --     Xs = (keys H)/(k -> H#k#0); -- only take one triangulation that matches
    --     Xs = for k in keys H list (
    --         X := H#k#0;
    --         setToricMoriConeCap(X, H#k);
    --         X
    --         )
    --     -- let's relabel these Xs?
    --     );
    F := openDatabaseOut dbfilename;
    for X in Xs do (
        setToricMoriConeCap X;
        computeIntersectionNumbers X; -- this should load all of the data we want
        F#(toString label X) = dump X;
        );
    close F;    
    )

addToCYDatabase(String, List) := opts ->(dbfilename, topes) -> (
    for tope in topes do addToCYDatabase(dbfilename, tope, opts);
    )

addToCYDatabase(String, String, List) := opts ->(dbfilename, dbQfilename, topeLabels) -> (
    for lab in topeLabels do (
        << "polytope " << lab << endl;
        Q := cyPolytope(dbQfilename, lab);
        addToCYDatabase(dbfilename, Q, opts);
        );
    )

-- Options => {
--         Limit => 100000,
--         "CYs" => true})

------------------------------
-- Top level functions -------
-- createCYDatabaseFiles  ----
-- combineCYDatabaseFiles ----
------------------------------
createCYDatabaseFiles = method()
combineCYDatabaseFiles = method()

createCYDatabaseFiles(String, String, ZZ) := (prefix, topesfile, ncores) -> (
    -- set file names based on prefix
    topes := value get topesfile;
    str := createM2Lines(prefix, topesfile, #topes, ncores);
    << "-- starting " << ncores << " jobs ------" << endl;
    run str; -- hard to tell when this is done!
    )

-- run this after all "ncores" processes have completed.
combineCYDatabaseFiles(String, String, ZZ) := (prefix, topesfile, ncores) -> (
    -- set file names based on prefix
    fileglob := prefix|"-range*";
    dbname := prefix|".dbm";
    topes := value get topesfile;
    files := lines(get("!ls "|fileglob)); -- uses output of `run str`.
    combineCYDatabases({dbname} | files);
    run ("rm "|fileglob); -- remove all these constructed files (but not dbname!)

   -- now we test if all the polytopes were actually included.
   F := openDatabase dbname;
   handled := for x in sort keys F list (y := value x; if instance(y, ZZ) then y else continue);
   if handled =!= toList(0..#topes-1) then (
       << "----WARNING: not all polytopes were constructed for some reason----" << end;
       );
   close F;
   dbname
   )

----------------------------------------
-- Reading already existing databases --
----------------------------------------
readCYDatabase = method(Options => {Ring => null})
readCYDatabase String := Sequence => opts -> (dbname) -> (
    F := openDatabase dbname;
      labs := (keys F)/value;
      Qlabels := sort select(labs, lab -> instance(lab, ZZ));
      Xlabels := sort select(labs, lab -> instance(lab, Sequence));
      Qs := hashTable for lab in Qlabels list lab => reflexivePolytope F#(toString lab);
      Xs := hashTable for lab in Xlabels list lab => calabiYau(F#(toString lab), i -> Qs#i, opts);
    close F;
    (Qs, Xs)
    )

readCYPolytopes = method()
readCYPolytopes String := HashTable => dbname -> (
    F := openDatabase dbname;
      labs := (keys F)/value;
      Qlabels := sort select(labs, lab -> instance(lab, ZZ));
      Qs := hashTable for lab in Qlabels list lab => reflexivePolytope F#(toString lab);
    close F;
    Qs
    )

readCYs = method(Options => {Ring => null})
readCYs(String, HashTable) := HashTable => opts -> (dbname, Qs) -> (
    F := openDatabase dbname;
      labs := (keys F)/value;
      Xlabels := sort select(labs, lab -> instance(lab, Sequence));
      Xs := hashTable for lab in Xlabels list lab => calabiYau(F#(toString lab), i -> Qs#i, opts);
    close F;
    Xs
    )

-------------------------------------------------------
-- Read one example from a database or database file --
-------------------------------------------------------
cyPolytope(String, ZZ) := ReflexivePolytope => opts -> (dbfilename, topeid) -> (
    db := openDatabase dbfilename;
    Q := cyPolytope(db, topeid, opts);
    close db;
    Q
    )

cyPolytope(Database, ZZ) := ReflexivePolytope => opts -> (db, topeid) -> (
    k := toString topeid;
    if not db#?k then error("polytope with label "|k|" does not exist");
    reflexivePolytope(db#k, opts)
    )

-- Check: this is not quite correct.
calabiYau(Database, ReflexivePolytope, Sequence) := CalabiYauInToric => opts -> (db, Q, lab) -> (
    -- lab should be (polytopelab, triangulationlabel).
    -- polytopelab should match label of Q.
    if first lab =!= label Q then error "incorrect label";
    k := toString lab;
    if not db#?k then error("polytope with label "|k|" does not exist");
    calabiYau(db#k, lab -> Q, opts)
    )

calabiYau(Database, Sequence) := CalabiYauInToric => opts -> (db, lab) -> (
    -- lab should be (polytopelab, triangulationlabel).
    -- first retrieve ReflexivePolytope, and then CalabiYauInToric.
    if #lab < 2 then error "expected well-formed label";
    Q := cyPolytope(db, first lab);
    k := toString lab;
    if not db#?k then error("CY with label "|k|" does not exist");
    calabiYau(db#k, lab -> Q, opts)
    )

calabiYau(String, ReflexivePolytope, Sequence) := CalabiYauInToric => opts -> (dbfilename, Q, lab) -> (
    db := openDatabase dbfilename;
    X := calabiYau(db, Q, lab, opts);
    close db;
    X
    )

calabiYau(String, Sequence) := CalabiYauInToric => opts -> (dbfilename, lab) -> (
    db := openDatabase dbfilename;
    X := calabiYau(db, lab, opts);
    close db;
    X
    )


----------------------------------------------------------------
-- Reading and writing polytopes, simplices, cy_classes files --
----------------------------------------------------------------
readPolytopes = method()
readPolytopes String := filename -> (
    contents := lines get filename;
    hashTable for L in contents list (
        v := toList value L;
        lab := v#0;
        i := 1;
        pts := while i+3 <= #v list (
            ans := for j from 0 to 3 list v#(i+j);
            i = i + 4;
            ans);
        lab => pts
        )
    )

readSimplices = method()
readSimplices String := filename -> (
    contents := lines get filename;
    hashTable for L in contents list (
        v := toList value L;
        labX := v#0;
        labQ := v#1;
        i := 2;
        simplices := while i+3 <= #v list (
            ans := for j from 0 to 3 list v#(i+j)-1;
            i = i + 4;
            ans);
        labX => {labQ, simplices}
        )
    )

readEquivalences = method()
readEquivalences String := List => filename -> (
    contents := lines get filename;
    for L in contents list (
        v := value ("{"|L|"}");
        i := 0;
        equivsets := while i < #v list (
            ans := while i < #v and v#i != -1 list (a := v#i; i=i+1; a);
            if i < #v then i = i+1;
            ans
            );
        equivsets
        )
    )

cyPolytope(HashTable, ZZ):= CYPolytope => opts -> (vertexData, ind) -> (
    cyPolytope(transpose matrix vertexData#ind, ID => ind)
    )
-- This function should not change the order of points?  But it does.


makeCY(ZZ, Sequence, MutableHashTable, MutableHashTable) := opts -> (labX, PSs, Qs, Xs) -> (
    error "deprecated, either rewrite me, or remove me";
    if Xs#?labX then return Xs#labX;
    (Ps, Ss) := PSs;
    polytopeid := Ss#labX#0;
    P := Ps#polytopeid;
    S := Ss#labX;
    if not Qs#?polytopeid then (
        Qs#polytopeid = cyPolytope(P, ID => polytopeid);
        );
    Q := Qs#polytopeid;
    -- now place into the cache the translation for rays?
    ---translate1 := hashTable for i from 0 to #Ps - 1 list Ps#i => i;
    translate2 := hashTable for i from 0 to #(rays Q) - 1 list (rays Q)#i => i;
    fromOldToNew := for i from 0 to #P-1 list translate2#(P#i);
    tri := sort for T in S#1 list sort for t1 in T list fromOldToNew#t1;
    X := calabiYau(Q, tri, ID => labX, Ring => opts#Ring);
    Xs#labX = X;
    X
    )

end--
-----------------------------------------
-- mongodb databases? -------------------
-----------------------------------------
-- Here is some experimentation
  restart
  needsPackage "StringTorics"
  polytopes6 = kreuzerSkarke(6, Limit => 100000);
  #polytopes6 == 17101

  polytopes7 = kreuzerSkarke(7, Limit => 100000);
  #polytopes7 == 50376

  polytopes8 = kreuzerSkarke(8, Limit => 500000);
  #polytopes8 == 128165

  polytopes9 = kreuzerSkarke(9, Limit => 1000000);
  #polytopes9 == 285929

  polytopes10 = kreuzerSkarke(10, Limit => 2000000);
  #polytopes10 == 568078

  elapsedTime mats6 = polytopes6/matrix;
  mats6/(m -> numcols m)//tally  

  elapsedTime mats7 = polytopes7/matrix;  -- 212 sec!  hmm doing it again: 95 sec, 87 sec
  elapsedTime mats7/(m -> numcols m)//tally  

  elapsedTime mats7a = mats7/(m -> entries transpose m); -- 9 sec, now 1.9 sec, now 68 sec!!
  elapsedTime mats7b = mats7/(m -> entries m); -- 1.6 sec, 1.7 sec
  elapsedTime mats7c = mats7/(m -> transpose m); -- 15 sec, .05 sec

  elapsedTime mats7a/(m -> transpose matrix m); -- 438 sec!!
  
  
  elapsedTime mats7a/(m -> matrix transpose m); -- 848 sec!!
  elapsedTime mats7a/(m -> matrix m); --  sec!!

  restart
  needsPackage "StringTorics"
  polytopes7 = kreuzerSkarke(7, Limit => 100000);
  #polytopes7 == 50376
  elapsedTime mats7 = polytopes7/matrix;  -- 212 sec!  hmm doing it again: 95 sec, 87 sec, 80 sec
  elapsedTime mats7/(m -> numcols m)//tally  
  elapsedTime mats7a = mats7/(m -> entries transpose m); -- 9 sec, now 1.9 sec, now 68 sec!! 56 sec

  GC_INITIAL_HEAP_SIZE=40G M2 ....
  polytopes7 = kreuzerSkarke(7, Limit => 100000);
  #polytopes7 == 50376
  elapsedTime mats7 = polytopes7/matrix;  -- 11 sec
  elapsedTime mats7/(m -> numcols m)//tally  
  elapsedTime mats7a = mats7/(m -> entries transpose m); -- 3 sec
  elapsedTime mats7a/(m -> matrix m); -- 3 sec

