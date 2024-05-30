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

-- This function adds the CYPolytope 'ks' to the database, if it is not there yet.
-- Actually, it only looks at the ID label in the 'ks' entry, not at the polytope itself.
-- Under default conditions, all NTFE triangulations are found, and all corresponding CY's
-- are placed into the data base.
-- This function returns the CYPolytope found or created.
addToCYDatabase(String, KSEntry) := CYPolytope => opts -> (dbfilename, ks) -> (
    lab := label ks;
    F := openDatabaseOut dbfilename;
    if not F#?(toString lab) then (
        << "computing for polytope " << lab << endl;
        Q := cyPolytope(ks, ID => lab); -- note that the polytope data is really that of the dual to topes#i.
        -- now fill it with data we want
        basisIndices Q; -- compute them
        isFavorable Q; -- compute h11, h21, favorability.
        annotatedFaces Q; -- compute annotated faces
        automorphisms Q;
        automorphismsAsPermutations Q;
        findAllFRSTs Q;
        -- now write it
        F#(toString lab) = dump Q;
        )
    else (
        Q = cyPolytope F#(toString lab);
        );
    close F;
    if opts#"CYs" then addToCYDatabase(dbfilename, Q, NTFE => opts.NTFE);
    Q
    )

processCYPolytopes = method()
processCYPolytopes(String, ZZ, Sequence) := (dbfilenamePrefix, h11, lohi) -> (
    elapsedTime topes := kreuzerSkarke(h11, Limit => 200000);
    mytopes := take(topes, toList lohi);
    dbname := dbfilenamePrefix | "-" | lohi#0 | "-" | lohi#1 | ".dbm";
    elapsedTime createCYDatabase(dbname, mytopes);
    )

--addToCYDatabase = method(Options => {NTFE => false})

-- Delete this older version (which doesn't compute toric mori cone caps)
-- This only adds the CY's coming from Q.
-- addToCYDatabase(String, CYPolytope) := opts -> (dbfilename, Q) -> (
--     elapsedTime Xs := findAllCYs Q; -- TODO: check: is findALlCYs still correct.
--     << "  " << #Xs << " triangulations total" << endl;
--     if opts.NTFE then (
--         elapsedTime H := partition(restrictTriangulation, Xs);
--         << "  " << #(keys H) << " NTFE triangulations" << endl;
--         Xs = (keys H)/(k -> H#k#0); -- only take one triangulation that matches
--         -- let's relabel these Xs
--         );
--     F := openDatabaseOut dbfilename;
--     for X in Xs do (
--         computeIntersectionNumbers X; -- this should load all of the data we want
--         F#(toString label X) = dump X;
--         );
--     close F;    
--     )

addToCYDatabase(String, CYPolytope) := opts -> (dbfilename, Q) -> (
    -- This version also finds "moriConeCap" which is a cone containing the actual mori cone: it is the
    -- intersection of all mori cones coming from triangulations equivalent to the given one.
    elapsedTime Xs := findAllCYs Q; -- TODO: check: is findALlCYs still correct.
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
        elapsedTime addToCYDatabase(dbfilename, Q, opts);
        );
    )

createCYDatabase = method(Options => {
        Limit => 100000,
        NTFE => true,
        "CYs" => true})
createCYDatabase(String, ZZ, List) := opts -> (dbfileprefix, h11, range) -> (
    topes := kreuzerSkarke(h11, Limit => opts.Limit);
    (lo, hi) := toSequence range;
    hi = hi-1;
    filename := dbfileprefix | "-range-"|lo|"-"|hi|".dbm";
    addToCYDatabase(filename, topes_{lo..hi}, NTFE => opts.NTFE, "CYs" => opts#"CYs")
    )
-- createCYDatabase = method()

-- createCYDatabase(String, List) := (dbfilename, topes) -> (
--     -- open data base file
--     F := openDatabaseOut dbfilename;
--     -- loop through topes, create CYPolytope, populate it, write it to data base.
--     elapsedTime for i from 0 to #topes - 1 do elapsedTime (
--         lab := label topes_i;
--         if lab === null then lab = i; -- else print "using label";
--         << "computing for polytope " << lab << endl;
--         V := cyPolytope(topes#i, ID => lab); -- note that the polytope data is really that of the dual to topes#i.
--         -- now fill it with data we want
--         basisIndices V; -- compute them
--         isFavorable V; -- compute h11, h21, favorability.
--         annotatedFaces V; -- compute annotated faces
--         automorphisms V;
--         -- now write it
--         F#(toString lab) = dump V;
--         );
--     close F;
--     )


-- addToCYDatabase(String, CYPolytope) := opts -> (dbfilename, Q) -> (
--     elapsedTime Xs := findAllCYs Q; -- TODO: check: is findALlCYs still correct.
--     << "  " << #Xs << " triangulations total" << endl;
--     if opts.NTFE then (
--         elapsedTime H := partition(restrictTriangulation, Xs);
--         << "  " << #(keys H) << " NTFE triangulations" << endl;
--         Xs = (keys H)/(k -> H#k#0); -- only take one triangulation that matches
--         -- let's relabel these Xs
--         );
--     F := openDatabaseOut dbfilename;
--     for X in Xs do (
--         computeIntersectionNumbers X; -- this should load all of the data we want
--         F#(toString label X) = dump X;
--         );
--     close F;    
--     )

-- addToCYDatabase(String, CYPolytope) := opts -> (dbfilename, Q) -> (
--     -- This version also finds "moriConeCap" which is a cone containing the actual mori cone: it is the
--     -- intersection of all mori cones coming from triangulations equivalent to the given one.
--     elapsedTime Xs := findAllCYs Q; -- TODO: check: is findALlCYs still correct.
--     << "  " << #Xs << " triangulations total" << endl;
--     if opts.NTFE then (
--         elapsedTime H := partition(restrictTriangulation, Xs);
--         << "  " << #(keys H) << " NTFE triangulations" << endl;
--         Xs = (keys H)/(k -> H#k#0); -- only take one triangulation that matches
--         Xs = for k in keys H list (
--             X := H#k#0;
--             setToricMoriConeCap(X, H#k);
--             X
--             )
--         -- let's relabel these Xs?
--         );
--     F := openDatabaseOut dbfilename;
--     for X in Xs do (
--         computeIntersectionNumbers X; -- this should load all of the data we want
--         F#(toString label X) = dump X;
--         );
--     close F;    
--     )

-- addToCYDatabase(String, Database, ZZ) := opts -> (dbfilename, topesDB, i) -> (
--     <<  "-- doing polytope " << i << endl;
--     Q := cyPolytope(topesDB#(toString i), ID => i);
--     addToCYDatabase(dbfilename, Q, opts);
--     )

readCYDatabase = method(Options => {Ring => null})
readCYDatabase String := Sequence => opts -> (dbname) -> (
    F := openDatabase dbname;
      labs := (keys F)/value;
      Qlabels := sort select(labs, lab -> instance(lab, ZZ));
      Xlabels := sort select(labs, lab -> instance(lab, Sequence));
      Qs := hashTable for lab in Qlabels list lab => cyPolytope F#(toString lab);
      Xs := hashTable for lab in Xlabels list lab => cyData(F#(toString lab), i -> Qs#i, opts);
    close F;
    (Qs, Xs)
    )

readCYPolytopes = method()
readCYPolytopes String := HashTable => dbname -> (
    F := openDatabase dbname;
      labs := (keys F)/value;
      Qlabels := sort select(labs, lab -> instance(lab, ZZ));
      Qs := hashTable for lab in Qlabels list lab => cyPolytope F#(toString lab);
    close F;
    Qs
    )

readCYs = method(Options => {Ring => null})
readCYs(String, HashTable) := HashTable => opts -> (dbname, Qs) -> (
    F := openDatabase dbname;
      labs := (keys F)/value;
      Xlabels := sort select(labs, lab -> instance(lab, Sequence));
      Xs := hashTable for lab in Xlabels list lab => cyData(F#(toString lab), i -> Qs#i, opts);
    close F;
    Xs
    )

-------------------------------------------------------
-- Read one example from a database or database file --
-------------------------------------------------------
cyPolytope(String, ZZ) := CYPolytope => opts -> (dbfilename, topeid) -> (
    db := openDatabase dbfilename;
    Q := cyPolytope(db, topeid, opts);
    close db;
    Q
    )

cyPolytope(Database, ZZ) := CYPolytope => opts -> (db, topeid) -> (
    k := toString topeid;
    if not db#?k then error("polytope with label "|k|" does not exist");
    cyPolytope(db#k, opts)
    )

-- Check: this is not quite correct.
calabiYau(Database, CYPolytope, Sequence) := CalabiYauInToric => opts -> (db, Q, lab) -> (
    -- lab should be (polytopelab, triangulationlabel).
    -- polytopelab should match label of Q.
    if first lab =!= label Q then error "incorrect label";
    k := toString lab;
    if not db#?k then error("polytope with label "|k|" does not exist");
    calabiYau(db#k, lab -> Q, opts)
    )

calabiYau(Database, Sequence) := CalabiYauInToric => opts -> (db, lab) -> (
    -- lab should be (polytopelab, triangulationlabel).
    -- first retrieve CYPolytope, and then CalabiYauInToric.
    if #lab < 2 then error "expected well-formed label";
    Q := cyPolytope(db, first lab);
    k := toString lab;
    if not db#?k then error("CY with label "|k|" does not exist");
    calabiYau(db#k, lab -> Q, opts)
    )

calabiYau(String, CYPolytope, Sequence) := CalabiYauInToric => opts -> (dbfilename, Q, lab) -> (
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


///
  -- h11=4 database use, 19 June 2023.
  -- XXX In construction
-*
  restart
  needsPackage "StringTorics"
*-
  R = ZZ[a,b,c,d]
  RQ = QQ (monoid R);
  (Qs, Xs) = readCYDatabase("mike-ntfe-h11-4.dbm", Ring => R);
  assert(#keys Qs == 1197) -- includes torsions and nonfavorables.
  assert(#keys Xs == 1994) -- note, none of the torsion Qs are in here yet.
  
  peek Xs#(20,0).cache
  ByH12 = partition(k -> hh^(1,2) Xs#k, keys Xs);
  -- by H12 value, 1994 examples are split into 86 groups.
  -- largest group is h12=64, at 195 in that group.
  86 == # hashTable for x in keys ByH12 list x => #ByH12#x

  -- Now let's divide by invariants to see how to separate them all.
  debug StringTorics -- invariantsAll isn't exported!
  elapsedTime IHall = partition(x -> elapsedTime invariantsAll x, values Xs);
  -- 1126 different groups here.
  assert(#keys IHall == 1126)
  (values IHall)/(x -> #x)//tally
  -- 723 different classes have exactly one element in them.
  -- largest class is 51 elements.
  -- of course, these all might be equivalent!  (Probably not, but who knows...)
  --  Tally{1 => 723}
            2 => 221
            3 => 109
            4 => 31
            5 => 8
            6 => 13
            7 => 3
            8 => 5
            9 => 1
            10 => 5
            11 => 1
            12 => 2
            13 => 2
            28 => 1
            51 => 1
  -- Now leave off inverse system invariant: get the same numbers.
  -- How many of these can be determined to be equivalent?
  -- Well, the 723 that are by themselves we can ignore.
  set2 = select(values IHall, k -> #k > 1);
  set3 = set2/(x -> (x/label//sort))
  count = 0;
  set4 = for Ls in set3 list (
      << "--- doing " << count << " with " << Ls << endl;
      count = count + 1;
      ans := elapsedTime partitionByTopology(Ls, Xs, 15);
      print ans;
      ans
      )
  set5 = for x in set4 list (
      for k in keys x list {k} | (x#k / first)
      )
  set5len = for x in set5 list (x/length)
  #set5len
  #select(set5len, x -> #x == 1) -- 340 of the 403 have one class.
  -- 53 sets have 2 classes
  --  8 sets have 3 classes
  --  2 sets have 4 classes
  -- so total number of topologies is likely:   723 + 340 + 53*2 + 8*3 + 2*4 = 1201
  set6 = for x in set5 list (x/sort/first//sort)
  set7 = sort select(set6, x -> #x > 1)
  #set7 == 63  -- these are sets we still would like to separate by invariants of that is possible
  -- range on number of topologies:
  -- low end: 723 + 340 + 63 == 1126
  --  hi end: 723 + 340 + 138 == 1201
  for ks in set7 list netList transpose {for k in ks list factor det hessian cubicForm Xs#k}
  for ks in set7 list netList transpose {for k in ks list factor cubicForm Xs#k}

  -- Here we just play some and try to separate these
  
  -- Let's try to separate some of these now, and then we can try to automate it
  -- XXX 19 June 2023.
  (L1, F1) = (c2Form Xs#(1182,0), cubicForm Xs#(1182,0))
  (L2, F2) = (c2Form Xs#(1183,2), cubicForm Xs#(1183,2))

  (L1, F1) = (c2Form Xs#(1143,0), cubicForm Xs#(1143,0))
  (L2, F2) = (c2Form Xs#(1145,0), cubicForm Xs#(1145,0))

  (L1, F1) = (c2Form Xs#(1123,0), cubicForm Xs#(1123,0))
  (L2, F2) = (c2Form Xs#(1124,0), cubicForm Xs#(1124,0))

  (L1, F1) = (c2Form Xs#(1123,0), cubicForm Xs#(1123,0))
  (L2, F2) = (c2Form Xs#(1124,0), cubicForm Xs#(1124,0))
  
  (A,phi) = genericLinearMap RQ
  T = source phi;
  I0 = sub(ideal last coefficients (phi sub(L1,T) - sub(L2,T)), ring A)
  A0 = A % I0
  phi0 = map(T,T,A0)
  trim(I0 + sub(ideal last coefficients (phi0 sub(F1,T) - sub(F2,T)), ring A))


  for k from 0 to #set7-1 list (
      (lab1,lab2) = toSequence set7#k_{0,1}; -- only do the first 2.
      I := getEquivalenceIdeal(lab1, lab2, Xs);
      << "k = " << k << " ideal " << netList I_* << endl;
      I
      )
  
  (lab1,lab2) = toSequence set7#0_{0,1}
  getEquivalenceIdeal(lab1,lab2, Xs)
  getEquivalenceIdealHelper((L1,F1),(L2,F2),A,phi)
///
    
///
  -- Example of construction of database for: h11=3, all h12's.
  -- Note, #232 is not favorable.
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  elapsedTime createCYDatabase("foo-ntfe-h11-3.dbm", topes)
  
  -- Now let's add in all the CY's total, including all triangulations.
  Qs = readCYPolytopes "foo-ntfe-h11-3.dbm";
  elapsedTime for Q in values Qs do addToCYDatabase("foo-ntfe-h11-3.dbm", Q, NTFE => true);

  -- How to access all of the polytopes and CY's at once.
  -- We create a hashtable for each, keys are their labels, and values are the CYPolytope and CalabiYauInToric's.
  Qs = readCYPolytopes "foo-ntfe-h11-3.dbm";
  for k in sort keys Qs do assert instance(Qs#k, CYPolytope)

  Xs = readCYs("foo-ntfe-h11-3.dbm", Qs);
  for k in sort keys Xs do assert instance(Xs#k, CalabiYauInToric)

  -- or both at the same time..
  (Qs1, Xs1) = readCYDatabase "foo-ntfe-h11-3.dbm";
  assert(Qs1 === Qs)
  assert(Xs1 === Xs)

  -- given all the3 Qs, how to access just some Xs (e.g. for a specific Q).
  F = openDatabase "foo-ntfe-h11-3.dbm"
    Qs = readCYPolytopes "foo-ntfe-h11-3.dbm";
    for k in sort keys F list (
      if not match("\\(6,", k) then continue else cyData(F#(toString k), i -> Qs#i)
      )
    oo/label
  close F  
///

/// 
  -- Example of creation of h11=4 database of those triangulations which are NTFE (not 2-face equivalent).
  -- One data base per (h11,h12) pair too.
-*
  restart
  needsPackage "StringTorics"
*-

  createNTFEDatabase = (dbname, topes) -> (
      createCYDatabase(dbname, topes);
      Qs := readCYPolytopes dbname;
      elapsedTime for Q in values Qs do addToCYDatabase(dbname, Q, NTFE => true);
      )
  createAllNTFEDatabases = (h12s, topes) -> (
      for i in h12s do (
          dbname = "cys-ntfe-h11-4-h12-"|i|".dbm";
          << "starting on " << dbname << endl;
          ourtopes = select(topes, ks -> last hodgeNumbers ks == i);
          createNTFEDatabase(dbname, ourtopes);
          );
      )


  topes = kreuzerSkarke(4, Limit => 10000);
  createAllNTFEDatabases({28, 34}, topes)
  assert(#topes == 1197)
  h12s = topes/hodgeNumbers/last//unique//sort
  h12s == {28, 34, 36, 37, 40, 42, 44, 46, 48, 49, 50, 52, 54, 55, 
      56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
      71, 72, 73, 74, 76, 77, 78, 79, 80, 81, 82, 84, 85, 86, 88, 
      89, 90, 91, 92, 93, 94, 96, 97, 98, 100, 101, 102, 104, 106, 
      108, 109, 110, 112, 114, 116, 118, 120, 121, 122, 124, 126, 
      128, 130, 136, 142, 144, 148, 154, 162, 166, 178, 190, 194, 
      202, 208, 214, 226, 238}

  createAllNTFEDatabases(h12s, topes)

///

///
-*
  restart
  debug needsPackage "StringTorics"
*-

  -- Example analysis of topologies for h11=4.
  topes = kreuzerSkarke(4, Limit => 10000);
  assert(#topes == 1197)
  h12s = topes/hodgeNumbers/last//unique//sort
  h12s == {28, 34, 36, 37, 40, 42, 44, 46, 48, 49, 50, 52, 54, 55, 
      56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
      71, 72, 73, 74, 76, 77, 78, 79, 80, 81, 82, 84, 85, 86, 88, 
      89, 90, 91, 92, 93, 94, 96, 97, 98, 100, 101, 102, 104, 106, 
      108, 109, 110, 112, 114, 116, 118, 120, 121, 122, 124, 126, 
      128, 130, 136, 142, 144, 148, 154, 162, 166, 178, 190, 194, 
      202, 208, 214, 226, 238}

  db4name = (h12) -> "cys-ntfe-h11-4-h12-"|h12|".dbm"
  RZ = ZZ[a,b,c,d];
  DB = hashTable for h in h12s list h => readCYDatabase(db4name h, Ring => RZ);

  (Qs, Xs) = DB#56
  H = partition(invariants2, values Xs)
  tops = hashTable for k in keys H list k => ((H#k)/label)
  INV = for k in sort keys tops list k => partitionByTopology(tops#k, Xs, 15)

  keys Xs
  F1 = cubicForm Xs#(73,0)
  F2 = cubicForm Xs#(80,3)
  L1 = c2Form Xs#(73,0)
  L2 = c2Form Xs#(80,3)
  

  (Qs, Xs) = DB#97
  (Qs, Xs) = DB#94
  H = partition(invariants2, values Xs)
  tops = hashTable for k in keys H list k => ((H#k)/label)
  INV = for k in sort keys tops list k => partitionByTopology(tops#k, Xs, 15)
  hashTable INV
  
  unfavorables = sort flatten for h in h12s list (
      Qs = readCYPolytopes(db4name h);
      for Q in values Qs list if not isFavorable Q then label Q else continue
      )
  unfavorables === {796, 800, 803, 1059, 1060, 1064, 1065, 1134, 1135, 1151, 1153, 1155}
  nontorsionfrees = sort flatten for h in h12s list (
      (Qs, Xs) = readCYDatabase(db4name h);
      for X in values Xs list (
          V := normalToricVariety(rays X, max X); 
          if classGroup V =!= ZZ^4 then label X else continue
          )
      )
  -- all unfavorables have torsion toric class group:
  nontorsionfrees == {(0, 0), (3, 0), (4, 0), (5, 0), (12, 0), (15, 0), 
      (796, 0), (796, 1), (800, 0), (803, 0), (1059, 0), (1060, 0), (1064, 0), 
      (1065, 0), (1134, 0), (1135, 0), (1151, 0), (1153, 0), (1155, 0), (1155, 1)}
  

  -- analyze one h12
  findDistincts= (h12) -> (
      (Qs, Xs) = DB#h12;
      H = partition(invariants2, values Xs);
      tops = hashTable for k in keys H list k => ((H#k)/label);
      INV = for k in sort keys tops list k => partitionByTopology(tops#k, Xs, 15);
      (INV, Qs, Xs)
      )
  (INV, Qs, Xs) = findDistincts 28;
  netList INV

  (INV, Qs, Xs) = findDistincts 34; -- only one here
  netList INV 

  (INV, Qs, Xs) = findDistincts 36; -- polytopes 3,4,5: all have class group torsion.
  netList INV 

  (INV, Qs, Xs) = findDistincts 37; 
  netList INV  -- only one

  (INV, Qs, Xs) = findDistincts 40;  -- 2 diff topologies, different invariants2
  netList INV

  (INV, Qs, Xs) = findDistincts 44;  -- has a torsion polytope.
  netList INV

  (INV, Qs, Xs) = findDistincts 46;
  netList INV -- (17,0), (19,1) are seemingly different.
  F1 = cubicForm Xs#(17,0)
  F2 = cubicForm Xs#(19,1)
  sing = (F) -> ideal F + ideal jacobian F
  linears = (I) -> (J := ideal select(I_*, f -> part(1, f) != 0); if J == 0 then trim ideal 0_(ring I) else J)
  linearcontent = (I) -> trim sum for f in (linears I)_* list content f
  betti res sub(sing F1, RQ)
  betti res sub(sing F2, RQ)

  betti res sub(saturate sing F1, RQ)
  betti res sub(saturate sing F2, RQ)
  primaryDecomposition sub(saturate sing F1, RQ)
  primaryDecomposition sub(saturate sing F2, RQ)
  see ideal gens gb saturate sing F1
  see ideal gens gb saturate sing F2 -- these have different linear parts (252a, 1512a).  SHows they are different.

  select(sort keys DB, h12 -> # keys last DB#h12 > 1)

  (INV, Qs, Xs) = findDistincts 48; 
  netList INV -- all 3 are distinct
  
  (INV, Qs, Xs) = findDistincts 49; 
  netList INV -- all 6 are distinct

  (INV, Qs, Xs) = findDistincts 50;
  netList INV -- only 1.

  (INV, Qs, Xs) = findDistincts 52; -- alot here, it seems
  netList INV -- one group of 6, one of 4, one of 3, 17 of 1 each.
  F1 = cubicForm Xs#(33,0)
  F2 = cubicForm Xs#(41,0)
  F3 = cubicForm Xs#(41,1)
  F4 = cubicForm Xs#(44,0)
  F5 = cubicForm Xs#(48,0)
  F6 = cubicForm Xs#(49,4)
  see linears ideal gens gb saturate sing F1
  see linears ideal gens gb saturate sing F2
  see linears ideal gens gb saturate sing F3
  see linears ideal gens gb saturate sing F4
  see linears ideal gens gb saturate sing F5
  see linears ideal gens gb saturate sing F6    
  trim content F1
  trim content F2  
  trim content F3
  trim content F4
  trim content F5
  trim content F6
  see linearcontent ideal gens gb saturate sing F1
  see linearcontent ideal gens gb saturate sing F2
  see linearcontent ideal gens gb saturate sing F3
  see linearcontent ideal gens gb saturate sing F4
  see linearcontent ideal gens gb saturate sing F5
  see linearcontent ideal gens gb saturate sing F6    

  F1 = cubicForm Xs#(35,0)
  F2 = cubicForm Xs#(38,0)
  F3 = cubicForm Xs#(39,0) -- F1 and F3 look pretty similar? TODO: I can't prove they are the same or different yet!
  F4 = cubicForm Xs#(53,0)

  for p in {3,5,7,11, 13} list (pointCount(F1, p), pointCount(F3, p))
  decompose sub(sing F1, RQ)
  decompose sub(ideal gens gb saturate sing F1, RQ)
  decompose sub(ideal gens gb saturate sing F3, RQ)

  invariants2 Xs#(35,0)
  invariants2 Xs#(39,0)
  partitionGVConeByGV(Xs#(35, 0), DegreeLimit => 25) -- simplicial (4 generators, GV's: -2, 8, 10, 64).
  partitionGVConeByGV(Xs#(39, 0), DegreeLimit => 25) -- 5 gens, GV: -2,-2,10,64,128.
  -- Question: how to distinguish F1, F3?  (F2, F4 are different and diff from F1, F3).

  L1 = c2Form Xs#(35,0)
  L3 = c2Form Xs#(39,0)
 
  decompose sub((ideal(L1) + sing F1), RQ)
  decompose sub((ideal(L3) + sing F3), RQ)
  decompose sub(ideal(L3, F3), RQ)

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 54;
  netList INV 
  F1 = cubicForm Xs#(64,0)
  F2 = cubicForm Xs#(64,1)
  see linears ideal gens gb saturate sing F1
  see linears ideal gens gb saturate sing F2 -- very different...

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 55;
  netList INV 
  F1 = cubicForm Xs#(67,0)
  F2 = cubicForm Xs#(69,0)
  see linears ideal gens gb saturate sing F1
  see linears ideal gens gb saturate sing F2 -- different linear part contents.

  -- Next one to try.  All ones that cannot be matched are distinct.
  (INV, Qs, Xs) = findDistincts 56;
  netList INV -- lots.  2: 2 diff, 1: 3 diff, 18: 1 only...
  -- set1
  F1 = cubicForm Xs#(73,0)
  F2 = cubicForm Xs#(80,3)
  see linears ideal gens gb saturate sing F1
  see linears ideal gens gb saturate sing F2 -- different linear part contents.
  -- set2
  F1 = cubicForm Xs#(75,0)
  F2 = cubicForm Xs#(75,1)
  F3 = cubicForm Xs#(75,3)
  see linears ideal gens gb saturate sing F1
  see linears ideal gens gb saturate sing F2 
  see linears ideal gens gb saturate sing F3 -- all different via linear part contents
  -- set3
  F1 = cubicForm Xs#(72,0)
  F2 = cubicForm Xs#(80,2)
  see linears ideal gens gb saturate sing F1
  see linears ideal gens gb saturate sing F2 -- different linear part contents.

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 57; -- 1 only.
  netList INV 

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 58; 
  netList INV -- 1: 9 different!, 1: 4 diff, 26: 1 diff
  -- set1
  F1 = cubicForm Xs#(92,0)
  F2 = cubicForm Xs#(93,0)
  F3 = cubicForm Xs#(96,0)
  F4 = cubicForm Xs#(97,0)
  F5 = cubicForm Xs#(105,0)
  F6 = cubicForm Xs#(105,2)
  F7 = cubicForm Xs#(112,0)
  F8 = cubicForm Xs#(120,0)
  F9 = cubicForm Xs#(120,2)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2
  linearcontent ideal gens gb saturate sing F3
  linearcontent ideal gens gb saturate sing F4
  linearcontent ideal gens gb saturate sing F5
  linearcontent ideal gens gb saturate sing F6
  linearcontent ideal gens gb saturate sing F7
  linearcontent ideal gens gb saturate sing F8
  linearcontent ideal gens gb saturate sing F9
  -- all are distinct except possibly F3, F7
  decompose sub(ideal gens gb saturate sing F3, RQ)
  decompose sub(ideal gens gb saturate sing F7, RQ) -- 1 point vs 3 points.  So F3, F7, therefore all, are distinct.
  -- set2
  F1 = cubicForm Xs#(98,0)
  F2 = cubicForm Xs#(109,0)
  F3 = cubicForm Xs#(110,0)
  F4 = cubicForm Xs#(110,1)
  see linears ideal gens gb saturate sing F1
  see linears ideal gens gb saturate sing F2 
  see linears ideal gens gb saturate sing F3
  see linears ideal gens gb saturate sing F4 -- all different via linear part contents

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 59; 
  netList INV -- 3, all distinct

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 60; 
  netList INV -- 
  -- set1
  F1 = cubicForm Xs#(133,0)
  F2 = cubicForm Xs#(137,0)
  F3 = cubicForm Xs#(140,0)
  F4 = cubicForm Xs#(146,0)
  F5 = cubicForm Xs#(148,0)
  F6 = cubicForm Xs#(148,1)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2
  linearcontent ideal gens gb saturate sing F3
  linearcontent ideal gens gb saturate sing F4
  linearcontent ideal gens gb saturate sing F5
  linearcontent ideal gens gb saturate sing F6 -- all distinct linear content!
  -- set2
  F1 = cubicForm Xs#(144,0)
  F2 = cubicForm Xs#(158,0)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2 -- all distinct linear content!
  -- set3
  F1 = cubicForm Xs#(140,1)
  F2 = cubicForm Xs#(142,0)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2 -- all distinct linear content!
  -- set4
  F1 = cubicForm Xs#(129,0)
  F2 = cubicForm Xs#(135,0)
  F3 = cubicForm Xs#(139,0)
  F4 = cubicForm Xs#(149,0)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2 
  linearcontent ideal gens gb saturate sing F3 
  linearcontent ideal gens gb saturate sing F4 -- all distinct linear content!
  -- set5
  F1 = cubicForm Xs#(135,1)
  F2 = cubicForm Xs#(138,0)
  F3 = cubicForm Xs#(149,2)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2 -- how to tell F1, F2 apart?
  linearcontent ideal gens gb saturate sing F3 
  decompose sub(ideal gens gb saturate sing F1, RQ)
  decompose sub(ideal gens gb saturate sing F2, RQ) -- stil todo: are F1, F2 equivalent?

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 61;
  netList INV -- 
  F1 = cubicForm Xs#(163,2)
  F2 = cubicForm Xs#(164,0)
  F3 = cubicForm Xs#(166,0)
  F4 = cubicForm Xs#(166,1)
  F5 = cubicForm Xs#(166,2)
  F6 = cubicForm Xs#(169,0)
  F7 = cubicForm Xs#(169,1)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2
  linearcontent ideal gens gb saturate sing F3
  linearcontent ideal gens gb saturate sing F4
  linearcontent ideal gens gb saturate sing F5
  linearcontent ideal gens gb saturate sing F6
  linearcontent ideal gens gb saturate sing F7 -- all distinct linear content!

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 62;
  netList INV -- 
  -- set1
  F1 = cubicForm Xs#(172,0)
  F2 = cubicForm Xs#(173,0)
  F3 = cubicForm Xs#(179,0)
  F4 = cubicForm Xs#(179,2)
  F5 = cubicForm Xs#(180,0)
  F6 = cubicForm Xs#(192,1)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2
  linearcontent ideal gens gb saturate sing F3
  linearcontent ideal gens gb saturate sing F4
  linearcontent ideal gens gb saturate sing F5
  linearcontent ideal gens gb saturate sing F6  -- all distinct linear content!
  -- set2
  F1 = cubicForm Xs#(181,1)
  F2 = cubicForm Xs#(191,0)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2 -- how to tell F1, F2 apart?
  -- set3
  F1 = cubicForm Xs#(182,0)
  F2 = cubicForm Xs#(182,1)
  linearcontent ideal gens gb saturate sing F1
  linearcontent ideal gens gb saturate sing F2 -- how to tell F1, F2 apart?

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 63;
  netList INV -- only 1

  -- Next one to try
  (INV, Qs, Xs) = findDistincts 64;
  netList INV -- LOTS!!  TODO: start here.  Actually, maybe I should start over and incorporate the linear content into the invariants.

  

  H = partition(invariants2, values Xs)
  tops = hashTable for k in keys H list k => ((H#k)/label)
  INV = for k in sort keys tops list k => partitionByTopology(tops#k, Xs, 15)
  hashTable INV

  h12s == {28, 34, 36, 37, 40, 42, 44, 46, 48, 49, 50, 52, 54, 55, 
      56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
      71, 72, 73, 74, 76, 77, 78, 79, 80, 81, 82, 84, 85, 86, 88, 
      89, 90, 91, 92, 93, 94, 96, 97, 98, 100, 101, 102, 104, 106, 
      108, 109, 110, 112, 114, 116, 118, 120, 121, 122, 124, 126, 
      128, 130, 136, 142, 144, 148, 154, 162, 166, 178, 190, 194, 
      202, 208, 214, 226, 238}

///

///
  -- An attempt to automate the following:
  -- 1. create one hash table with all of the Qs, Xs.
  -- 2. separate the Xs via invariants.
  -- 3. for each invariant, use topology to separate.
  
-*
  restart
  debug needsPackage "StringTorics"
*-

  -- Example analysis of topologies for h11=4.
  --topes = kreuzerSkarke(4, Limit => 10000);
  --assert(#topes == 1197)
  --h12s = topes/hodgeNumbers/last//unique//sort
  h12s = {28, 34, 36, 37, 40, 42, 44, 46, 48, 49, 50, 52, 54, 55, 
      56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
      71, 72, 73, 74, 76, 77, 78, 79, 80, 81, 82, 84, 85, 86, 88, 
      89, 90, 91, 92, 93, 94, 96, 97, 98, 100, 101, 102, 104, 106, 
      108, 109, 110, 112, 114, 116, 118, 120, 121, 122, 124, 126, 
      128, 130, 136, 142, 144, 148, 154, 162, 166, 178, 190, 194, 
      202, 208, 214, 226, 238}

  nontorsionfrees = {(0, 0), (3, 0), (4, 0), (5, 0), (12, 0), (15, 0), 
      (796, 0), (796, 1), (800, 0), (803, 0), (1059, 0), (1060, 0), (1064, 0), 
      (1065, 0), (1134, 0), (1135, 0), (1151, 0), (1153, 0), (1155, 0), (1155, 1)}

  db4name = (h12) -> "cys-ntfe-h11-4-h12-"|h12|".dbm"
  RZ = ZZ[a,b,c,d];
  RQ = QQ (monoid RZ);
  DB = hashTable for h in h12s list h => readCYDatabase("../m2-examples/"|db4name h, Ring => RZ);
  Qs = new MutableHashTable
  Xs = new MutableHashTable
  for a in keys DB do (
      (Q1s, X1s) = DB#a;
      for q in sort keys Q1s do Qs#q = Q1s#q;
      for k in sort keys X1s do if not member(k, nontorsionfrees) then Xs#k = X1s#k;
      );
  Xs = new HashTable from Xs;
  Qs = new HashTable from Qs;
  sort keys Qs
  sort keys Xs
  #oo == 1994 -- 2014 total for all
  assert(# sort keys Xs == 2014 - #nontorsionfrees)
  
  elapsedTime H = partition(invariants3, values Xs);
  H1 = hashTable for k in keys H list k => (H#k)/label;
  elapsedTime INV = for k in sort keys H1 list k => partitionByTopology(H1#k, Xs, 15);
  -- number of different topologies
  tally for k in sort keys H1 list #H1#k
  elapsedTime INV = for k in sort keys H1 list k => partitionByTopology(H1#k, Xs, 15);
  tally for k in INV list #(keys last k)
  INV1 = select(INV, k -> #(keys last k) >= 2);
  netList INV1
  
  -- separate topologies:
  -- output is a list of lists {S0, S1, S2, ..., Sr}
  -- each CY in Si is not equivalent to any in Sj, i!=j.
  -- each pair of CY's in a specific Si are not proved to be distinct (but are likely distinct?)
  -- every CY in h11=4 database is knoqn to be equivalent to one of these 
  --  (Question: or are they all there in one of these lists?)
  TOP1 = for k in INV list for a in keys last k list (
      prepend(a, ((last k)#a)/first)
      );
  netList(TOP1/(a -> a/(a1->#a1)))

  REPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then continue
      else K)

  ONEONLY = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then first K else continue)

  -- So: 1994 total different NTFE h11=4 CY's, with torsion free class groups (all of these are favorable)
  -- Of these:
  #TOP1
  TOP1/(a -> a/(a1->#a1)//sum)//sum == 1994
  #ONEONLY == 904 -- the number that are by themselves.
  #REPS == 118 
  # flatten REPS == 287

  REPS = {{(14, 1), (14, 0)}, 
      {(38, 4), (53, 0), (30, 0), (38, 0)}, 
      {(35, 0), (39, 0)}, 
      {(80, 6), (80, 2)}, 
      {(96, 0), (126, 0)}, 
      {(148, 1), (141, 0)}, 
      {(140, 0), (144, 0)}, 
      {(135, 0), (134, 3)}, 
      {(135, 1), (138, 0)}, 
      {(180, 0), (179, 2)}, 
      {(241, 0), (210, 0), (258, 0)}, 
      {(265, 0), (250, 0), (265, 1), (236, 0), (228, 1)}, 
      {(290, 1), (258, 1)}, 
      {(235, 0), (262, 0)}, 
      {(225, 2), (260, 2)}, 
      {(265, 2), (235, 1)}, 
      {(265, 3), (228, 0)}, 
      {(250, 1), (254, 2)}, 
      {(249, 1), (246, 0)}, 
      {(227, 0), (261, 0)}, 
      {(209, 3), (284, 1), (261, 1)}, 
      {(240, 3), (278, 1)}, 
      {(319, 0), (313, 0)}, 
      {(316, 0), (322, 0)}, 
      {(334, 2), (328, 1), (329, 0)}, 
      {(331, 0), (337, 1), (339, 3)}, 
      {(334, 1), (329, 1)}, 
      {(339, 5), (337, 0)}, 
      {(385, 2), (385, 3)}, 
      {(364, 1), (377, 0)}, 
      {(356, 0), (395, 6)}, 
      {(378, 0), (378, 2), (379, 1), (395, 2)}, 
      {(379, 0), (379, 3)}, 
      {(383, 0), (356, 1)}, 
      {(397, 0), (350, 0)}, 
      {(387, 4), (387, 0)}, 
      {(377, 2), (364, 0)}, 
      {(416, 0), (441, 0), (438, 0), (450, 4)}, 
      {(400, 0), (403, 0), (405, 0), (414, 0)}, 
      {(436, 0), (433, 0)}, 
      {(419, 0), (419, 2)}, 
      {(436, 1), (433, 1)}, 
      {(451, 0), (455, 1)}, 
      {(475, 1), (481, 1)}, 
      {(473, 0), (478, 0)}, 
      {(458, 5), (458, 1)}, 
      {(458, 3), (458, 0)}, 
      {(508, 0), (514, 0)}, 
      {(512, 1), (505, 0)}, 
      {(507, 0), (519, 3)}, 
      {(526, 1), (527, 0), (518, 3), (522, 0)}, 
      {(524, 0), (518, 0)}, 
      {(529, 0), (510, 0)}, 
      {(550, 0), (538, 0)}, 
      {(532, 0), (573, 0)}, 
      {(556, 0), (562, 0)}, 
      {(553, 0), (548, 1), (554, 0)}, 
      {(552, 0), (551, 2)}, 
      {(559, 0), (577, 0)}, 
      {(540, 0), (545, 0)}, 
      {(588, 0), (586, 0)}, 
      {(616, 0), (606, 0)}, 
      {(650, 0), (650, 1), (667, 1), (650, 3), (656, 0), (666, 9), (630, 0)}, 
      {(643, 0), (649, 0), (626, 0)}, 
      {(628, 0), (653, 0)}, 
      {(647, 0), (669, 0)}, 
      {(688, 1), (685, 0)}, 
      {(687, 0), (684, 1)}, 
      {(693, 0), (694, 0)}, 
      {(712, 0), (714, 0), (706, 0), (707, 0), (709, 0)}, 
      {(716, 3), (705, 0)}, 
      {(715, 0), (704, 0)}, 
      {(738, 3), (737, 1)}, 
      {(731, 0), (740, 1)}, 
      {(764, 0), (764, 1)}, 
      {(773, 7), (773, 0)}, 
      {(770, 1), (774, 5)}, 
      {(771, 0), (771, 1)}, 
      {(781, 2), (781, 1)}, 
      {(815, 0), (810, 0), (851, 0), (884, 0), (820, 0), (869, 0)}, 
      {(884, 7), (881, 0), (884, 2)}, 
      {(835, 2), (806, 0), (822, 0)}, 
      {(878, 0), (885, 1), (886, 0)}, 
      {(879, 9), (808, 0), (845, 0), (896, 0), (865, 0), (879, 3), (896, 3)}, 
      {(825, 0), (862, 0)}, 
      {(863, 0), (887, 0)}, 
      {(876, 0), (854, 0)}, 
      {(900, 0), (901, 0)}, 
      {(913, 2), (913, 1)}, 
      {(935, 17), (938, 0)}, 
      {(938, 3), (933, 0)}, 
      {(919, 0), (938, 4)}, 
      {(930, 1), (934, 0), (922, 0)}, 
      {(932, 0), (916, 0), (917, 0)}, 
      {(922, 1), (930, 0)}, 
      {(938, 11), (933, 5)}, 
      {(935, 0), (937, 0), (915, 0), (924, 0), (935, 5)}, 
      {(943, 0), (953, 0), (958, 0)}, 
      {(960, 0), (962, 0)}, 
      {(1002, 0), (987, 0), (988, 0), (996, 1)}, 
      {(987, 5), (1004, 0), (994, 0)}, 
      {(976, 0), (989, 0)}, 
      {(993, 0), (997, 0)}, 
      {(983, 0), (998, 2)}, 
      {(979, 0), (1005, 3), (1002, 2)}, 
      {(1000, 0), (986, 0)}, 
      {(1023, 0), (1027, 2)}, 
      {(1039, 0), (1028, 0)}, 
      {(1051, 0), (1048, 5)}, 
      {(1078, 1), (1077, 3), (1082, 3)}, 
      {(1067, 0), (1068, 0)}, 
      {(1085, 0), (1082, 0), (1086, 0)}, 
      {(1090, 2), (1094, 4)}, 
      {(1090, 0), (1094, 0)}, 
      {(1121, 0), (1122, 0)}, 
      {(1123, 0), (1124, 0)}, 
      {(1148, 0), (1149, 0)}, 
      {(1183, 2), (1182, 0)}}
  REPS = REPS/sort

  for S in REPS list (#S - # (S/(lab -> invariants4 Xs#lab)//unique))
  for S in REPS list {#S,  # (S/(lab -> invariants4 Xs#lab)//unique)}
  
  #REPS == 118 -- this is the number of sets that we don't know yet how to separate
               -- without use of GV's (which is heuristic).

  sort REPS#0
  X1 = Xs#(14,0)
  X2 = Xs#(14,1)
  invariants4 X1
  invariants4 X2
  F1 = cubicForm X1
  F2 = cubicForm X2
  sing = (F) -> trim saturate(F + ideal jacobian F)
  decompose(sub(sing F1, RQ))
  sing F2
  
///

///
  -- 2 Jan 2023.
  -- An attempt to automate the following:
  -- 1. create one hash table with all of the Qs, Xs.
  -- 2. separate the Xs via invariants4 (all the ones except point counts)
  -- 3. for each invariant, use topology to separate.
-*
  restart
  debug needsPackage "StringTorics"
*-

  -- Example analysis of topologies for h11=4.
  --topes = kreuzerSkarke(4, Limit => 10000);
  --assert(#topes == 1197)
  --h12s = topes/hodgeNumbers/last//unique//sort
  h12s = {28, 34, 36, 37, 40, 42, 44, 46, 48, 49, 50, 52, 54, 55, 
      56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
      71, 72, 73, 74, 76, 77, 78, 79, 80, 81, 82, 84, 85, 86, 88, 
      89, 90, 91, 92, 93, 94, 96, 97, 98, 100, 101, 102, 104, 106, 
      108, 109, 110, 112, 114, 116, 118, 120, 121, 122, 124, 126, 
      128, 130, 136, 142, 144, 148, 154, 162, 166, 178, 190, 194, 
      202, 208, 214, 226, 238}

  nontorsionfrees = {(0, 0), (3, 0), (4, 0), (5, 0), (12, 0), (15, 0), 
      (796, 0), (796, 1), (800, 0), (803, 0), (1059, 0), (1060, 0), (1064, 0), 
      (1065, 0), (1134, 0), (1135, 0), (1151, 0), (1153, 0), (1155, 0), (1155, 1)}

  db4name = (h12) -> "cys-ntfe-h11-4-h12-"|h12|".dbm"
  RZ = ZZ[a,b,c,d];
  RQ = QQ (monoid RZ);
  DB = hashTable for h in h12s list h => readCYDatabase("../m2-examples/"|db4name h, Ring => RZ);
  Qs = new MutableHashTable
  Xs = new MutableHashTable
  for a in keys DB do (
      (Q1s, X1s) = DB#a;
      for q in sort keys Q1s do Qs#q = Q1s#q;
      for k in sort keys X1s do if not member(k, nontorsionfrees) then Xs#k = X1s#k;
      );
  Xs = new HashTable from Xs;
  Qs = new HashTable from Qs;
  #(keys Xs) == 1994 -- 2014 total for all, including torsions.
  assert(# sort keys Xs == 2014 - #nontorsionfrees)
  
  elapsedTime H = partition(invariants4, values Xs); --  135 sec
  H1 = hashTable for k in keys H list k => (H#k)/label;
  elapsedTime INV = for k in sort keys H1 list k => partitionByTopology(H1#k, Xs, 15); -- 751 sec
  -- number of different topologies

-----------------------------------
  -- new code: let's first see which ones are identical.
  -- XXXX
  # keys Xs == 1994
  Diff1 = partition(lab -> (X := Xs#lab; {c2Form X, cubicForm X}), sort keys Xs);
  select(pairs Diff1, kv -> # kv#1 > 1)

  elapsedTime H = partition(lab -> (
          X := Xs#lab;
          elapsedTime join({hh^(1,1) X, hh^(1,2) X}, invariantsAll Xs#lab)
          ), sort keys Xs); --  295 sec  
  #(keys H) == 1294 -- this leaves out h11, h12.
  #(keys H) == 1300 -- this has h11, h12.  Interesting.

  elapsedTime INV = for k in sort keys H list k => partitionByTopology(H#k, Xs, 15); -- 751 sec
  -- number of different topologies: is at least 1300, at most XXX

  elapsedTime H4 = partition(lab -> elapsedTime invariants4 Xs#lab, sort keys  Xs); 
  #(keys H4) == 1040

  REPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then continue
      else K)

  # select(INV, k -> #k#1 == 1) == 1248
  # select(INV, k -> #k#1 == 2) == 47
  # select(INV, k -> #k#1 == 3) == 5
  # select(INV, k -> #k#1 > 3) == 0

  #REPS==52
  -- # different topologies is bounded above by 1248 + 2*47 + 3*5 == 1357
  --                           bounded below by 1248 + 52 = 1300
  
  REPS#0
  hashTable invariantsAll Xs#(254,5)
  hashTable invariantsAll Xs#(228,1)
  hashTable invariantsAll Xs#(REPS#2#0)
----------------------------------

  
  -- separate topologies:
  -- output is a list of lists {S0, S1, S2, ..., Sr}
  -- each CY in Si is not equivalent to any in Sj, i!=j.
  -- each pair of CY's in a specific Si are not proved to be distinct (but are likely distinct?)
  -- every CY in h11=4 database is knoqn to be equivalent to one of these 
  --  (Question: or are they all there in one of these lists?)
  TOP1 = for k in INV list for a in keys last k list (
      prepend(a, ((last k)#a)/first)
      );
  netList(TOP1/(a -> a/(a1->#a1)))

  REPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then continue
      else K)

  ONEONLY = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then first K else continue)
  
  #INV == 1040
  #ONEONLY == 932
  #REPS = 108
  
  netList REPS
  REPS = {
    {(14, 1), (14, 0)}, 
    {(53, 0), (38, 0)}, 
    {(35, 0), (39, 0)}, 
    {(30, 0), (38, 4)}, 
    {(80, 6), (80, 2)}, 
    {(96, 0), (126, 0)}, 
    {(135, 1), (138, 0)}, 
    {(180, 0), (179, 2)}, 
    {(265, 0), (250, 0), (265, 1), (236, 0), (228, 1)}, 
    {(235, 0), (262, 0)}, 
    {(225, 2), (260, 2)}, 
    {(265, 2), (235, 1)}, 
    {(265, 3), (228, 0)}, 
    {(250, 1), (254, 2)}, 
    {(249, 1), (246, 0)}, 
    {(210, 0), (258, 0)}, 
    {(209, 3), (284, 1), (261, 1)}, 
    {(240, 3), (278, 1)}, 
    {(319, 0), (313, 0)}, 
    {(316, 0), (322, 0)}, 
    {(334, 2), (328, 1), (329, 0)}, 
    {(331, 0), (337, 1), (339, 3)}, 
    {(334, 1), (329, 1)}, 
    {(339, 5), (337, 0)}, 
    {(385, 2), (385, 3)}, 
    {(364, 1), (377, 0)}, 
    {(356, 0), (395, 6)}, 
    {(379, 1), (395, 2), (378, 0)}, 
    {(383, 0), (356, 1)}, 
    {(397, 0), (350, 0)}, 
    {(387, 4), (387, 0)}, 
    {(377, 2), (364, 0)}, 
    {(416, 0), (441, 0), (438, 0)}, 
    {(400, 0), (403, 0), (405, 0), (414, 0)}, 
    {(436, 0), (433, 0)}, 
    {(419, 0), (419, 2)}, 
    {(436, 1), (433, 1)}, 
    {(451, 0), (455, 1)}, 
    {(473, 0), (478, 0)}, 
    {(458, 5), (458, 1)}, 
    {(458, 3), (458, 0)}, 
    {(508, 0), (514, 0)}, 
    {(512, 1), (505, 0)}, 
    {(507, 0), (519, 3)}, 
    {(526, 1), (527, 0), (518, 3), (522, 0)}, 
    {(524, 0), (518, 0)}, 
    {(529, 0), (510, 0)}, 
    {(550, 0), (538, 0)}, 
    {(556, 0), (562, 0)}, 
    {(553, 0), (554, 0)}, 
    {(552, 0), (551, 2)}, 
    {(559, 0), (577, 0)}, 
    {(540, 0), (545, 0)}, 
    {(588, 0), (586, 0)}, 
    {(616, 0), (606, 0)}, 
    {(650, 0), (650, 1), (667, 1), (650, 3), (656, 0), (666, 9), (630, 0)}, 
    {(643, 0), (649, 0), (626, 0)}, 
    {(628, 0), (653, 0)}, 
    {(647, 0), (669, 0)}, 
    {(688, 1), (685, 0)}, 
    {(687, 0), (684, 1)}, 
    {(693, 0), (694, 0)}, 
    {(712, 0), (714, 0), (706, 0), (707, 0), (709, 0)}, 
    {(716, 3), (705, 0)}, 
    {(715, 0), (704, 0)}, 
    {(738, 3), (737, 1)}, 
    {(731, 0), (740, 1)}, 
    {(773, 7), (773, 0)}, 
    {(770, 1), (774, 5)}, 
    {(771, 0), (771, 1)}, 
    {(781, 2), (781, 1)}, 
    {(815, 0), (810, 0), (851, 0), (884, 0), (820, 0), (869, 0)}, 
    {(884, 7), (881, 0), (884, 2)}, 
    {(835, 2), (806, 0), (822, 0)}, 
    {(878, 0), (885, 1), (886, 0)}, 
    {(879, 9), (808, 0), (845, 0), (896, 0), (865, 0), (879, 3), (896, 3)}, 
    {(825, 0), (862, 0)}, 
    {(863, 0), (887, 0)}, 
    {(876, 0), (854, 0)}, 
    {(900, 0), (901, 0)}, 
    {(935, 17), (938, 0)}, 
    {(938, 3), (933, 0)}, 
    {(930, 1), (934, 0), (922, 4)}, 
    {(932, 0), (916, 0), (917, 0)}, 
    {(922, 1), (930, 0)}, 
    {(938, 11), (933, 5)}, 
    {(935, 0), (937, 0), (915, 0), (924, 0), (935, 5)}, 
    {(943, 0), (953, 0), (958, 0)}, 
    {(960, 0), (962, 0)}, 
    {(987, 0), (988, 0), (1002, 0)}, 
    {(987, 5), (1004, 0), (994, 0)}, 
    {(976, 0), (989, 0)}, 
    {(993, 0), (997, 0)}, 
    {(983, 0), (998, 2)}, 
    {(979, 0), (1005, 3), (1002, 2)}, 
    {(1000, 0), (986, 0)}, 
    {(1023, 0), (1027, 2)}, 
    {(1039, 0), (1028, 0)}, 
    {(1051, 0), (1048, 5)}, 
    {(1078, 1), (1077, 3), (1082, 3)}, 
    {(1067, 0), (1068, 0)}, 
    {(1085, 0), (1082, 0), (1086, 0)}, 
    {(1090, 2), (1094, 4)}, 
    {(1090, 0), (1094, 0)}, 
    {(1121, 0), (1122, 0)}, 
    {(1123, 0), (1124, 0)}, 
    {(1148, 0), (1149, 0)}, 
    {(1183, 2), (1182, 0)}}  
  REPS = REPS/sort;
  netList REPS
  -- The ones in REPS are the ones we need to separate (if they are truly different!).

  sing = F -> trim(ideal F + ideal jacobian F)
  toQQ = I -> sub(I, RQ)

  --- REPS#0
  LFs = for lab in REPS#0 list (X := Xs#lab; {toQQ c2Form X, toQQ cubicForm X})
  netList LFs
  
  F0 = LFs#0#1
  F1 = LFs#1#1
  for a in {-2,-2,-2,-2}..{2,2,2,2} list if sub(F0, matrix{a}) == 0 then a else continue
  for a in {-2,-2,-2,-2}..{2,2,2,2} list if sub(F1, matrix{a}) == 0 then a else continue
  saturate sing F0
  saturate sing F1
  factor(F0 - F1)

  findLinearMaps(List, List) := List => (LF1, LF2) -> (
      -- not complete...
      (L1,F1) := toSequence LF1;
      (L2,F2) := toSequence LF2;
      RQ :=ring L1;
      n := numgens RQ;
      t := symbol t;
      T := QQ[t_(1,1)..t_(n,n)];
      TR := T (monoid RQ);
      M := genericMatrix(T, n, n);
      phi := map(TR, TR, M)
      )
  
  phi = findLinearMaps(LFs#0, LFs#1)
  T = target phi
  TC = coefficientRing T
  phi

  (L1,F1) = toSequence (LFs#0/(g -> sub(g, T)))
  (L2,F2) = toSequence (LFs#1/(g -> sub(g, T)))
  L1 = L1/2
  L2 = L2/2
      -- now we make the ideals for each key, and each permutation.
      I1 = trim sub(ideal last coefficients (phi L1 - L2), TC)
      I2 = trim sub(ideal last coefficients (phi F1 - F2), TC)
      I = trim(I1 + I2);
      Tp = ZZ/101[t_(1,1), t_(1,2), t_(1,3), t_(1,4), t_(2,1), t_(2,2), t_(2,3), t_(2,4), t_(3,1), t_(3,2), t_(3,3), t_(3,4), t_(4,1), t_(4,2), t_(4,3), t_(4,4)]
      psi = map(Tp, TC, gens Tp)
      J = psi I
      gbTrace=3
      gens gb(J, DegreeLimit => 10);
      ids := for k in keys gv1 list (
        perms := permutations(#gv1#k);
        mat1 := transpose matrix gv1#k;
        mat2 := transpose matrix gv2#k;
        for p in perms list (
            I := trim ideal (M * mat1 - mat2_p); 
            if I == 1 then continue else I
            )
        );
    topval := ids/(x -> #x - 1);
    zeroval := ids/(x -> 0);
    fullIdeals := for a in zeroval .. topval list (
        J := trim sum for i from 0 to #ids-1 list ids#i#(a#i);
        if J == 1 then continue else J
        );
    Ms := for i in fullIdeals list M % i;
    --newMs := select(Ms, m -> (d := det m; d == 1 or d == -1));
    --if any(newMs, m -> support m =!= {}) then << "some M is not reduced to a constant" << endl;
    Ms
    )

  hessian = (F) -> diff(vars ring F, diff(transpose vars ring F, F))

  positions(REPS, x -> #x == 5) == {8, 62, 86}

  LFs = for lab in REPS#0 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1)
  netList for a in LFs list (betti res sing toQQ a#1) 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1)

  LFs = for lab in REPS#1 list (X := Xs#lab; {c2Form X, cubicForm X});
  for a in LFs list (decompose sing toQQ a#1)
  for a in LFs list (betti res sing toQQ a#1) -- different.

  LFs = for lab in REPS#2 list (X := Xs#lab; {c2Form X, cubicForm X});
  for a in LFs list (decompose sing toQQ a#1)
  for a in LFs list (betti res sing toQQ a#1) 
  for a in LFs list (saturate sing2 ideal a) -- different

  LFs = for lab in REPS#3 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1)
  netList for a in LFs list (betti res sing toQQ a#1) 
  netList for a in LFs list (saturate sing2 ideal a) -- could be the same
  netList for a in LFs list (factor det hessian toQQ a#1) -- contents of hessians are distinct, so I think these must be distinct.

  LFs = for lab in REPS#4 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1)
  netList for a in LFs list (betti res sing toQQ a#1) 
  netList for a in LFs list (saturate sing2 ideal a) 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same.

  LFs = for lab in REPS#5 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) -- actually all 4 invariants are different!
  netList for a in LFs list (saturate sing2 ideal a) 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#6 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) 
  netList for a in LFs list (factor det hessian toQQ a#1) -- seem to have different contents, so different...

  LFs = for lab in REPS#7 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#8 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,2,4 might be the same? 1,3 are diff.

  LFs = for lab in REPS#9 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) --

  LFs = for lab in REPS#10 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) --

  LFs = for lab in REPS#11 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) --

  LFs = for lab in REPS#12 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- hessians have different content

  LFs = for lab in REPS#13 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same.

  LFs = for lab in REPS#14 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same.

  LFs = for lab in REPS#15 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- different content

  LFs = for lab in REPS#16 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- (0,2), 1
  netList for a in LFs list (betti res sing toQQ a#1) -- 0,1,2 all different
  netList for a in LFs list (saturate sing2 ideal a) -- all different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#17 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different content
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#18 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#19 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#20 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 0, (1,2)
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 1,2 could be same

  LFs = for lab in REPS#21 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1,2 could all be same

  LFs = for lab in REPS#22 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be the same

  LFs = for lab in REPS#23 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be the same

  LFs = for lab in REPS#24 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- different

  LFs = for lab in REPS#25 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same.

  LFs = for lab in REPS#26 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- different
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- different

  LFs = for lab in REPS#27 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- (0,1), 2
  netList for a in LFs list (saturate sing2 ideal a) -- all different
  netList for a in LFs list (factor det hessian toQQ a#1) --

  LFs = for lab in REPS#28 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same (can check)

  LFs = for lab in REPS#29 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#30 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#31 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#32 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- (0,2), 1
  netList for a in LFs list (saturate sing2 ideal a) -- all different
  netList for a in LFs list (factor det hessian toQQ a#1) -- all different

  LFs = for lab in REPS#33 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- (0,1,2), 3
  netList for a in LFs list (saturate sing2 ideal a) -- (0,1), 2, 3
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be the same.

  LFs = for lab in REPS#34 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same.

  LFs = for lab in REPS#35 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same.

  LFs = for lab in REPS#36 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same.

  LFs = for lab in REPS#37 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same.

  LFs = for lab in REPS#38 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same.

  LFs = for lab in REPS#39 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- different

  LFs = for lab in REPS#40 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- different
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#41 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) -- different
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#42 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) -- different
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#43 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) -- different
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#44 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- (0,1,2), 3
  netList for a in LFs list (saturate sing2 ideal a) -- all 4 different.
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#45 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- different contents

  LFs = for lab in REPS#46 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#47 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- different contents

  LFs = for lab in REPS#48 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- different 
  netList for a in LFs list (factor det hessian toQQ a#1) -- different contents

  LFs = for lab in REPS#49 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) --
  netList for a in LFs list (factor det hessian toQQ a#1) -- different contents

  LFs = for lab in REPS#50 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) --
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#51 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) --
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#52 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) --
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#53 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- different contents

  LFs = for lab in REPS#54 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#55 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- (0,1,3,4), 2, (5,6)
  netList for a in LFs list (betti res sing toQQ a#1) -- (0,3,4), 1, 2, 5, 6
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,3,4 could be same

  LFs = for lab in REPS#56 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1,2 could all be same

  LFs = for lab in REPS#57 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could all be same

  LFs = for lab in REPS#58 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could all be same

  LFs = for lab in REPS#59 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#60 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) --different
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#61 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#62 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- 0, (1,2,3,4)
  netList for a in LFs list (factor det hessian toQQ a#1) -- (1,2,3,4) could be same

  LFs = for lab in REPS#63 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 could be same

  LFs = for lab in REPS#64 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- different contents

  LFs = for lab in REPS#65 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#66 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) --
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#67 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different 
  netList for a in LFs list (betti res sing toQQ a#1) -- dofferent
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#68 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#69 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#70 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#71 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- (0,1,3), (2,5), 4
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1,3) could be same, (2,5) could be same.

  LFs = for lab in REPS#72 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- all 3 different
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- all different
  netList for a in LFs list (factor det hessian toQQ a#1) --

  LFs = for lab in REPS#73 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- (0,2), 1
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,2) could be same.

  LFs = for lab in REPS#74 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- all 3 different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#75 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- (0,4,6), (1,2,3), 5.
  netList for a in LFs list (betti res sing toQQ a#1) -- (0,4,6), 5, (1,2), 3, 5.
  netList for a in LFs list (saturate sing2 ideal a) -- (0,6) rest distinct.
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,6) possibly same

  LFs = for lab in REPS#76 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#77 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 possibly same.  ACTUALLY: they are identical!!  How did that get through?  Totally different polytopes, exact same L, F...

  LFs = for lab in REPS#78 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 possibly same.

  LFs = for lab in REPS#79 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 possibly same.

  LFs = for lab in REPS#80 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#81 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 possibly same

  LFs = for lab in REPS#82 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- (0,1), 2
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 possibly same

  LFs = for lab in REPS#83 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --(0,2), 1
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,2 possibly same

  LFs = for lab in REPS#84 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 possibly same

  LFs = for lab in REPS#85 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- 0,1 possibly same

  LFs = for lab in REPS#86 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- (0,3), (1,4), 2.
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,3) and (1,4) possibly same

  LFs = for lab in REPS#87 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- (0,2), 1
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,2) possibly same

  LFs = for lab in REPS#88 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#89 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- all 3 different
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#90 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 0, (1,2)
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (1,2) could be same.

  LFs = for lab in REPS#91 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#92 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#93 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#94 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1,2) could be same.

  LFs = for lab in REPS#95 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) --
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#96 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- different
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#97 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- different
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#98 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- different
  netList for a in LFs list (factor det hessian toQQ a#1) -- 

  LFs = for lab in REPS#99 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- (0,2), 1
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,2) could be same.

  LFs = for lab in REPS#100 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#101 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1,2) could be same.

  LFs = for lab in REPS#102 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#103 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#104 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#105 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#106 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  LFs = for lab in REPS#107 list (X := Xs#lab; {c2Form X, cubicForm X});
  netList for a in LFs list (decompose sing toQQ a#1) -- 
  netList for a in LFs list (betti res sing toQQ a#1) -- 
  netList for a in LFs list (saturate sing2 ideal a) -- 
  netList for a in LFs list (factor det hessian toQQ a#1) -- (0,1) could be same.

  
  -- REP#8
    LFs = for lab in REPS#8 list (X := Xs#lab; {c2Form X, cubicForm X})
    F0 = LFs#0#1
    F1 = LFs#1#1
    F2 = LFs#2#1
    F3 = LFs#3#1
    F4 = LFs#4#1
    L0 = LFs#0#0
    L1 = LFs#1#0
    L2 = LFs#2#0
    L3 = LFs#3#0
    L4 = LFs#4#0

    -- STEP 1 singular locus over QQ.  All have 1 singular point over QQ.
      decompose sing toQQ F0
      decompose sing toQQ F2
      decompose sing toQQ F4
      
      decompose sing toQQ F1 -- different

      decompose sing toQQ F3 -- different

      betti res inverseSystem toQQ F0
      betti res inverseSystem toQQ F2
      betti res inverseSystem toQQ F4

      sing2 = I -> trim(I + minors(2, jacobian I))
      s0 = saturate sing2 ideal LFs#0
      s1 = saturate sing2 ideal LFs#2
      s2 = saturate sing2 ideal LFs#4
      
      factor det hessian toQQ F0
      factor det hessian toQQ F2
      factor det hessian toQQ F4

      -- promising: each of these have 5 ppints.
      decompose minors(3, hessian toQQ F0)
      decompose minors(3, hessian toQQ F2)
      decompose minors(3, hessian toQQ F4)

  ----------------------------------------------------
  -- NOT DONE YET
    positions(REPS, x -> #x == 6) == {71}
    LFs = for lab in REPS#71 list (X := Xs#lab; {c2Form X, cubicForm X})
    F0 = LFs#0#1
    F1 = LFs#1#1
    F2 = LFs#2#1
    F3 = LFs#3#1
    F4 = LFs#4#1
    F5 = LFs#5#1
    L0 = LFs#0#0
    L1 = LFs#1#0
    L2 = LFs#2#0
    L3 = LFs#3#0
    L4 = LFs#4#0
    L5 = LFs#5#0

    -- STEP 1 singular locus over QQ.  All have 1 singular point over QQ.
      decompose sing toQQ F0
      decompose sing toQQ F1
      decompose sing toQQ F2
      decompose sing toQQ F3
      decompose sing toQQ F4
      decompose sing toQQ F5 

      -- (F0,F1,F2,F3,F4,F5) doesn't separate them at all.
   -- STEP 2: inverse systems
      betti res inverseSystem toQQ F0
      betti res inverseSystem toQQ F1
      betti res inverseSystem toQQ F2
      betti res inverseSystem toQQ F3
      betti res inverseSystem toQQ F4
      betti res inverseSystem toQQ F5

      -- all the same still!

   -- STEP 3. jacobian over ZZ
      sing2 = I -> trim(I + minors(2, jacobian I))
      s0 = saturate sing2 ideal LFs#0
      s1 = saturate sing2 ideal LFs#1
      s2 = saturate sing2 ideal LFs#2
      s3 = saturate sing2 ideal LFs#3
      s4 = saturate sing2 ideal LFs#4
      s5 = saturate sing2 ideal LFs#5

      s0_0, s1_0, s3_0
      s2_0, s5_0
      s4_0
      -- (F0,F1,F3), (F2,F5), F4.  
      -- NOT DONE: the stuff below this doesn't yet separate these 3.

      Rp = ZZ/467[a,b,c,d]
      saturate sub(sing2 ideal LFs#2, Rp)
      saturate sub(sing2 ideal LFs#5, Rp)

      factor det hessian toQQ F0
      factor det hessian toQQ F1
      factor det hessian toQQ F3

      factor det hessian toQQ F2
      factor det hessian toQQ F5

      factor det hessian toQQ F4
      
      -- TODO: still need to separate both of these groups.
      decompose minors(2, hessian toQQ F0)
      decompose minors(2, hessian toQQ F1)
      decompose minors(2, hessian toQQ F3)
      
      decompose minors(2, hessian toQQ F2)
      decompose minors(2, hessian toQQ F5)
      s0 = gens gb saturate sing F0
      s1 = gens gb saturate sing F1
      s1 = saturate sing2 ideal LFs#1
      s2 = saturate sing2 ideal LFs#2
      s3 = saturate sing2 ideal LFs#3
      s4 = saturate sing2 ideal LFs#4
      s5 = saturate sing2 ideal LFs#5

  ---------------------------------------------  
  -- REPS#55
    positions(REPS, x -> #x == 7) == {55, 75}
    LFs = for lab in REPS#55 list (X := Xs#lab; {c2Form X, cubicForm X})
    F0 = LFs#0#1
    F1 = LFs#1#1
    F2 = LFs#2#1
    F3 = LFs#3#1
    F4 = LFs#4#1
    F5 = LFs#5#1
    F6 = LFs#6#1
    L0 = LFs#0#1
    L1 = LFs#1#0
    L2 = LFs#2#0
    L3 = LFs#3#0
    L4 = LFs#4#0
    L5 = LFs#5#0
    L6 = LFs#6#0

    -- STEP 1 singular locus over QQ
      decompose sing toQQ F0
      decompose sing toQQ F1
      decompose sing toQQ F3
      decompose sing toQQ F4
      
      decompose sing toQQ F2

      decompose sing toQQ F5 -- 2 pts, irred over QQ
      decompose sing toQQ F6
      -- (F0,F1,F3,F4), F2, (F5,F6)
      
   -- STEP 2 inverse system
      betti res inverseSystem toQQ F0
      betti res inverseSystem toQQ F1 -- distinct from F0,F3,F4.
      betti res inverseSystem toQQ F3
      betti res inverseSystem toQQ F4
      
      betti res inverseSystem toQQ F2
      
      betti res inverseSystem toQQ F5
      betti res inverseSystem toQQ F6 -- F6 distinct from F5.
      -- (F0,F3,F4), F1, F2, F5, F6.

   -- STEP 3. jacobian over ZZ
      sing2 = I -> trim(I + minors(2, jacobian I))
      s0 = saturate sing2 ideal LFs#0
      s3 = saturate sing2 ideal LFs#3
      s4 = saturate sing2 ideal LFs#4

      s5 = saturate sing2 ideal LFs#5

      s6 = saturate sing2 ideal LFs#6

      s1 = saturate sing2 ideal LFs#1 -- different from F0

      s2 = saturate sing2 ideal LFs#2

      s0_0, s3_0, s4_0, s1_0, s2_0, s5_0, s6_0
      -- (F0,F3,F4), F1, F2, F5, F6.  Same as step 1+2.

   -- STEP 4. Separate F0, F3, F4
      factor det hessian toQQ F0
      factor det hessian toQQ F3
      factor det hessian toQQ F4 -- all similar structure...

      -- first compare F0, F3: F0, F3 are equiv over QQ, not ZZ.
      (A, phi) = genericLinearMap RQ
      use target phi
      start1 = {{b+3*c, a+d},
          {3*c+2*d, a},
          {a+2*c, a+b-c+d}}
      start2 = {{b+3*c, a+d},
          {3*c+2*d, a+b-c+d},
          {a+2*c, a}}
      signs = (toList((set{-1,1}) ** set{-1,1} ** set{-1,1}))/splice/toList -- note: not all these signs are needed...
      chsigns = (L, sgn) -> for i from 0 to #L-1 list {L#i#0, sgn#i * L#i#1}
      set1 = for sgn in signs list chsigns(start1, sgn)
      set2 = for sgn in signs list chsigns(start2, sgn)
      allsets = join(set1, set2)
      netList for L in allsets list (
          (A0, phi0) = linearEquationConstraints(A, phi, append(L, {sub(F0, target phi), sub(F3, target phi)}), {});
          if A0 == 0 then continue;
          phi1 = map(RQ, RQ, transpose sub(A0, QQ));
          {A0, det A0, phi1 toQQ(LFs#0#0) - toQQ LFs#3#0, phi1 toQQ LFs#0#1 - toQQ LFs#3#1}
          )

      -- second compare F0, F4: -- to see yet: F0, F4 are equiv over ZZ, but X0, X4 are equiv over QQ.
      (A, phi) = genericLinearMap RQ
      use target phi
      factor det hessian toQQ F0
      factor det hessian toQQ F4
      start1 = {{b+3*c, c+d},
          {3*c+2*d, d},
          {a+2*c, a-b-c}}
      start2 = {{b+3*c, c+d},
          {3*c+2*d, a-b-c},
          {a+2*c, d}}
      signs = (toList((set{-1,1}) ** set{-1,1} ** set{-1,1}))/splice/toList -- note: not all these signs are needed...
      chsigns = (L, sgn) -> for i from 0 to #L-1 list {L#i#0, sgn#i * L#i#1}
      set1 = for sgn in signs list chsigns(start1, sgn)
      set2 = for sgn in signs list chsigns(start2, sgn)
      allsets = join(set1, set2)
      netList for L in allsets list (
          (A0, phi0) = linearEquationConstraints(A, phi, append(L, {sub(F0, target phi), sub(F4, target phi)}), {});
          if A0 == 0 then continue;
          phi1 = map(RQ, RQ, transpose sub(A0, QQ));
          {A0, det A0, phi1 toQQ(LFs#0#0) - toQQ LFs#4#0, phi1 toQQ LFs#0#1 - toQQ LFs#4#1}
          )

      -- second compare F0 to F3 -- only equiv over QQ.
      (A, phi) = genericLinearMap RQ
      use target phi
      factor det hessian toQQ F3
      factor det hessian toQQ F4
      start1 = {{a+d, c+d},
          {a, d},
          {a+b-c+d, a-b-c}}
      start2 = {{a+d, c+d},
          {a, a-b-c},
          {a+b-c+d, d}}
      signs = (toList((set{-1,1}) ** set{-1,1} ** set{-1,1}))/splice/toList -- note: not all these signs are needed...
      chsigns = (L, sgn) -> for i from 0 to #L-1 list {L#i#0, sgn#i * L#i#1}
      set1 = for sgn in signs list chsigns(start1, sgn)
      set2 = for sgn in signs list chsigns(start2, sgn)
      allsets = join(set1, set2)
      netList for L in allsets list (
          (A0, phi0) = linearEquationConstraints(A, phi, append(L, {sub(F3, target phi), sub(F4, target phi)}), {});
          if A0 == 0 then continue;
          phi1 = map(RQ, RQ, transpose sub(A0, QQ));
          {A0, det A0, phi1 toQQ(LFs#3#0) - toQQ LFs#4#0, phi1 toQQ LFs#3#1 - toQQ LFs#4#1}
          )
  
    -- Upshot: these 7 are all distinct.

  -- REPS#75
    positions(REPS, x -> #x == 7) == {55, 75}
    LFs = for lab in REPS#75 list (X := Xs#lab; {toQQ c2Form X, toQQ cubicForm X})
    F1 = LFs#0#1
    F2 = LFs#1#1
    F3 = LFs#2#1
    F4 = LFs#3#1
    F5 = LFs#4#1
    F6 = LFs#5#1
    F7 = LFs#6#1

    LFZs = for lab in REPS#75 list (X := Xs#lab; {c2Form X, cubicForm X})
    FZ1 = LFZs#0#1
    FZ2 = LFZs#1#1
    FZ3 = LFZs#2#1
    FZ4 = LFZs#3#1
    FZ5 = LFZs#4#1
    FZ6 = LFZs#5#1
    FZ7 = LFZs#6#1

    -- STEP 1: consider singular loci over QQ.
      -- one singular point
      decompose sing F1
      decompose sing F5
      decompose sing F7
      -- two singular points
      decompose sing F4
      decompose sing F2
      decompose sing F3
      -- three singular points
      decompose sing F6
      -- at this point: (F1,F5,F7), (F2,F3,F4), F6.  ones in parens might be equivalent.
      
    -- STEP 2: inverse systems (over QQ)
      betti res inverseSystem F1            
      betti res inverseSystem F5
      betti res inverseSystem F7

      betti res inverseSystem F4 -- F4 on its own

      betti res inverseSystem F2
      betti res inverseSystem F3 -- F2, F3 have different singular loci mod 19.

      betti res inverseSystem F6 -- F6
      -- at this point: (F1,F5,F7), F4, (F2,F3), F6.  ones in parens might be equivalent.

    -- STEP 3: singular locus of (L,F) over ZZ
      sing2 = I -> trim(I + minors(2, jacobian I))
      s1 = saturate sing2 ideal LFZs#0
      s5 = saturate sing2 ideal LFZs#4 -- this shows that F5 is different from (F1, F7).
      s7 = saturate sing2 ideal LFZs#6

      s4 = saturate sing2 ideal LFZs#3
      s2 = saturate sing2 ideal LFZs#1
      s3 = saturate sing2 ideal LFZs#2 -- this shows that F2 is different from F3.
        -- (it also show separately that F4 is distinct from F2 and F3).
      -- at this point: (F1,F7), F5, F4, F2, F3, F6.
    
    -- STEP 4: separate F1, F7.  This one is a bit tricky, as in fact F1 and F7 are equivalent over ZZ!
    -- But they are not equivalent when the c2 form is taken into account.
    -- For this one, we compute the Hessian of F1, F7.
      factor det hessian F1 -- (b+2*d)*(b+3*d)^2*(a+b+2*d)*(-20736)
      factor det hessian F7 -- (a-d)*(a+c)^2*(a-b)*(-20736)
      -- this shows that we seek a matrix over ZZ (or QQ, if we are interested)
      -- with b+3d --> \pm (a+c)
      -- with b+2d --> \pm (a-d)  OR \pm (a-b)
      -- with a+b+2d --> \pm (a-d)  OR \pm (a-b) [but for the other linear form].
      -- note: this fixes phi(a), phi(b), phi(d), leaving phi(c) not known.
      -- and F1 is linear in c.
      -- all the choices:
      (A, phi) = genericLinearMap RQ
      use target phi
      start1 = {{b+3*d, a+c}, {b+2*d, a-d}, {a+b+2*d, a-b}}
      start2 = {{b+3*d, a+c}, {b+2*d, a-b}, {a+b+2*d, a-d}}
      signs = (toList((set{-1,1}) ** set{-1,1} ** set{-1,1}))/splice/toList -- note: not all these signs are needed...
      chsigns = (L, sgn) -> for i from 0 to #L-1 list {L#i#0, sgn#i * L#i#1}
      set1 = for sgn in signs list chsigns(start1, sgn)
      set2 = for sgn in signs list chsigns(start2, sgn)
      allsets = join(set1, set2)
      netList for L in allsets list (
          (A0, phi0) = linearEquationConstraints(A, phi, append(L, {sub(F1, target phi), sub(F7, target phi)}), {});
          if A0 == 0 then continue;
          phi1 = map(RQ, RQ, transpose sub(A0, QQ));
          {A0, det A0, phi1 (LFs#0#0) - LFs#6#0, phi1 (LFs#0#1) - LFs#6#1}
          )
      -- this shows that F1, F7 are in fact equivalent over ZZ, but X1, X7 are distinct topologically.    
      -- it also shows there is a matrix over QQ, with det 1, which does map L1 to L7, F1 to F7.
///

///
  -- Example use of a constructed data base.
  restart
  debug needsPackage "StringTorics"
  RZ = ZZ[a,b,c]
  F = openDatabase "../m2-examples/foo-h11-3.dbm" -- note: 
    Xlabels = sort select(keys F, k -> (a := value k; instance(a, Sequence)))
    Qlabels = sort select(keys F, k -> (a := value k; instance(a, ZZ)))
    assert(#Xlabels == 526)
    assert(#Qlabels == 244)
    elapsedTime Qs = for k in Qlabels list cyPolytope(F#k, ID => value k);
    elapsedTime Xs = for k in Xlabels list cyData(F#k, i -> Qs#i, Ring => RZ);
    assert(Xs/label === Xlabels/value)
  close F

  Xs = hashTable for x in Xs list (label x) => x;  

  -- these are the ones we will consider.
  torsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V === ZZ^3 then lab else continue
      )

  -- 6 polytopes are not torsion free (i.e. classGroup is not torsion free)
  nontorsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V != ZZ^3 then lab else continue
      )
  nontorsionfrees == {(0, 0), (9, 0), (9, 1), (10, 0), (10, 1), (10, 2), (55, 0), (62, 0), (232, 0)}

  nonFavorables = select(sort keys Xs, lab -> not isFavorable cyPolytope Xs#lab)
  favorables = select(sort keys Xs, lab -> isFavorable cyPolytope Xs#lab)

  elapsedTime H = partition(lab -> (
          X := Xs#lab;
          elapsedTime netList join({"h11" => hh^(1,1) X, "h12" => hh^(1,2) X}, invariantsAll Xs#lab)
          ), sort torsionfrees); --  18 sec
  #(keys H) == 195
  H


  elapsedTime INV = for k in sort keys H list k => partitionByTopology(H#k, Xs, 15); --34 sec

  REPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then continue
      else k#0 => K)
  #REPS == 9

  positions(INV, k -> #(keys (k#1)) > 1)
  netList INV_{26, 56, 63, 64, 76, 120, 121, 156, 182}
  REPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then continue
      else K)

  X1 = Xs#(34,0)
  X2 = Xs#(37,0)
  L1 = c2Form X1
  L2 = c2Form X2
  F1 = cubicForm X1
  F2 = cubicForm X2
  
  hashTable invariantsAll X1, hashTable invariantsAll X2



  -- XXXX
  elapsedTime H1 = partition(lab -> topologicalData Xs#lab, torsionfrees); -- take one from each.
  netList (values H1)
  SAME = hashTable for k in keys H1 list (first H1#k) => drop(H1#k, 1)

  labelYs = sort keys SAME -- these have non-equal topological data 
  #labelYs == 286 -- if torsionfrees is used
  #labelYs == 291 -- all are equivalent to one of these, so there are 291 possible different topologies
    -- (but probably less than this).
    
  --elapsedTime INV = partition(lab -> invariants Xs#lab, labelYs); -- 98 sec (FIXME: about 200 seconds now, with slow pointCount code)
  elapsedTime INV = partition(lab -> invariants2 Xs#lab, labelYs); -- 6 seconds

  #INV == 165 -- torsionfrees, invariants2.
  -- #INV == 168 -- for both invariants, invariants2
  
  -- This selects the different topologies, except it can't tell about (52,0), (53,0).
  (keys INV)/(x -> #INV#x)//tally

  RESULT1 = hashTable for k in sort keys INV list (
      labs := INV#k;
      k => partitionByTopology(labs, Xs, 15)
      )

  labs = INV#{3, 99, 0, 2, 1, 1, 1, 1}
  gvInvariants(Xs#(labs#0), DegreeLimit => 15)
  gvCone(Xs#(labs#0), DegreeLimit => 15)
  
  TODO = select(keys RESULT1, k -> #(keys RESULT1#k) > 1)
  DONE = select(keys RESULT1, k -> #(keys RESULT1#k) == 1)
  hashTable for k in DONE list k => RESULT1#k
  TODO = hashTable for k in TODO list k => RESULT1#k

  -- try comparing (26,0), (37,0)
  partitionByTopology({(26,0), (37,0)}, Xs, 25)
  partitionGVConeByGV(Xs#(26,0), DegreeLimit => 20)
  partitionGVConeByGV(Xs#(37,0), DegreeLimit => 20) -- this makes them appear distinct.
  F26 = cubicForm Xs#(26,0)
  F37 = cubicForm Xs#(37,0)
  L26 = c2Form Xs#(26,0)
  L37 = c2Form Xs#(37,0)
  
  RQ = QQ[a,b,c]
  F1 = sub(F26, RQ)
  L1 = sub(L26, RQ)
  F2 = sub(F37, RQ)
  L2 = sub(L37, RQ)
  decompose ideal(F1, L1)
  decompose ideal(F2, L2)

  mons3 = basis(3, RQ)
  coeffsF1 = flatten entries last coefficients(F1, Monomials => mons3)
  coeffsF2 = flatten entries last coefficients(F2, Monomials => mons3)
  needsPackage "EllipticCurves"
  toWeierstrass(RingElement, List) := (F, pt) -> (
      R := ring F;
      kk := coefficientRing R;
      if numgens R =!= 3 or #pt != 3 then error "expected ring in 3 variables and point in 3-space";
      mons3 := basis(3, R);
      coeffsF := flatten entries last coefficients(F, Monomials => mons3);
      toWeierstrass(coeffsF, pt, kk)
      )
  jInvariant toWeierstrass(F1, {2, 3, -4})
  jInvariant toWeierstrass(F2, {1, 1, -1})

  TOANALYZE = for k in sort keys RESULT1 list (
      if #RESULT1#k == 1 then continue;
      k => (FINISH ME !!!!!!!!!!!!!!!!!!!!!!!!
      
  select(sort keys RESULT1, k -> #(keys RESULT1#k) > 1)

  tally for k in keys RESULT1 list #RESULT1#k -- all 1's, meaning that in each case, all elements of 
  -- INV#k are equivalent.
  
  ALLTOPS = (keys RESULT1)/(k -> (
          for x in keys RESULT1#k list x => join(RESULT1#k#x, SAME#x)
          ))//flatten//hashTable

///

----------------------------------------------------------
-- Newer version of last example -------------------------
----------------------------------------------------------


///
  -- Example use of a constructed data base.
  restart
  debug needsPackage "StringTorics"
  RZ = ZZ[a,b,c]
  F = openDatabase "../m2-examples/foo-ntfe-h11-3.dbm" -- note: 
    Xlabels = sort select(keys F, k -> (a := value k; instance(a, Sequence)))
    Qlabels = sort select(keys F, k -> (a := value k; instance(a, ZZ)))
    assert(#Xlabels == 306)
    assert(#Qlabels == 244)
    elapsedTime Qs = for k in Qlabels list cyPolytope(F#k, ID => value k);
    elapsedTime Xs = for k in Xlabels list cyData(F#k, i -> Qs#i, Ring => RZ);
    assert(Xs/label === Xlabels/value)
  close F

  Xs = hashTable for x in Xs list (label x) => x;  

  -- these are the ones we will consider.
  torsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V === ZZ^3 then lab else continue
      )

  -- 6 polytopes are not torsion free (i.e. classGroup is not torsion free)
  nontorsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V != ZZ^3 then lab else continue
      )
  nontorsionfrees == {(0, 0), (9, 0), (10, 0), (55, 0), (62, 0), (232, 0)}

  nonFavorables = select(sort keys Xs, lab -> not isFavorable cyPolytope Xs#lab)
  favorables = select(sort keys Xs, lab -> isFavorable cyPolytope Xs#lab)

  -- XXXX
  elapsedTime H1 = partition(lab -> topologicalData Xs#lab, torsionfrees); -- take one from each.
  netList (values H1)
  SAME = hashTable for k in keys H1 list (first H1#k) => drop(H1#k, 1)

  labelYs = sort keys SAME -- these have non-equal topological data 
  #labelYs == 286 -- if torsionfrees is used
  #labelYs == 291 -- all are equivalent to one of these, so there are 291 possible different topologies
    -- (but probably less than this).
    
  --elapsedTime INV = partition(lab -> invariants Xs#lab, labelYs); -- 98 sec (FIXME: about 200 seconds now, with slow pointCount code)
  elapsedTime INV = partition(lab -> invariants2 Xs#lab, labelYs); -- 6 seconds

  #INV == 165 -- torsionfrees, invariants2.
  -- #INV == 168 -- for both invariants, invariants2
  
  -- This selects the different topologies, except it can't tell about (52,0), (53,0).
  (keys INV)/(x -> #INV#x)//tally

  RESULT1 = hashTable for k in sort keys INV list (
      labs := INV#k;
      k => partitionByTopology(labs, Xs, 15)
      )

  labs = INV#{3, 99, 0, 2, 1, 1, 1, 1}
  gvInvariants(Xs#(labs#0), DegreeLimit => 15)
  gvCone(Xs#(labs#0), DegreeLimit => 15)
  
  TODO = select(keys RESULT1, k -> #(keys RESULT1#k) > 1)
  DONE = select(keys RESULT1, k -> #(keys RESULT1#k) == 1)
  hashTable for k in DONE list k => RESULT1#k
  TODO = hashTable for k in TODO list k => RESULT1#k

  -- try comparing (26,0), (37,0)
  partitionByTopology({(26,0), (37,0)}, Xs, 25)
  partitionGVConeByGV(Xs#(26,0), DegreeLimit => 20)
  partitionGVConeByGV(Xs#(37,0), DegreeLimit => 20) -- this makes them appear distinct.
  F26 = cubicForm Xs#(26,0)
  F37 = cubicForm Xs#(37,0)
  L26 = c2Form Xs#(26,0)
  L37 = c2Form Xs#(37,0)
  
  RQ = QQ[a,b,c]
  F1 = sub(F26, RQ)
  L1 = sub(L26, RQ)
  F2 = sub(F37, RQ)
  L2 = sub(L37, RQ)
  decompose ideal(F1, L1)
  decompose ideal(F2, L2)

  mons3 = basis(3, RQ)
  coeffsF1 = flatten entries last coefficients(F1, Monomials => mons3)
  coeffsF2 = flatten entries last coefficients(F2, Monomials => mons3)
  needsPackage "EllipticCurves"
  toWeierstrass(RingElement, List) := (F, pt) -> (
      R := ring F;
      kk := coefficientRing R;
      if numgens R =!= 3 or #pt != 3 then error "expected ring in 3 variables and point in 3-space";
      mons3 := basis(3, R);
      coeffsF := flatten entries last coefficients(F, Monomials => mons3);
      toWeierstrass(coeffsF, pt, kk)
      )
  jInvariant toWeierstrass(F1, {2, 3, -4})
  jInvariant toWeierstrass(F2, {1, 1, -1})

  TOANALYZE = for k in sort keys RESULT1 list (
      if #RESULT1#k == 1 then continue;
      k => (FINISH ME !!!!!!!!!!!!!!!!!!!!!!!!
      
  select(sort keys RESULT1, k -> #(keys RESULT1#k) > 1)

  tally for k in keys RESULT1 list #RESULT1#k -- all 1's, meaning that in each case, all elements of 
  -- INV#k are equivalent.
  
  ALLTOPS = (keys RESULT1)/(k -> (
          for x in keys RESULT1#k list x => join(RESULT1#k#x, SAME#x)
          ))//flatten//hashTable

///


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

///

  restart
  debug needsPackage "StringTorics"
  DIRNAME = "~/Dropbox/Collaboration/Physics-Liam/Inequivalent CYs/cy_classes/"

  Ps = readPolytopes(DIRNAME|"polytopes_h11=2.dat")
  Ss = readSimplices(DIRNAME|"simplices_h11=2.dat")
  Es = readEquivalences(DIRNAME|"cy_classes_h11=2.dat")
  cyPolytope(Ps, 3)

  readPolytopes(DIRNAME|"polytopes_h11=3.dat")
  readSimplices(DIRNAME|"simplices_h11=3.dat")
  readEquivalences(DIRNAME|"cy_classes_h11=3.dat")

  Ps = readPolytopes(DIRNAME|"polytopes_h11=4.dat");
  Ss = readSimplices(DIRNAME|"simplices_h11=4.dat");
  Es = readEquivalences(DIRNAME|"cy_classes_h11=4.dat");

  Ps = readPolytopes(DIRNAME|"polytopes_h11=5.dat");
  Ss = readSimplices(DIRNAME|"simplices_h11=5.dat");
  Es = readEquivalences(DIRNAME|"cy_classes_h11=5.dat");

  Qs = new MutableHashTable
  Xs = new MutableHashTable
  RZ = ZZ[a,b,c,d]
  elapsedTime X = makeCY(13, (Ps, Ss), Qs, Xs, Ring => RZ) -- this sets Qs#899, Xs#13.
  elapsedTime X = makeCY(12, (Ps, Ss), Qs, Xs, Ring => RZ) -- this sets Qs#899, Xs#13.
  peek Qs#899  .cache
  
  elapsedTime Qs = for lab in sort keys Ps list cyPolytope(Ps#lab, ID => lab);

  Ps#3,  rays Q -- these are in different orders!
  Q = cyPolytope(Ps#3, ID => 3)
  latticePointList polytope(Q, "N")
  rays Q
  calabiYau(13, Ps, Ss) -- 
  Ss#13
  Q = cyPolytope(Ps#(Ss#13#0), ID => Ss#13#0)
  label Q === 899
  rays Q
  Q.cache#"face dimensions"
  Ps#899  
  netList annotatedFaces Q
  get (DIRNAME|"simplices_h11=2.dat")
  get (DIRNAME|"cy_classes_h11=4.dat")
  Q = cyPolytope oo
  hh^(1,2) Q
  isReflexive polytope Q
  vertices polytope Q
  netList annotatedFaces Q
///

makeCY(ZZ, Sequence, MutableHashTable, MutableHashTable) := opts -> (labX, PSs, Qs, Xs) -> (
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
    
