---------------------------------------------------
-- Creating databases for ReflexivePolytope type --
---------------------------------------------------
doc ///
  Key
    "Creating a CYDatabase file for h11=2"
  Headline
    How to create a CYDatabase of reflexive polytopes and corresponding Calabi-Yau 3-folds
  Description
    Text
      A CYDatabase file is a @TO Database@, which contains pre-computed information about
      4-dimensional reflexive polytopes, and their triangulations.  The fact that it contains
      precomputed features allows us to spend much less time scanning over many polytopes.

      Here, we describe how to construct such a file.  We will construct the file for all examples
      with $h^{1,1} = 2$.  There are not very many of these, but the same method works for
      larger cases.  We limit the number below, in order to make the
      example run faster.  You would generally not place a limit.
    Example
      topes = kreuzerSkarke(2, Limit => 10);
      assert(#topes == 10) -- without the limit, it would be 36.
      elapsedTime addToCYDatabase("can-delete-me-ntfe-h11-2.dbm", topes)
    Text
      Let's test that this was created correctly.  We see that in particular the
      annotated faces (@TO (annotatedFaces, ReflexivePolytope)@) has been computed (this is
      one of the things that seems to take the longest.  That, and the list of triangulations.
      For higher $h^{1,1}(X)$, we must arrange to not compute these, as there are too many triangulations.
      This has not been done yet.

      The ring of a @TO CalabiYauInToric@ is the intersection ring: a ring over $\ZZ$ in
      $h^{1,1}(X)$ variables.  This is the same for all Calabi-Yaus with the same $h^{1,1}(X)$, so we
      create one such ring, and use it for all constructed Calabi-Yau's.
    Example
      RZ = ZZ[a,b]
      (Qs, Xs) = readCYDatabase("can-delete-me-ntfe-h11-2.dbm", Ring => RZ);
      assert(keys Qs === toList(0..9))
      sort keys Xs
      #oo == 10
    Text
      The line above confirms that there are 10 constructed Calabi-Yau 3-fold hypersurfaces of
      $h^{1,1} = 2$.  Note that each polytope only gives one triangulation.  As  $h^{1,1}$ increases,
      the number of triangulations increases, and at some point, becomes astronomical.

      For now, we take one, and look at some of its invariants.
    Example
      X = Xs#(0,0)
      label X
      hh^(1,1) X
      hh^(1,2) X
      cubicForm X
      c2Form X
    Text
      In general you would not run the next line, but since this is a test, we will do so here.
    Example
      removeFile "can-delete-me-ntfe-h11-2.dbm"
    Text
      Now let's deal with a larger case, when we want to have multiple processes working on creating
      databases, and then we will merge them together.
    Pre
      debug needsPackage "StringTorics"
      DBNAME = "can-delete-me-cys-ntfe-h11-5.dbm"
      topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
      assert(#topes == 4990)
      "topes-h11-5.txt" << toExternalString topes << endl << close;
      topes2 = value get "topes-h11-5.txt";
      assert(topes === topes2)
      createM2Lines("cys-h11-5", "topes-h11-5.txt", 4990, 22)
      elapsedTime addToCYDatabase(DBNAME, topes) -- on Apple M4 Max, Jan 2025: 3.14 hours to create.
    Text
      Here we try to run lots of M2's to do this.
    Pre
      -- TODO: this is WRONG!!
      needsPackage "StringTorics"
      R = ZZ[a..f]
      filename = elapsedTime processCYPolytopes("cys-h11-6", "topes-h11-6", 2, 3000)
      filename = elapsedTime addToCYDatabase("cys-h11-6", "topes-h11-6", 2, 3000)
      filename = elapsedTime addToCYDatabase("cys-h11-6", "topes-h11-6", 2999, 3000)
      (Qs, Xs) = readCYDatabase(filename, Ring => ZZ[a_0..a_5]);

       M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5' -e 'hi=8' -e 'processCYPolytopes("cys-h11-6", "topes-h11-6", (lo,hi))' -e 'exit 0'
       M2 --silent --stop -e 'needsPackage "StringTorics"' -e 'lohi = (5,8)' -e 'processCYPolytopes("cys-h11-6", "topes-h11-6", lohi)' -e 'exit 0'

      filename = elapsedTime processCYPolytopes("cys-h11-6", "topes-h11-6", (5,8))
      (Qs, Xs) = readCYDatabase(filename, Ring => ZZ[a_0..a_5]);
      
      topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
      assert(#topes == 4990)
      elapsedTime("topes-h11-5" << toExternalString topes << close)
      elapsedTime get "topes-h11-5";
      elapsedTime value oo;
      oo === topes

      topes = kreuzerSkarke(6, Limit => 20000); -- 17101 of these
      assert(#topes == 17101)
      elapsedTime("topes-h11-6" << toExternalString topes << close);
      topes2 = elapsedTime value get "topes-h11-6";


      
      
      DBNAME = "dbm-h11-6-range-5-20"
      elapsedTime addToCYDatabase(DBNAME, topes_{5..20})
      M2 --stop -e 'needsPackage "StringTorics" -e 'lo=5' -e 'hi=10' -e 'addToCYDatabase(DBNAME|"-range-"|lo|"-"|hi, topes

      cyDatabase(DBNAMEPREFIX, topesFile, lo, hi);
  SeeAlso
    addToCYDatabase
    readCYDatabase
///

///
   -- working on whether our creation of Q, X from the data base is
   -- doing more work than should be.  i.e. we should just grab data
   -- from the database, not doing any further computation (the N
   -- polytope in particular should not be created, I would hope...
   restart
   needsPackage "StringTorics"
   F = openDatabase "can-delete-me-ntfe-h11-2.dbm"
   F#"0"
   Q = reflexivePolytope F#"0"
   peek Q
   peek Q.cache
   read

   Q = cyPolytope(F, 3)
   peek Q.cache
   X = calabiYau(F, (3,0))
 ///

/// -- working on this one 12 Mar 2025.
  Key
    readCYDatabase
    (readCYDatabase, String)
    [readCYDatabase, Ring]
  Headline
    read in all ReflexivePolytope's and CalabiYauInToric's
  Usage
    (Qs, Xs) = readCYDatabase dbname
    (Qs, Xs) = readCYDatabase(dbname, Ring => ZZ[t_1..t_n])
  Inputs
    dbname:String
      a file name, containing a CY Database
  Outputs
    :Sequence
      of two item: the first is a hash table of all ReflexivePolytope's.
      The second is a hash table of all CalabiYauInToric's
  Description
    Text
    Example
      (Qs, Xs) = readCYDatabase("cy3-h11-3.dbm", Ring => (R = ZZ[a,b,c]));
      #(keys Qs) == 244
      #(keys Xs) == 275
    Example
      readCYPolytopes "cy3-h11-3.dbm";
      readCYs("cy3-h11-3.dbm", Qs, Ring => R);
  SeeAlso
    
///



///
restart
needsPackage "StringTorics"
-- code to create a data base for a specific h11.
-- 1. first creation.
-- 2. adding new data when we make it later.

    ------------
    -- Step 1 --
    -- create text files of the h11=d polytopes, 1 <= d <= 7 (only 5 so far).
    ------------
    topes = kreuzerSkarke(2, Limit => 20000);
    assert(#topes == 36)
    "topes-h11-2.txt" << toExternalString topes << endl << close;
    
    topes = kreuzerSkarke(3, Limit => 20000);
    assert(#topes == 244)
    "topes-h11-3.txt" << toExternalString topes << endl << close;

    topes = kreuzerSkarke(4, Limit => 20000);
    assert(#topes == 1197)
    "topes-h11-4.txt" << toExternalString topes << endl << close;

    topes = kreuzerSkarke(5, Limit => 20000);
    assert(#topes == 4990)
    "topes-h11-5.txt" << toExternalString topes << endl << close;

    topes = kreuzerSkarke(6, Limit => 200000);
    assert(#topes == 17101)
    "topes-h11-6.txt" << toExternalString topes << endl << close;

    topes = kreuzerSkarke(7, Limit => 200000);
    assert(#topes == 50376)
    "topes-h11-7.txt" << toExternalString topes << endl << close;
    
    elapsedTime value get "topes-h11-7.txt"; -- < .2 sec

   ------------
   -- Step 2 --
   ------------
   -- h11=2
   createCYDatabaseFiles("cy3-h11-2", "topes-h11-2.txt", 15)
   -- now wait for all 15 of these M2 processes to stop (they display when the are done).
   combineCYDatabaseFiles("cy3-h11-2", "topes-h11-2.txt", 15)

   -- h11=3
   createCYDatabaseFiles("cy3-h11-3", "topes-h11-3.txt", 15)
   -- now wait for all 15 of these M2 processes to stop (they discplay when the are done).
   combineCYDatabaseFiles("cy3-h11-3", "topes-h11-3.txt", 15)

   -- h11=4
   createCYDatabaseFiles("cy3-h11-4", "topes-h11-4.txt", 15)    
   combineCYDatabaseFiles("cy3-h11-4", "topes-h11-4.txt", 15)

   -- h11=5
   createCYDatabaseFiles("cy3-h11-5", "topes-h11-5.txt", 15)    
   combineCYDatabaseFiles("cy3-h11-5", "topes-h11-5.txt", 15)

   -- h11=6
   createCYDatabaseFiles("cy3-h11-6", "topes-h11-6.txt", 15)    
   combineCYDatabaseFiles("cy3-h11-6", "topes-h11-6.txt", 15)

   time (Qs, Xs) = readCYDatabase "cy3-h11-6.dbm"; -- this is
   Q = cyPolytope("cy3-h11-6.dbm", 11)
   keys Q.cache -- "N polytope" is not here!
   member("N polytope", keys Q.cache)
   X = calabiYau("cy3-h11-6.dbm", Q, (11,0))
   member("N polytope", keys Q.cache) -- now it is here!
///

doc ///
  Key
    addToCYDatabase
    (addToCYDatabase, String, List)
  Headline
    create or append to a database file and populate it with CYPolytope's and possibly CalabiYauInToric's
  Usage
    addToCYDatabase(filename, topes)
  Inputs
    filename:String
      the desired name of the data base file.  If the file doesn't exist it is created,
      otherwise the name should be the name of an existing data base file, and this file
      is modified
    topes:List
      of @ofClass KSEntry@'s, a list of Kreuzer-Skarke type entries for some polytopes
    "CYs" => Boolean
      if true, then also all Calabi Yau hypersurfaces are computed and added to the database.
    NTFE => Boolean
      if true, then triangulations which are identical on the set of 2-faces are considered the
      same, and only one is placed into the data base.
  Consequences
    Item
      For each polytope corresponding to an entry in the {\tt topes} list,
      a @ofClass CYPolytope@ is created, and various information about it is computed
      and then stored in the data base file for later use
  Description
    Text
      A CYDatabase file is a database file whose contents are precomputed data
      about some @TO CYPolytope@'s and @TO CalabiYauInToric@'s.  Since some information takes
      non-trivial time to construct, we precompute this data, and then we can later pull up
      this data via the functions @TO readCYDatabase@, @TO "readCYPolytopes"@, and @TO readCYs@.
    Text
      The CYPolytope corresponding to each item of the {\tt topes} list is constructed
      and some basic data is computed (e.g. information about the faces of the polytopes, whether the
      polytope is favorable, and degree information about it.  This data is then stored in the
      database for later retrieval.
    Example
      filename = "foo-remove-me.dbm"
      if fileExists filename then removeFile filename
      topes = kreuzerSkarke(2, Limit => 4)
      addToCYDatabase(filename, topes_{1,2,3})
    Example
      F = openDatabase filename
      F#"1"
      Q = reflexivePolytope F#"1"
      hh^(1,1) Q
      hh^(1,2) Q
      isFavorable Q
    Text
      As a data base file, all keys of {\tt F} are strings, and the values are strings too.
    Example
      sort keys F
    Text
      Close the database file when done with it.
    Example
      close F
    Text
      For this example, we also delete this database file.
    Example
      removeFile filename
  SeeAlso
    addToCYDatabase
    readCYDatabase
    readCYs
///
