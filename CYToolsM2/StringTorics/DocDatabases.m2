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
      larger cases.
    Example
      topes = kreuzerSkarke(2, Limit => 1000);
      assert(#topes == 36)
      elapsedTime addToCYDatabase("can-delete-me-ntfe-h11-2.dbm", topes)
    Text
      Let's test that this was created correctly.  We see that in particular the
      annotated faces (@TO (annotatedFaces, CYPolytope)@) has been computed (this is
      one of the things that seems to take the longest.  That, and the list of triangulations.
      For higher $h^{1,1}(X)$, we must arrange to not compute these, as there are too many triangulations.
      This has not been done yet.

      The ring of a @TO CalabiYauInToric@ is the intersection ring: a ring over $\ZZ$ in
      $h^{1,1}(X)$ variables.  This is the same for all Calabi-Yaus with the same $h^{1,1}(X)$, so we
      create one such ring, and use it for all constructed Calabi-Yau's.
    Example
      RZ = ZZ[a,b]
      (Qs, Xs) = readCYDatabase("can-delete-me-ntfe-h11-2.dbm", Ring => RZ);
      assert(keys Qs === toList(0..35))
      sort keys Xs
      #oo == 36
    Text
      The line above confirms that there are 36 constructed Calabi-Yau 3-fold hypersurfaces of
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
   
///
