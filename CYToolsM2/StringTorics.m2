newPackage(
        "StringTorics",
        Version => "0.7", -- bumped on 12 April.
        Date => "12 April 2024",
        Authors => {
            {Name => "Mike Stillman", 
            Email => "mike@math.cornell.edu", 
            HomePage => "http://pi.math.cornell.edu/~mike"}
        },
        Headline => "toric variety functions for string theory",
        DebuggingMode => true,
        AuxiliaryFiles => true,
        PackageExports => {
            "FourTiTwo", -- where is this used?
            "SimplicialComplexes",
            "NormalToricVarieties",
            "Schubert2",
            "Polyhedra",
            "ReflexivePolytopesDB",
            "CohomCalg",
            "Topcom",
            "Triangulations",
            "InverseSystems",
            "IntegerEquivalences"
            },
        PackageImports => {
            "LLLBases"
            }
        )

export {
    -- Types defined here
    "CYPolytope", -- rename to CYReflexivePair?  How about CYPolytope?
    "CalabiYauInToric",
    "CYToolsCY3", -- we should have a superclass for CalabiYauInToric, CYToolsCY3
    "TopologicalDataOfCY3",

    -- CYPolytope, CalabiYauInToric
    "InteriorFacets",
    "ID",
    "cyPolytope",
    "dump",
    "label",
    
    "cyData",
    "makeCY",
    "calabiYau", -- versions might include:
       -- calabiYau(CYPolytope, Triangulation, Ring => RZ, Label => (a,i))
       -- The following are all taking the CY3 from a database:
       -- all take Ring as optional argument.
       -- calabiYau(String, CYPolytope, Label) -- from database, given Q
       -- calabiYau(String, Label) -- from database, uses Q from same database.
       -- calabiYau(String, Function, Label) -- from database, uses Q from same database.
       -- calabiYau(Database, CYPolytope, Label) -- from database, given Q.
       -- calabiYau(Database, Label) -- from database, uses Q from same database.
       -- calabiYau(Database, Function, Label) -- takes Q from Function, label.
       
    "basisIndices",
    "restrictTriangulation", -- restrict triangulation to each 2-face
    "picardRing",

    "isFavorable",

    -- Extra polyhedral facilities, for lattice points and faces of a Polyhedron
    -- how much of this shoiuld be exported??
    "vertexMatrix",
    "vertexList",
    "faceDimensionHash",
    "faceList",
    "dualFace",
    "minimalFace",
    "latticePointList",
    "latticePointHash",
    "interiorLatticePointList",
    "annotatedFaces",
    "automorphisms",

    -- current triangulation code
    "findAllFRSTs",
    "findAllCYs",
    "findAllConnectedStarFine",
    "findStarFineGraph",
    
    -- older triangulation code (still useful?)
    "Origin",
    "pointConfiguration",
    "regularStarTriangulation",
    "reflexiveToSimplicialToricVariety",
    "reflexiveToSimplicialToricVarietyCleanDegrees",
    "allZeros",
    "augmentWithOrigin",
    
    -- This set maybe should be included in NormalToricVarieties?
    "singularCones",
    "singularLocusInToric",
    "normalToricVarietyFromGLSM",

    -- Interfacing with Sage triangulations, triangulations from elsewhere
    "readSageTriangulations",
    "sortTriangulation",
    "matchNonZero",
    "applyPermutation",
    "checkFan",
    "sageTri",

    -- IntersectionNumbers
    "intersectionNumbers",
    "toricIntersectionNumbers",
    "c2",
    "intersectionForm",
    "cubicForm",
    "c2Form",
    "intersectionNumbersOfCY", -- possibly not for export

    -- Cones of curves
    "toricMoriCone",
    "toricMoriConeCap",
        
    -- gvInvariants
    "gvInvariants",
    "gvCone",
    "partitionGVConeByGV",
    "classifyExtremalCurves",
    "extremalRayGVs",

    -- remove these gvInvariant functions?
    "classifyExtremalCurve",
    "gvInvariantsAndCone",
    --    "gvRay",
    "findLinearMaps",

    
    -- Invariants
    "hubschInvariants",
    "PointCounter",
    "pointCounter",
    "pointCounts",
    "invariantsH11H12",
    "invariantContents",
    "hessianInvariants",
    "singularContents",
    "cubicConductorInvariants", -- rename?
    "cubicLinearConductorInvariants", -- rename?
    "singularContentsQuartic", -- rename?
    "aronhold",
    
    -- Topology
    "topologicalData",
    "isEquivalent",
    "invariants",
    "mapIsIsomorphism",
    "partitionByTopology",

    -- TopologySet's -- for determining equivalence
    "TopologySet",
    "topologySet",
    "representatives",
    "equivalences",
    "IgnoreSingles",
    "separateIfDifferent",
    
        
    -- CompleteIntersectionInToric's
    "completeIntersection",
    "CompleteIntersectionInToric",
        "Ambient",
        "CI",
        "Equations",
    "equations",
    "LineBundle",
    "lineBundle",
    "LinearForm",
    "Basis",
    

    -- Creating databases of polytopes (with precomputed data).
    "hodgeNumbers", -- of KSEntry: gives (h11, h12) from KSEntry.  Should be in ReflexivePolytopesDB?
    "createCYDatabase",
    "addToCYDatabase",
    "readCYDatabase",
    "readCYs",
    "readCYPolytopes",
    "combineCYDatabases",
    "processCYPolytopes",

    -- Cohomology
    "toricCohomologySetup",
    "cohomologyBasis",
    "CohomologySetup",
    "cohomologyMatrix",
    "cohomologyMatrixRank",
    "genericCohomologyMatrix",
    
    -- new attempt at faster cohomology of line bundles in toric varieties.
    -- want the Laurent monomials.
    "toricOrthants",
    "H1orthants",
    "normalDegrees",
    "normalDegree",
    "setOfColumns",
    "getFraction",
    "findTope",
    "cohomologyFractions",
    
    -- cohomology for complete intersections in torics
    --TODO clean up this interface!
    "cohomologyVector",  -- currently calls CohomCalg, but shouold take an optional argument?
    "cohomologyFromLES",
    "collectLineBundles",
    "cohomologyOmega1",
    "removeUnusedVariables",
    "example0912'3524",
    "exampleP111122'44",
    "cohom",
    "nextBundle",
    "basicCohomologies",
    "hodgeDiamond",
    "hodgeVector",
    
    -- CY hypersurface cohomology info from polyhedral information
    "subcomplex",
    "complexWithInteriorFacesRemoved",
    -- The following are the ones we care about
    "hodgeCY", -- debugging only?
    "hodgeOfCYToricDivisor",
    "hodgeOfCYToricDivisors",
    "h11OfCY",
    "h21OfCY",

    
    -- new formula
    "hodgeVectorViaTheorem", -- TODO: is likely not correct currently.
    "tentativeHodgeVector", -- deprecated
    
    -- Topology of a CY3-fold.
    -- input: reflexive poytope, and triangulation.
    -- Some functions to have here:
    -- type: CY3ToricHypersurface.
    -- h11 X
    -- h12, h21 X
    -- intersectionNumbers X
    -- topology X (returns an object of TopologicalDataOfCY3)
    
    -- Gopakumar-Vafa invariants
    -- computed using Andres' computeGV code in C++
    
    -- Flop chains, Mori cones

    "IntersectionNumbers",
    "MoriHilbertGens",
    "Automorphisms",
    "FilePrefix",
    "Executable",
    "Mori",
    "Count",
    "Hodge",
    "NTFE",
    "PicardRing"
    }

--- kludge to access parts of the 'Core'
hasAttribute = value Core#"private dictionary"#"hasAttribute";
getAttribute = value Core#"private dictionary"#"getAttribute";
ReverseDictionary = value Core#"private dictionary"#"ReverseDictionary";

------------------------------------
-- New types -----------------------
------------------------------------

TopologicalDataOfCY3 = new Type of List
  -- contains c2, cubic intersection form, h11, h12

CompleteIntersectionInToric = new Type of HashTable
  
CYPolytope = new Type of HashTable
CalabiYauInToric = new Type of HashTable
CYToolsCY3 = new Type of HashTable

LineBundle = new Type of HashTable

lineBundle = method()
equations = method()

dump = method(Options => true)
isFavorable = method();
h11OfCY = method() -- deprecate this
h21OfCY = method() -- deprecate this
findAllFRSTs = method()

-- Utility function
dotProduct = method()
dotProduct(List, List) := (v,w) -> (
    if #v =!= #w then error "expected vectors of the same length";
    sum for i from 0 to #v-1 list v#i * w#i
    )

findEquivalence(CalabiYauInToric, CalabiYauInToric) := (X1, X2) -> (
    findEquivalence({c2Form X1, cubicForm X1}, {c2Form X2, cubicForm X2})
    )

load (currentFileDirectory | "StringTorics/MyPolyhedra.m2")
load (currentFileDirectory | "StringTorics/CYPolytope.m2")
load (currentFileDirectory | "StringTorics/CalabiYauInToric.m2")
load (currentFileDirectory | "StringTorics/IntersectionNumbers.m2")
load (currentFileDirectory | "StringTorics/Invariants.m2")
load (currentFileDirectory | "StringTorics/Topology.m2")
load (currentFileDirectory | "StringTorics/ToricCompleteIntersections.m2") -- has some util code, but not much.  TODO: clean that up.
load (currentFileDirectory | "StringTorics/DatabaseCreation.m2")
load (currentFileDirectory | "StringTorics/Extras.m2")

  findAllConnectedStarFine = method()
  findAllConnectedStarFine Triangulation := (T) -> (
      stars := new MutableHashTable;
      stars#T = 0;
      starcount := 1;
      oldTODO := {T};
      radius := 0;
      while #oldTODO > 0 do (
          radius = radius + 1;
          newTODO := flatten for t1 in oldTODO list for n in neighbors t1 list (
              oldone := stars#t1;
              t := n#1;
              if isStar t and not stars#?t and isRegularTriangulation t then(
                  << "adding new triangulation " << starcount << " at radius " << radius << " with circuit " << n#0 << " from " << oldone << endl;
                  stars#t = starcount;
                  starcount = starcount + 1;
                  t
                  )
              else continue
              );
          oldTODO = newTODO;
          );
      keys stars
      )

  findStarFineGraph = method()
  findStarFineGraph Triangulation := (T) -> (
      -- this version returns the determined graph of the FRST's.
      -- 3 things are returned:
      --   1. a list of triangulations
      --   2. a hash table: for each triangulation index: key is a list of {tri#, affine circuit used to get to that}
      stars := new MutableHashTable;
      edges := new MutableList;      
      stars#T = 0;
      starcount := 1;
      oldTODO := {T};
      radius := 0;
      while #oldTODO > 0 do (
          radius = radius + 1;
          newTODO := flatten for t1 in oldTODO list for n in neighbors t1 list (
              newOneIsNew := false;
              oldone := stars#t1;
              t := n#1;
              alreadyThere := stars#?t;
              if not alreadyThere then (
                  if isStar t and isRegularTriangulation t then (
                      stars#t = starcount;
                      starcount = starcount + 1; 
                      newOneIsNew = true;
                      )
                  else continue
                  ); -- this is the case when we don't need to add an edge, nor place tri onto the newTODO list.
              -- at this point both t and oldone are good. So let's add an edge.
              edges#(#edges) = {oldone, stars#t, n#0};
              if newOneIsNew then t else continue
              );
          oldTODO = newTODO;
          );
      starsInv := hashTable for k in keys stars list stars#k => k;
      starsList := for i from 0 to starcount-1 list starsInv#i;
      (starsList, new List from edges)
      )

protect nextVar
protect nextFinalVar
---------------------------------------------------
singularLocusInToric = method()
singularLocusInToric(Ideal, Ideal) := (I, B) -> (
    -- I: an ideal in the Cox ring, where B is the irrelevant ideal.
    -- return: a list of (codim Sing locus, groebner basis) for each max cone
    -- TODO: this I think only works in the smooth case so far.  FIX THIS.
    -- find singular locus on each open set.
    c := codim I;
    Is := for b in B_* list (
      sub(I, for x in support b list x => 1)
      );
    for J in Is list (
        singJ := J + minors(c, jacobian J);
        (codim singJ, gens gb singJ)
        )
    )

---------------------------------------------------

allzeros0 = (F, topvar, pt, lo, hi) -> (
  -- F is a polynomial in a poly ring (e.g. QQ[a,b,c])
  -- topvar
  if topvar == -1 then (return if F == 0 then {pt} else {});
  R := ring F;
  flatten for v from lo to hi list (
      G := sub(F, R_topvar => v);
      pt1 := prepend(v, pt);
      allzeros0(G, topvar-1, pt1, lo, hi)
      )
  )
allzeros1 = (F, topvar, pt, lo, hi, f) -> (
  -- F is a polynomial in a poly ring (e.g. QQ[a,b,c])
  -- topvar
  if topvar == -1 then (if F == 0 then f pt; return null);
  R := ring F;
  for v from lo to hi do (
      G := sub(F, R_topvar => v);
      pt1 := prepend(v, pt);
      allzeros1(G, topvar-1, pt1, lo, hi, f)
      )
  )
allZeros = method()
allZeros(RingElement, Sequence) := (F, lohi) -> (
    (lo,hi) := lohi;
    allzeros0(F, numgens ring F - 1, {}, lo, hi)
    )
allZeros(RingElement, ZZ, Sequence) := (F, topvar, lohi) -> (
    (lo,hi) := lohi;
    allzeros0(F, topvar, {}, lo, hi)
    )
allZeros(RingElement, ZZ, Sequence, Function) := (F, topvar, lohi, f) -> (
    (lo,hi) := lohi;
    allzeros1(F, topvar, {}, lo, hi, f)
    )

reflexiveToSimplicialToricVariety = method(Options => {
        CoefficientRing => QQ, 
        Variable =>  getSymbol "x"
        })

reflexiveToSimplicialToricVariety Polyhedron := opts -> (P1) -> (
    -- P1 is a reflexive polytope in the M lattice.
    -- Creates a simplicial toric variety via a triangulation
    -- of the normal fan of P1, using all of lattice points of (polar P1).
    -- Returns: a normal toric variety, whch uses all
    -- lattice points in the triangulation.
    P2 := polar P1;
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    normalToricVariety(LP,tri,opts)
    )

-- Subroutines for reflexiveToSimplicialToricVarietyCleanDegrees
-- These might be generally useful too?
  findFirstUnitVectors = method()
  findFirstUnitVectors Matrix := List => (M) -> (
      -- returns the indices of the the columns which are unit vectors,
      -- if there are two or more columns corresponding to the same unit vector, take the first.
      -- the result is in increasing list of integers.
      e := entries transpose M;
      unitposition := (col) -> if all(col, a -> a >= 0) and sum col === 1 then position(col, a -> a === 1) else null;
      H := partition(c -> unitposition e_c, toList(0..numcols M-1));
      ks := select(keys H, k -> k =!= null);
      sort for k in ks list min(H#k)
      )

  -- extend the nice unit vectors to a list p of n columns (n = numrows M)
  -- so that det M_p = 1 or -1.
  -- compute the inverse A of this n x n matrix.
  -- then compute A^-1 * M, and return this list of columns.
  -- TODO: don't use subsets...
  findInvertibleSubmatrix = method(Options => {Limit => 10000})
  findInvertibleSubmatrix(Matrix, List) := List => opts -> (M, p) -> (
      -- M is a matrix over the integers
      -- p is a list of column indices of M (indices run from 0 to numcols M - 1)
      -- Find a list q containing p (if possible), of column indices,
      --   s.t. det(M_q) = 1 or -1.
      -- Return q, or null, if one cannot be found.
      S := set toList(0..numcols M - 1);
      others := sort toList(S - set p);
      needed := numrows M - #p;
      if binomial(#others, needed) > opts.Limit then return null;
      -- really: do some random choices of q containing p.
      trythese := subsets(others, needed);
      tried := 0;
      for q1 in trythese do (
          q := join(q1,p);
          tried = tried+1;
          if abs det(M_q) == 1 then (
              if debugLevel >= 1 then (
                  << "found suitable set of columns after " << tried << " step" 
                  << if tried == 1 then "" else "s" << endl;
                  );
              return sort q;
              );
          );
      << "tried all determinants: none had unit determinant" << endl;
      null
      )

reflexiveToSimplicialToricVarietyCleanDegrees = method(Options => options reflexiveToSimplicialToricVariety)
reflexiveToSimplicialToricVarietyCleanDegrees Polyhedron := Sequence => opts -> (P1) -> (
    -- returns (V, q) where
    -- V is essentially reflexiveToSimplicialToricVariety P1
    --   except that the degrees of the ring have been cleaned up
    -- AND 
    -- a list of column indices (or of rays) in ascending order, 
    -- whose degrees form a basis, in fact the identity matrix.
    -- null is returned if no such q can be found, or the algorithm
    -- doesn't find it.
    -- Note: we don't really need ring V0 here...
    P2 := polar P1;
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    D := transpose syz matrix transpose LP;
    print D;
    p := findFirstUnitVectors D;
    q := findInvertibleSubmatrix(D, p);
    if q === null then null
    else (
        A := D_q;
        D' := A^-1 * D;
        (normalToricVariety(LP, tri, WeilToClass => D', opts), q)
        )
    )

augmentWithOrigin = method()
augmentWithOrigin Matrix := (A) -> (
    -- A is a matrix over ZZ
    -- add column of 0's, then add in a first row of 1's.
    n := numColumns A;
    zeros := matrix {numRows A : {0}};
    ones := matrix {{(n+1) : 1}};
    ones || (A | zeros)
    )

augmentBack = (A) -> (
    -- A is a matrix over ZZ
    -- add in a last row of 1's.  Also add in one column of zeros
    n := numColumns A;
    zeros1 := matrix{numRows A : {0}};
    zeros0 := matrix {1+numRows A : {0}};
    ones := matrix {{n+1 : 1}};
    ((A | zeros1) || ones) --| zeros0
    )

readSageTriangulations = method()
readSageTriangulations String := (str) -> (
   --  "     [[array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 5]), array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 5]), array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 5]), array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])]] "
    s1 := replace("array", "", str);
    s2 := replace("\\(", "", s1);    
    s3 := replace("\\)", "", s2);
    s4 := replace("\\[", "{", s3);
    s5 := replace("\\]", "}", s4);
    value s5
    )

sortTriangulation = method()
sortTriangulation List := (T) -> sort for t in T list sort t

matchNonZero = method()
matchNonZero(Matrix, Matrix) := (A,B) -> (
    eA := entries transpose A;
    Ah := hashTable for i from 0 to #eA-1 list (
        if all(eA#i, x -> x == 0) then continue else eA#i => i
        );
    eB := entries transpose B;
    Bh := hashTable for i from 0 to #eB-1 list (
        if all(eB#i, x -> x == 0) then continue else eB#i => i
        );
    if set keys Ah =!= set keys Bh then (
        error "non-zero columns are different";
        );
    AtoB := for e in eA list if Bh#?e then Bh#e else -1;
    BtoA := for e in eB list if Ah#?e then Ah#e else -1;
    (AtoB, BtoA)
    )

applyPermutation = method()
applyPermutation(List, ZZ) := (P, i) -> if i >= 0 and i < #P then P#i else -1
applyPermutation(List, List) := (P, L) -> sort for f in L list applyPermutation(P,f)

---------------------------------------
-- Code for toric varieties packages --
---------------------------------------

  singularCones = method()
  singularCones(ZZ, NormalToricVariety) := (coneDim, X) -> (
      rys := transpose matrix rays X;
      which := orbits(X, coneDim);
      select(which, f -> (
              I := trim minors(#f, rys_f);
              I != 1))
      )
  singularCones(ZZ, NormalToricVariety, Thing) := (coneDim, X, notused) -> (
      -- This version also returns the multiplicity, for non smooth cones.
      rys := transpose matrix rays X;
      which := orbits(X, coneDim);
      for f in which list (
          I := trim minors(#f, rys_f);
          if I == 1 then continue;
          f => I_0
          )
      )

  QQCartierCoefficients = method ()
  QQCartierCoefficients ToricDivisor := List => D -> (
    X := variety D;
    rayMatrix := QQ ** matrix rays X;
    coeffs := QQ ** transpose (matrix {entries D});
    apply (max X, sigma -> coeffs^sigma // rayMatrix^sigma)
    );

  monomialsToLatticePoints = method()
  monomialsToLatticePoints(ToricDivisor, List) := List => (D, monoms) -> (
      -- given a list of monomials in the Cox ring, returns their corresponding lattice point in the M lattice.
      -- Maybe we do not need D?  Just the toric variety itself?
    X := variety D;    
    if not isProjective X then 
        error "--expected the underlying toric variety to be projective";
    degs := QQ ** transpose matrix rays X;
    deginverse := transpose(id_(QQ^(numrows degs)) // degs);
    coeff := matrix vector D;
    for mon in monoms list (
        flatten entries sub(deginverse * (transpose matrix{first exponents mon} - coeff), ZZ)
	)
    )

  normalToricVarietyFromGLSM = method(Options=>options normalToricVariety)
  normalToricVarietyFromGLSM(Matrix, List) := opts -> (GLSM, maxCones) -> (
      -- create rays from GLSM degrees (which will then be the output of 
      rays := entries syz GLSM;
      normalToricVariety(rays, maxCones, opts, WeilToClass => GLSM)
      )
  normalToricVarietyFromGLSM(Matrix, MonomialIdeal) := opts -> (GLSM, SRIdeal) -> (
      -- this switch to dual should be much simpler than this...
      n := numColumns GLSM; -- also should be the number of variables of the ring SRIdeal
      allvars := set toList(0..n-1);
      maxcones := (dual SRIdeal)_*/(f -> sort toList(allvars - (support f)/index//set));
      normalToricVarietyFromGLSM(GLSM, maxcones, opts)
      )

  -- This function isn't used anywhere??
  pointConfiguration = method(Options => {Homogenize=>false, Origin => false});
  pointConfiguration(ZZ,Polyhedron) := opts -> (maxdim, P2) -> (
      LP := select(latticePointList P2, lp -> dim(P2, minimalFace(P2, lp)) <= maxdim);
      A := transpose matrix LP;
      if opts.Origin then 
        A = A | matrix for i from 1 to numRows A list {0};
      if opts.Homogenize then
        A = A || matrix {for i from 1 to numColumns A list 0};
      A
      )

-- TODO: what should this be returning??  Probably a Triangulation object.
-- But with or without the cone vertex??
regularStarTriangulation = method()
regularStarTriangulation Polyhedron := (P2) -> (
    LP := drop(latticePointList P2, -1);
    A := transpose matrix LP;
    tri := max regularFineTriangulation A;
    -- Now, this is not a star triangulation...  But we think we can make it so in the
    -- following way.
    facetlist := (faceList(dim P2-1, P2))/(f -> latticePointList(P2,f));
    newtri := sort flatten for f in tri list (
        fac := select(facetlist, g -> #((set f) * (set g)) == dim P2);
        for g in fac list sort toList ((set f) * (set g))
        );
    (LP,newtri)
    )

regularStarTriangulation(ZZ,Polyhedron) := (maxdim, P2) -> (
    LP := select(latticePointList P2, lp -> dim(P2, minimalFace(P2, lp)) <= maxdim);
    A := transpose matrix LP;
    A = A | matrix for i from 1 to numRows A list {0};
    tri := topcomRegularFineTriangulation A;
    -- Now, this is not a star triangulation...  But we think we can make it so in the
    -- following way.
    H := latticePointHash P2;
    fullset := set for lp in LP list H#lp;
    facetlist := faceList(dim P2-1, P2);
    facetlist = facetlist/(f -> (
            sort toList(set latticePointList(P2,f) * fullset)
            ));
    newtri := sort flatten for f in tri list (
        fac := select(facetlist, g -> #((set f) * (set g)) == dim P2);
        for g in fac list sort toList ((set f) * (set g))
        );
    (LP,newtri)
    )

  sageTri = "     [[array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), 
array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), 
array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), 
array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), 
array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 3]), 
array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), 
array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), 
  array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]),   
  array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), 
  array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
  array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
  array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), 
  array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), 
  array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), 
  array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), 
  array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), 
  array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
  array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
  array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), 
  array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), 
  array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), 
  array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), 
  array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), 
  array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 3]), 
  array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), 
  array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), 
  array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
  array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], 
[array([0, 1, 2, 5]), array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), 
    array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), 
    array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), 
    array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), 
    array([1, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), 
    array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), 
    array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), 
    array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
    array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), 
    array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), 
    array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), 
    array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), 
    array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 5]), array([0, 1, 4, 5]), 
    array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), 
    array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), 
    array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), 
    array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 5]), 
    array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), 
    array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), 
    array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
    array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], 
[array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), 
    array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), 
    array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), 
    array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), 
    array([2, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), 
    array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), 
    array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), 
    array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
    array([1, 2, 4, 6]), array([2, 4, 5, 6])]] "

load (currentFileDirectory | "StringTorics/LineBundleCohomology.m2")

-----------------------
-- Examples -----------
-----------------------
example0912'3524 = () -> (value /// () -> (
    R := QQ[v1, v2, v3, v4, v5, v6, v1s, v7, v8, v9, v10];
    SR := {{v3,v9},{v5,v9},{v7,v10},{v1,v2,v3},
          {v4,v1s,v8},{v4,v7,v8},{v4,v8,v9},
          {v5,v6,v1s},{v5,v6,v10},{v1,v2,v6,v1s}};
    SR = monomialIdeal(SR/product);
    GLSM := transpose matrix {{3, 3, 3, 3, 0}, {2, 2, 2, 2, 0},
      {1, 0, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0},
      {0, 1, 0, 0, 0}, {0, 1, 1, 0, 0}, {0, 0, 1, 0, 1},
      {0, 0, 1, 0, 0}, {0,-1,-1, 1,-1}, {0, 0, 0, 0, 1}};
    normalToricVarietyFromGLSM(GLSM, SR)
    )
  ///
  )

exampleP11222'8 = () ->  (value /// () -> (
        -- P^4_{11222}[8]
        -- Actually: we desingularize the ZZ/2 - singularity of this
        rays := {{-1,-2,-2,-2}, 
                 {1,0,0,0},
                 {0,1,0,0},
                 {0,0,1,0},
                 {0,0,0,1},
                 {0,-1,-1,-1}};
        maxcones := {
            {0, 2, 3, 4}, 
            {0, 2, 3, 5}, 
            {0, 2, 4, 5}, 
            {0, 3, 4, 5}, 
            {1, 2, 3, 4}, 
            {1, 2, 3, 5}, 
            {1, 2, 4, 5}, 
            {1, 3, 4, 5}};
        GLSM := transpose matrix {
            {1,0}, {1,0}, {2,1}, {2,1}, {2,1}, {0,1}};
        normalToricVariety(rays, maxcones, WeilToClass => GLSM)
        )
    ///
    )

exampleP111122'44 = () -> (value /// () -> (
        rays := {{-1,-1,-1,-2,-2}, 
                 {1,0,0,0,0},
                 {0,1,0,0,0},
                 {0,0,1,0,0},
                 {0,0,0,1,0},
                 {0,0,0,0,1}};
        maxcones := {
            {0, 1, 2, 3, 4},
            {0, 1, 2, 3, 5},
            {0, 1, 2, 4, 5},
            {0, 1, 3, 4, 5},
            {0, 2, 3, 4, 5},
            {1, 2, 3, 4, 5}};
        GLSM := transpose matrix {
            {1}, {1}, {1}, {1}, {2}, {2}};
        normalToricVariety(rays, maxcones, WeilToClass => GLSM)
        )
    ///
    )

 findAllFRSTs NormalToricVariety := (V) -> findAllFRSTs transpose matrix rays V
 findAllFRSTs Matrix := List => (A) -> (
     A1 := A | map(target A, (ring source A)^1, 0);
     Ts := allTriangulations(A1, Fine => true, RegularOnly => true);
     if #Ts === 0 or #Ts#0 == 0 then (
         count := 0;
         while count < 10 and (#Ts === 0 or #Ts#0 == 0) do (
             Ts = allTriangulations(A1, Fine => true, RegularOnly => true);
             count = count + 1;
             );
         --if #Ts == 0 then error "no triangulation could be found";
         << "WARNING: TOPCOM failed to find triangulations, then found them after " << 
           count << " attempt(s)" << endl;
         );
     --<< "Ts = (before selection): " << netList Ts << endl;
     Ts1 := select(Ts, isStar);
     --<< "Ts1 = (after): " << netList Ts1 << endl;
     assert all(Ts1, tri -> all(max tri, s -> s#-1 == numcols A));
     Ts1/(t -> (entries transpose A, (max t)/(s -> drop(s, -1))))
     )
 findAllFRSTs Polyhedron := List => (P) -> (
     L := latticePointList P;
     assert all(L#-1, a -> a == 0);
     L = drop(L, -1);
     A := transpose matrix L;
     findAllFRSTs A
     )


-- Being rewritten 22 Aug 2023.
-- findAllCYs = method(Options => {Ring => null}) -- opts.Ring: ZZ[h11 variables].
-- findAllCYs CYPolytope := List => opts -> Q -> (
--     Ts := findAllFRSTs Q;
--     RZ := if opts#Ring === null then (
--         a := getSymbol "a";
--         h11 := hh^(1,1) Q;
--         ZZ[a_1 .. a_h11]
--         )
--     else (
--         opts#Ring
--         );
--     for i from 0 to #Ts - 1 list cyData(Q, Ts#i, ID => i, Ring => RZ)
--     )


-- keys: id, cypolytopedata, triangulation, cache.  The id is what? (id of polytope, which triangulation)
--  write date: for cypolytopedata, just writes the id.
--  read data: given id, need to be able to get at which polytope it is.
--    maybe a table with id => CYPolytope, or a function which takes an integer and returns 
--    the CYPolytope object to use, with this id.
--  construct one from a CYPolytope, id, triangulation.
--  what is in the cache?
--    ambient toric
--    CYInToric?
--    abstract toric variety (depends on base)
--    abstract variety for CY3 (depends on base)
--    intersectionNumbers
--    cubicForm (string or polynomial? or intersection numbers only?)
--    c2Form (string or list or polynomial?)
--    mori cone info?
--    gv invariants?


----------------------------------------------------------------


----------------------------------------------------------------
-- FRST Triangulations (Fine, regular, star triangulations) ----
----------------------------------------------------------------
  
reflexiveToSimplicialToricVariety Polyhedron := opts -> (P1) -> (
    -- P1 is a reflexive polytope in the M lattice.
    -- Creates a simplicial toric variety via a triangulation
    -- of the normal fan of P1, using all of lattice points of (polar P1).
    -- Returns: a normal toric variety, whch uses all
    -- lattice points in the triangulation.
    P2 := polar P1;
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    normalToricVariety(LP,tri,opts)
    )


load (currentFileDirectory | "StringTorics/GVInvariants.m2")

-- This file refers to many of the method names defined earlier, applied to CYToolsCY3
load (currentFileDirectory | "StringTorics/CYTools.m2")

beginDocumentation()

-- Problems
-- . how to determine favorable?
-- . which toric variety do we want?
-- . want a smooth one
-- . want to compute h^11(X) using cohomcalg, but
--   examples can be too big
-- . triangulations can be too big
-- . what else can be too big?

load (currentFileDirectory | "StringTorics/doc.m2")
load (currentFileDirectory | "StringTorics/test.m2")

end--

restart
  uninstallAllPackages()

restart
  installPackage "IntegerEquivalences" -- works, lots of warnings
  installPackage "DanilovKhovanskii"
  installPackage "StringTorics"


  check IntegerEquivalences -- 8 checks, finishes to completion.
  check DanilovKhovanskii -- 10 checks, finishes, 3 take some time
  check StringTorics -- 35 tests, finishes to completion.  3 tests take > 10 sec.

restart
needsPackage "StringTorics"
path = append(path, "......")
needsPackage("StringTorics", FileName => "/Users/.../StringTorics.m2")
restart
uninstallPackage "StringTorics"
restart
installPackage "StringTorics"
viewHelp oo

restart
needsPackage "StringTorics"
check oo -- all tests currently check.

-- Generation of some examples
L = kreuzerSkarke(20, Limit => 100, Access=>"wget")
L = kreuzerSkarke(15, Limit => 100)
L = kreuzerSkarke(13, Limit => 100)
L = kreuzerSkarke(10, Limit => 100)
A = matrix L_50
P = convexHull A
isReflexive P
X = reflexiveToSimplicialToricVariety P
  max X
  # rays X
  sr = dual monomialIdeal X
  for i from 0 to #L-1 list (
      elapsedTime X = reflexiveToSimplicialToricVariety convexHull matrix L_i;
      if # rays X > 64 then continue;
      elapsedTime S = try ring X else null;
      if S === null then continue;
      sr = elapsedTime dual monomialIdeal X;
      if numgens sr > 64 then continue;
      i
      )
  Xs = for i from 0 to #L-1 list (
      elapsedTime X = reflexiveToSimplicialToricVariety convexHull matrix L_i
      )

  tally apply(Xs, X -> # rays X)
  select(Xs, X -> elapsedTime try (ring X; true) else false) -- h11=13, h11=15: all ok here...
  tally apply(Xs, X -> numgens elapsedTime dual monomialIdeal X)
  tally apply(Xs, X -> # rays X)
  positions(Xs, X -> numgens elapsedTime dual monomialIdeal X < 64)
  max Xs_28
  debug CohomCalg
  print toCohomCalg Xs_28  -- h11 = 15
  
  -- h11 = 15:
  positions(Xs, X -> numgens elapsedTime dual monomialIdeal X == 39) -- h11=15 indices: 51, 66.
  max Xs_51      
  rays Xs_51
  print toCohomCalg Xs_51

  -- h11 = 15:
  positions(Xs, X -> numgens elapsedTime dual monomialIdeal X == 24) -- h11=10 indices: 6.
  max Xs_6
  rays Xs_6
  print toCohomCalg Xs_6
    

str = getKreuzerSkarke(4,100,Limit=>10)
time polytopes = parseKS str;
#polytopes
time for p in polytopes list (p_0, matrixFromString p_1);


str = getKreuzerSkarke(50,100,Limit=>4000);
str = getKreuzerSkarke(50,100,Limit=>0)
str = getKreuzerSkarke(4,4,Limit=>10)
str = getKreuzerSkarke(4,8,Limit=>10)
for i from 1 to 20 list parseKS getKreuzerSkarke(4,i,Limit=>10)
for i from 21 to 30 list parseKS getKreuzerSkarke(4,i,Limit=>10)
for i from 31 to 40 list parseKS getKreuzerSkarke(4,i,Limit=>10)
join oo
flatten oo
netList oo

-- Experiments for Batyrev formula
restart
needsPackage "StringTorics"
getKreuzerSkarke(3,57)
polytopes = parseKS oo;
netList polytopes
polystr = polytopes_2
A = matrixFromString polystr_1
P = convexHull A
P2 = polar P

-- Part 1:
# latticePoints(P2) == 8

-- Part2: facets of P2
faces1 = faces(1, P2)

F = faces1_0
vertices F
latticePoints F

for f in faces1 list ((# latticePoints f) - numColumns vertices F)
F = faces1_3
hyperplanes F
(first (hyperplanes F)) * vertices F
H = for f in faces(1,F) list hyperplanes f
latticePoints F
for p in oo list H_0_0 * p

-- Need a function: given a lattice point, and a polytope, is the point in the polytope?
  -- Is it on the boundary?
L = latticePoints F  
interiorLatticePoints F
for p in L list contains(F,p)
for f in faces(1,P2) list interiorLatticePoints f
for f in faces(2,P2) list interiorLatticePoints f
F22 = faces(2,P2)
first oo
polar oo
vertices oo

F2 = faces(2,P)
#F22
#F2
faces(4,P)
faces(1,P2)
# faces(2,P2)
# faces(3,P)
faces(3,P2)
faces(2,P)

-- Ben's code for morigenerators, translated into M2.
moriGenerators = method()
moriGenerators(List, List) := (rays, cones) -> (
    -- rays: a list of points in ZZ^n
    -- cones: a list of lists of size n, each element is a subset of 0..n-1,
    --   representing the simplicial cones of a fan with rays 'rays'
    n := #(rays#0);
    cones := cones/set;
    -- assert: all rays have length n, and all cones have n as well
    walls := select(subsets(n,2), x -> # toList (cones#(x_0) * cones#(x_1)) == n-1);
    commons := for x in walls list toList (cones#(x_0) + cones#(x_1));
    walls
    )
rays25 = {{2, -1, -1, -1}, {-1, 0, -1, -1}, {-1, 1, -1, -1}, {-1, 0, 
    4, -1}, {-1, 1, 0, -1}, {-1, 0, -1, 0}, {-1, 1, -1, 4}, {-1, 0, 4,
     0}, {-1, 1, 0, 4}, {-1, 1, 0, 3}, {-1, 1, -1, 3}, {-1, 1, 0, 
    2}, {-1, 1, -1, 2}, {-1, 1, 0, 1}, {-1, 1, -1, 1}, {-1, 1, 0, 
    0}, {-1, 1, -1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {0, 0, -1, 1}, {0,
     0, 1, 0}, {0, 0, -1, 0}, {0, 0, 1, -1}, {0, 0, 0, -1}, {0, 
    0, -1, -1}, {-1, 0, 3, 0}, {-1, 0, 2, 0}, {-1, 0, 1, 0}, {-1, 0, 
    0, 0}, {-1, 0, 3, -1}, {-1, 0, 2, -1}, {-1, 0, 1, -1}, {-1, 0, 
    0, -1}};
cones25 = {{0, 7, 9, 11}, {0, 6, 9, 11}, {6, 7, 9, 11}, {0, 4, 10, 
    12}, {0, 5, 10, 12}, {4, 5, 10, 12}, {0, 4, 7, 15}, {0, 4, 6, 
    15}, {4, 6, 7, 15}, {0, 2, 4, 16}, {0, 2, 5, 16}, {2, 4, 5, 
    16}, {6, 7, 8, 17}, {7, 8, 9, 17}, {6, 7, 8, 9}, {0, 6, 8, 
    17}, {0, 8, 9, 17}, {0, 6, 8, 9}, {0, 6, 7, 17}, {0, 7, 9, 
    17}, {0, 6, 7, 18}, {3, 4, 5, 6}, {5, 6, 10, 19}, {4, 5, 6, 
    10}, {5, 6, 18, 19}, {0, 4, 6, 10}, {0, 6, 10, 19}, {0, 6, 18, 
    19}, {0, 5, 10, 19}, {0, 5, 18, 19}, {7, 11, 13, 20}, {6, 7, 11, 
    13}, {7, 13, 15, 20}, {6, 7, 13, 15}, {0, 11, 13, 20}, {0, 6, 11, 
    13}, {0, 13, 15, 20}, {0, 6, 13, 15}, {0, 7, 11, 20}, {0, 7, 15, 
    20}, {5, 12, 14, 21}, {4, 5, 12, 14}, {5, 14, 16, 21}, {4, 5, 14, 
    16}, {0, 4, 12, 14}, {0, 12, 14, 21}, {0, 4, 14, 16}, {0, 14, 16, 
    21}, {0, 5, 12, 21}, {0, 5, 16, 21}, {2, 3, 4, 22}, {2, 3, 4, 
    5}, {3, 4, 7, 22}, {3, 4, 6, 7}, {0, 2, 4, 22}, {0, 4, 7, 22}, {0,
     2, 3, 22}, {0, 3, 7, 22}, {0, 2, 3, 23}, {1, 2, 5, 24}, {1, 2, 
    23, 24}, {0, 2, 5, 24}, {0, 2, 23, 24}, {0, 1, 5, 24}, {0, 1, 23, 
    24}, {6, 7, 18, 25}, {3, 6, 7, 25}, {0, 7, 18, 25}, {0, 3, 7, 
    25}, {6, 18, 25, 26}, {3, 6, 25, 26}, {0, 18, 25, 26}, {0, 3, 25, 
    26}, {6, 18, 26, 27}, {3, 6, 26, 27}, {0, 18, 26, 27}, {0, 3, 26, 
    27}, {6, 18, 27, 28}, {3, 6, 27, 28}, {0, 18, 27, 28}, {0, 3, 27, 
    28}, {5, 6, 18, 28}, {3, 5, 6, 28}, {0, 5, 18, 28}, {0, 3, 5, 
    28}, {2, 3, 23, 29}, {2, 3, 5, 29}, {0, 3, 5, 29}, {0, 3, 23, 
    29}, {2, 23, 29, 30}, {2, 5, 29, 30}, {0, 5, 29, 30}, {0, 23, 29, 
    30}, {2, 23, 30, 31}, {2, 5, 30, 31}, {0, 5, 30, 31}, {0, 23, 30, 
    31}, {2, 23, 31, 32}, {2, 5, 31, 32}, {0, 5, 31, 32}, {0, 23, 31, 
    32}, {1, 2, 23, 32}, {1, 2, 5, 32}, {0, 1, 5, 32}, {0, 1, 23, 
    32}};

-- Investigate cohomology cones
-- fano(4,10):
B = ideal V
S = (coefficientRing ring B)[gens ring B, DegreeRank=>numgens ring B]
B = sub(B,S)
Ext^2(comodule B, S)
  -- get x_0*x_1*x_6, x_0*x_5*x_6
  -- remember: allow these to be as negative as desired, rest of the variables must be positive.
  -- Still: get a cone in RR^n (n is 7 here).
  degs = degrees ring V
  for i from 0 to numgens ring V list if (x_0*x_1*x_6)%x_i == 0 then - degs_i else degs_i
  

---------------------------------------------------
-- remove code below this line (18 May 2019) !!! --
---------------------------------------------------

-- This changes the function from ReflexivePolytopesDB
parseKS String := (str) -> (
    -- result is a List of pairs of strings.
    locA := regex("<b>Result:</b>\n", str);
    locB := regex("#NF", str);
    if locB === null then locB = regex("Exceeded", str);
    --if locA === null or locB === null then error "data not in correct Kreuzer-Skarke format";
    firstloc := if locA === null then 0 else locA#0#0 + locA#0#1;
    lastloc := if locB === null then #str else locB#0#0;
    --firstloc := locA#0#0 + locA#0#1;
    --lastloc := locB#0#0;
    cys := substring(firstloc, lastloc-firstloc, str);
    cys = lines cys;
    cys = select(cys, s -> #s > 0);
    starts := positions(cys, s -> s#0 != " ");
    starts = append(starts, #cys);
    for i from 0 to #starts-2 list (
        cys_(starts#i) | " id:" | i | "\n" |  demark("\n", cys_{starts#i+1 .. starts#(i+1)-1})
        )
    )

-- expect that str is something like:
str = ///4 5  M:53 5 N:9 5 H:3,43 [-80] id:0
   1   0   2   4 -10
   0   1   3   5  -9
   0   0   4   0  -4
   0   0   0   8  -8
   ///

matrixFromKS = method()
matrixFromKS String := (str) -> (
    matrixFromString concatenate between("\n", drop(lines str, 1))
    )

matrixFromKSEntry = method()
matrixFromKSEntry String := (str) -> (
    matrixFromString concatenate between("\n", drop(lines str, 1))
    )

///
  assert(matrixFromKS str == matrix {{1, 0, 2, 4, -10}, {0, 1, 3, 5, -9}, {0, 0, 4, 0, -4}, {0, 0, 0, 8, -8}})
///
  
  
