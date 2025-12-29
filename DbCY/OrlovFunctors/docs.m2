doc ///
Node
  Key
    OrlovFunctors
  Headline
    functors between singularity and derived categories for Calabi-Yau varieties
  Description
    Text
      Blah
  --   Example
  --   Code
  -- Contributors
  -- References
  -- Caveat
  -- SeeAlso
  Subnodes
    singularityToDerived
    subcomplexByDegrees
    orlovTruncateGeqDualize
    supTruncate

Node
  Key
    singularityToDerived
   (singularityToDerived, ZZ, Module)
   (singularityToDerived, ZZ, Matrix)
   (singularityToDerived, ZZ, Complex)
   --(singularityToDerived, ZZ, ComplexMap) ??
  Headline
    apply Orlov's functor from the singularity category to the derived category
  -- Usage
  -- Inputs
  -- Outputs
  Description
    -- Text
    Example
      R = ZZ/101[x_0] / ideal(x_0^3)
      M = coker vars R
      singularityToDerived(1, M, LengthLimit => 1)
  --   Code
  -- ExampleFiles
  -- Contributors
  -- References
  -- Caveat
  -- SeeAlso

Node
  Key
    subcomplexByDegrees
   (subcomplexByDegrees, Complex,    Sequence)
   (subcomplexByDegrees, ComplexMap, Sequence)
  Headline
    produce the subcomplex of a graded free complex with degrees in a given range
  Usage
    subcomplexByDegrees(F, (mindeg, maxdeg))
  Inputs
    F:Complex -- graded free complex (or morphism of them)
    "(mindeg, maxdeg)":
      a degree range
  Outputs
    :Complex -- a subcomplex (or morphism of subcomplexes)
  Description
    Text
      This is equivalent to using @TT "submatrixByDegrees(f, (mindeg, maxdeg), (mindeg, maxdeg))"@ on each differential.
    -- Example
  SeeAlso
    submatrixByDegrees
  Subnodes
    orlovTruncateLess
    orlovTruncateGeq

Node
  Key
    orlovTruncateLess
   (orlovTruncateLess, ZZ, Complex)
   (orlovTruncateLess, ZZ, ComplexMap)
  Headline
    produce the subcomplex of a graded free complex with generators less than a given degree
  Usage
    orlovTruncateLess(i, F)
  Inputs
    i:ZZ -- upper bound degree
    F:Complex -- graded free complex (or morphism of them)
  Outputs
    :Complex
      a subcomplex (or morphism of subcomplexes) given by summands of the form
      $R(j)$ with $j>-i$ (so $R(j)$ is generated in degree $<i$)
  Description
    Text
      This is equivalent to using @TT "subcomplexByDegrees(F, (, i-1))"@,
      which in turn is equivalent to @TT "submatrixByDegrees(f, (, i-1), (, i-1))"@ on each differential.
    -- Example
  SeeAlso
    submatrixByDegrees
    orlovTruncateGeq

Node
  Key
    orlovTruncateGeq
   (orlovTruncateGeq, ZZ, Complex)
   (orlovTruncateGeq, ZZ, ComplexMap)
  Headline
    produce the subcomplex of a graded free complex with generators at or above a given degree
  Usage
    orlovTruncateGeq(i, F)
  Inputs
    i:ZZ -- lower bound degree
    F:Complex -- graded free complex (or morphism of them)
  Outputs
    :Complex
      a subcomplex (or morphism of subcomplexes) given by summands of the form
      $R(j)$ with $j\leq -i$ (so $R(j)$ is generated in degree $\geq i$)
  Description
    Text
      This is equivalent to using @TT "subcomplexByDegrees(F, (i, ))"@,
      which in turn is equivalent to @TT "submatrixByDegrees(f, (i, ), (i, ))"@ on each differential.
    -- Example
  SeeAlso
    submatrixByDegrees
    orlovTruncateLess
    orlovTruncateGeqDualize

Node
  Key
    orlovTruncateGeqDualize
   (orlovTruncateGeqDualize, ZZ, Module)
   (orlovTruncateGeqDualize, ZZ, Matrix)
   (orlovTruncateGeqDualize, ZZ, Complex)
   --(orlovTruncateGeqDualize, ZZ, ComplexMap) ??
  -- Headline
  -- Usage
  -- Inputs
  -- Outputs
  -- Description
  --   Text
  --   Example
  --   Code
  -- ExampleFiles
  -- Contributors
  -- References
  -- Caveat
  -- SeeAlso

Node
  Key
    supTruncate
   (supTruncate, ZZ, Module)
   --(supTruncate, ZZ, Matrix) ??
   (supTruncate, ZZ, Complex)
   --(supTruncate, ZZ, ComplexMap) ??
  -- Headline
  -- Usage
  -- Inputs
  -- Outputs
  -- Description
  --   Text
  --   Example
  --   Code
  -- ExampleFiles
  -- Contributors
  -- References
  -- Caveat
  -- SeeAlso
///

end--

uninstallPackage "OrlovFunctors"
restart
installPackage "OrlovFunctors"
viewHelp OrlovFunctors
