doc ///
Node
  Key
    OrlovFunctors
  Headline
    Orlov's functors between singularity and derived categories for Calabi-Yau varieties

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
  -- Description
  --   Text
  --   Example
  SeeAlso
    orlovTruncateLess

Node
  Key
    singularityToDerived

Node
  Key
    supTruncate
///
