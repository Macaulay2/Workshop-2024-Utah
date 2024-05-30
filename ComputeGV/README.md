# GV Invariants

This folder contains the C++ code used to compute GV invariants. The following dependencies are needed to build it:
- [GMP](https://gmplib.org/)
- [MPFR](https://www.mpfr.org/)
- [MPFR C++](http://www.holoborodko.com/pavel/mpfr/) the mpreal.h header is included in this repo, so it is not necessary to install anything.
- [Clang](https://clang.llvm.org/) ([gcc](https://gcc.gnu.org/) can be used instead, but Clang produces slightly a faster binary)

To install the dependencies on Ubuntu/Debian-based distros you can use

```bash
apt-get install libgmp-dev libmpfr-dev clang
```

To compile the code you use

```bash
clang -c -std=c++17 -march=native -O3 -fPIC computeGV.cpp
#g++ -c -std=c++17 -march=native -O3 computeGV.cpp   # If using gcc
g++ -o computeGV computeGV.o -lmpfr -lgmp -pthread   # I'm not sure why linking with Clang fails
```

An example input file is provided. The input data is organized as follows.
- The list of curves. These are points in the Mori cone in the chosen basis. Only the generating curves must be given, as it will check if there are any missing curves below the specified maximum degree.
- The set of curves whose past light cone will be used. Usually one would just input an empty list []. This is useful if one wants to compute the GV invariant of a specific curve with the minimum computation effort.
- The grading vector. This is the vector whose dot product with curves determines the degree of the curves.
- The GLSM charge matrix. This is an h11 by h11+4 matrix whose rows specify the basis.
- The nef-partition of the polytope. If an empty list is given then it assumes an anticanonical hypersurface.
- The nonzero intersection numbers in the format [[i,j,k,K_ijk],...]. For three-folds these are simply triple-itersection numbers between divisors, but for higher-dimensional CYs the first index labels an element of H_{2,2}.
- A vector containing four entries. The first one is the maximum degree, and the second one is the number of decimal digits of precision that will be used in the computations. If the maximum degree is set to a negative number, it will be inferred from the input list of curves. The third integer is the mode: (0) normal, (1) Hilbert and (2) verbatim. The last parameter is the free RAM threshold below which the computation will be terminated.

Note that the input vectors can be specified with square brackets, curly brackets or parenthesis.

The compiled binary takes the data from the standard input, so it can be used as follows.

```bash
./computeGV < example_input_3fold_hyper.txt
```

The computed GV invariants are printed to the standard output, whereas the process information is printed to the standard error. So the GV invariants can be saved as follows.

```bash
./computeGV < example_input_3fold_hyper.txt 1>results.txt
```
