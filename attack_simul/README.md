# Attack simulation

We have implemented serialized versions of some of the smoothness
algorithm of section 4.  Each attack picks a random value x, and
decompses it as a product (or fraction) of smooth numbers.  It does not
include the final step of the algorithm to compute the root.

The code requires clang 16 (with support for arbitrary precision integers).

When run for the first time, it will precompute a set of primes and
some associated data, used in the algorithm.

## Parallel smoothness test

The file `parallel_smooth.c` implement the smoothness algorithm of section 4.5,
with parallel smoothness test.

## Parallel smoothness test with rational reconstruction

The files `rational-ref.c` and `rational-opt.c` implement the smoothness
algorithm of section 4.6, with with parallel smoothness test and
rational reconstruction. `rational-ref.c` is a relatively direct
implementation, while `rational-opt.c` uses aditional precomputed tables
to compute modular reductions and inverses faster.
