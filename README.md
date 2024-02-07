# Smoothness simulation

This code implements a serialized version of the smoothness algorithm
of section 4.6.  It picks a random value x, and decompses it as a
fraction of two smooth numbers.  (It does not include the final step
of the algorithm to compute the root).

The code requires clang 16 (with support for arbitrary precision integers).

When run for the first time, it will precompute a set of primes and
some associated data, used in the algorthm.

The are two versions of the code: smooth-ref.c is a relatively direct
implementation, while smooth-opt.c uses precomputed tables to compute
modular reductions and inverses faster.
