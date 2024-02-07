#! /usr/bin/env python3

import bach_random_factored_numbers
import gmpy2
import random
import math

if __name__ == "__main__":
    import sys

    num = 2**16
    
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print("Usage: python sommth_prob.py bits B B' [num]")
        print(" bits: Number of bits (e.g. 512 generates numbers between 2^511 and 2^512)")
        print(" B: smoothness bound")
        print(" B': extra prime")
        print(" [num]: sample size (2^num)")
        print ()
        sys.exit(1)

    bits = int(sys.argv[1])
    B = 2**int(sys.argv[2])
    BB = 2**int(sys.argv[3])
    if len(sys.argv) == 5:
        num = 2**int(sys.argv[4])
    n = gmpy2.mpz(2**bits)

    print(f"Testing {num} {bits}-bit numbers for ({B},{BB})-almost smoothness")
    
    bach_random_factored_numbers.random_state = gmpy2.random_state(random.randint(1, 10**10))
    count = 0
    for i in range(num):
        (p, pi) = bach_random_factored_numbers.process_r(n)
        for x in pi:
            if x < B:
                assert((p % x) == 0)
                p //= x
        if p < BB:
            count += 1

    print(f"{count} almost smooth numbers found")
    if count:
        print(f"Experimental proba: {count}/{num} = 2^{math.log2(count/num)}")
    else:
        print(f"Experimental proba < 2^{math.log2(1/num)}")
