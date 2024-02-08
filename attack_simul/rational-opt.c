#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#include <unistd.h>
#include <sys/random.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "libdivide.h"

typedef unsigned _BitInt(256) word;
typedef signed _BitInt(256) signed_word;
typedef unsigned _BitInt(512) double_word;

#define PSIZE 256
// #define PSIZE 128

uint32_t MOD(word x, uint32_t n, uint32_t reduce32[8]) {
  uint64_t t = (uint32_t)x;
  for (int i=1; i<8; i++) {
    t += (uint64_t)reduce32[i]*(uint32_t)(x>>(32*i));
  }
  return t%n;
}

uint32_t SMOD(signed_word x, uint64_t n, uint32_t reduce32[8]) {
  if (x>0) {
    return MOD(x, n, reduce32);
  } else {
    uint32_t r = MOD(-x, n, reduce32);
    if (r == 0)
      return 0;
    else
      return n-r;
  }
}

  // #define SMOD(x,n) ((x)%(n)>0? (x)%(n) : n+((x)%(n)))// Safe mod with unsigned x

void print_word(word x) {
  char buffer[100];
  buffer[99] = 0;
  int i = 99;
  while (x) {
    assert(i>1);
    buffer[--i] = '0' + (x%10);
    x /= 10;
  }
  printf("%s", buffer+i);
}

float log_word(word x) {
  int z = 0;
  while (x >= ((word)1<<64)) {
    z += 32;
    x>>=32;
  }
  return z + log2((uint64_t)x);
}

word abs_word(signed_word x) {
  if (x<0)
    return -x;
  else
    return x;
}

typedef struct {
  uint32_t B;
  uint32_t N;
  struct {
    uint32_t factor;
    uint32_t prime;
    uint32_t lg;
    uint32_t phi;
    uint32_t reduce32[8];
    struct libdivide_u64_t fast_d;
  } data[0];
} prime_data;

// Fixed point arithmetic
#define FIXEDP(x) ((int32_t)((x)*(1<<16)))

uint64_t mod(uint64_t x, uint32_t p, struct libdivide_u64_t *fast_d) {
  uint32_t d = libdivide_u64_do(x, fast_d);
  return x - p*d;
}

uint32_t powermod(uint32_t x, uint32_t e, uint32_t p, struct libdivide_u64_t *fast_d) {
  if (e==0) {
    return 1;
  } else if (e==1) {
    return x;
  }

  uint32_t y = 1;
  while (e > 1) {
    if (e%2) {
      y = mod(((uint64_t)x*y), p, fast_d);
    }
    x = mod(((uint64_t)x*x), p, fast_d);
    e /= 2;
  }
  return mod(((uint64_t)x*y), p, fast_d);
}

prime_data *precompute(uint32_t B) {
  uint32_t *sieve = calloc(B/32, 4); // bit array
  assert(sieve);
#define   SET(i) sieve[(i)/32] |= ((uint32_t)1<<((i)%32))
#define   GET(i) ((sieve[(i)/32]>>((i)%32)) & 1)
#define UNSET(i) sieve[(i)/32] &= ~((uint32_t)1<<((i)%32))

  // Sieve prime numbers
  uint32_t N = 0;
  for (uint32_t i=2; i<B; i++) {
    if (!GET(i)) {
      N++;
      for (uint64_t j=2*(uint64_t)i; j<B; j+=i) {
	SET(j);
      }
    }
  }

  // Add prime powers
  for (uint32_t i=ceil(sqrt(B)); i>=2; i--) {
    if (!GET(i)) {
      for (uint64_t j=(uint64_t)i*i; j<B; j*=i) {
	assert(GET(j));
	UNSET(j);
	N++;
      }
    }
  }
  
  prime_data *primes = malloc(sizeof(prime_data) + N*sizeof(primes->data[0]));
  assert(primes);
  primes->B = B;
  primes->N = N;
  N = 0;

  for (uint32_t i=2; i<B; i++) {
    if (!GET(i)) {
      for (uint64_t j=i; j<B; j*=i) {
	primes->data[N].factor = j;
	primes->data[N].prime = i;
	primes->data[N].phi = j - j/i;
	N++;
	assert(!GET(j));
	SET(j);
      }
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (uint32_t i=0; i<primes->N; i++) {
    uint32_t p = primes->data[i].prime;
    uint32_t f = primes->data[i].factor;
    primes->data[i].lg = FIXEDP(log2f(p));
    primes->data[i].fast_d = libdivide_u64_gen(f);
    uint32_t t = 1;
    for (int j=0; j<8; j++) {
      primes->data[i].reduce32[j] = t;
      t = ((uint64_t)t<<32) % f;
    }
  }

  free(sieve);
  return primes;
}

int main (int argc, char *argv[]) {
#if PSIZE==256
  // 256-bit parameters
  uint32_t B = 1ULL<<27;
  uint64_t C = 1ULL<<21;

  // 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
  word P = ((((word)0x4000000000000000LL << 128) + 0x224698fc094cf91bLL) << 64) +0x992d30ed00000001LL;
#else
  // 128-bit parameters
  uint32_t B = 1ULL<<14;
  uint64_t C = 1ULL<<9;

  // 0x30000003000000010000000000000001
  word P = (((word)0x3000000300000001)<<64) + 1;
#endif

  printf ("Parameters:\nB = %#llx C = %#llx\n",
	  (unsigned long long) B, (unsigned long long) C);
  printf ("P = ");
  print_word(P);
  printf ("\n");

  char *file = "/tmp/root.data";
  if (argc == 2) {
    file = argv[1];
  }
  printf ("Using file '%s' for precomputation\n", file);  

  prime_data *primes;

  int fd = open(file, O_RDONLY);
  if (fd == -1) {
    perror("Error reading file");
    goto BAD_FILE;    
  }
  struct stat st;
  int r = fstat(fd, &st);
  if (r) {
    perror("Error stat-ing file");
    close(fd);
    goto BAD_FILE;    
  }
  primes = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (primes == MAP_FAILED) {
    perror("mmap failed");
    close(fd);
    goto BAD_FILE;    
  }

  if (primes->B != B || st.st_size != sizeof(*primes) + primes->N*sizeof(primes->data[0])) {
    printf ("Bad file.\n");
  BAD_FILE:
    printf ("Starting precomputation.\n");
    primes = precompute(B);
    FILE *F = fopen(file, "w");
    if (!F) {
      
    } else {
      fwrite(primes, sizeof(*primes) + primes->N*sizeof(primes->data[0]), 1, F);
      printf ("Precomputed data written to '%s'.\n", file);
    }
  } else {
    printf ("Reusing precomputation.\n");    
  }

  
  printf ("#primes: 2^%f\n", log2(primes->N));

  double work = 0;
  for (uint32_t i=0; i<primes->N; i++) {
    uint32_t p = primes->data[i].factor;
    work += ceil(1.0*C/p);
  }
  printf ("Expected cost: 2^%f\n", log2(work));
  
  word x;
  getrandom(&x, sizeof(x), 0);
  x %= P;
  
  printf ("x = ");
  print_word(x);
  printf ("\n");


  // Rational reconstruction using extended Euclid
  signed_word u1 = 0, v1 = P;
  signed_word u2 = 1, v2 = x;

  while (abs_word(v2) > abs_word(u2)) {
    signed_word div = v1/v2;
    signed_word ut = u2, vt = v2;
    u2 = u1 - div*u2;
    v2 = v1 - div*v2;
    u1 = ut; v1 = vt;
  }

  // Consider series of fractions
  // x = (v1 + i.v2) / (u1 + i.u2)

  uint64_t steps = 0;

  /**** Find i with smooth Vi ****/
  // Note: v1, v2 > 0
  assert(v1 > 0);
  assert(v2 > 0);
  
  int32_t *Vsmooth = calloc(C, sizeof(Vsmooth[0]));
  assert(Vsmooth);
  Vsmooth[0] = FIXEDP(log_word(v1));
  uint32_t lgV = FIXEDP(log_word(v2));
  for (uint64_t i=1; i<C; i++) {
    Vsmooth[i] = lgV+FIXEDP(log2f(i));
  }

  /* Main loop */
#pragma omp parallel for schedule(dynamic) reduction (+:steps)
  for (uint32_t i=0; i<primes->N; i++) {
    uint32_t p = primes->data[i].factor;
    uint32_t *r32 = primes->data[i].reduce32;
    uint32_t r = MOD(v1, p, r32);
    uint32_t v2p = MOD(v2, p, r32);
    uint32_t mod_inv = p - powermod(v2p, primes->data[i].phi-1, p, &primes->data[i].fast_d);
    if (((uint64_t)v2p*mod_inv)%p == p-1) {
      for (uint64_t j=((uint64_t)r*mod_inv)%p; j<C; j+=p) {
	__atomic_add_fetch(&Vsmooth[j], -primes->data[i].lg, __ATOMIC_RELAXED);
	// assert((v1+j*v2) % p == 0);
	steps++;
      }
    }
  }

  /**** Find i with smooth Ui ****/
  
  int32_t *Usmooth = calloc(C, sizeof(Usmooth[0]));
  assert(Usmooth);
  Usmooth[0] = FIXEDP(log_word(abs_word(u1)));
  uint32_t lgU = FIXEDP(log_word(abs_word(u2)));
  for (uint64_t i=1; i<C; i++) {
    Usmooth[i] = lgU+FIXEDP(log2f(i));
  }

  /* Main loop */
#pragma omp parallel for schedule(dynamic) reduction (+:steps)
  for (uint32_t i=0; i<primes->N; i++) {
    uint32_t p = primes->data[i].factor;
    uint32_t *r32 = primes->data[i].reduce32;
    uint32_t r = SMOD(u1, p, r32);
    uint32_t u2p = SMOD(u2, p, r32);
    uint32_t mod_inv = p - powermod(u2p, primes->data[i].phi-1, p, &primes->data[i].fast_d);
    if (((uint64_t)u2p*mod_inv)%p == p-1) {
      for (uint64_t j=((uint64_t)r*mod_inv)%p; j<C; j+=p) {
	__atomic_add_fetch(&Usmooth[j], -primes->data[i].lg, __ATOMIC_RELAXED);
	// assert((u1+j*u2) % p == 0);
	steps++;
      }
    }
  }

  printf ("Steps: 2^%f\n", log2(steps));
  int win = 0;
  
  int32_t min = Usmooth[0]+Vsmooth[0], minat = 0;
  for (uint64_t i=0; i<C; i++) {
    int32_t smooth = Usmooth[i] > Vsmooth[i]? Usmooth[i]: Vsmooth[i];
    if (smooth < FIXEDP(45)) {
      printf ("[%10lu]: %f %f\n", i,
	      Vsmooth[i]/((float) FIXEDP(1)), Usmooth[i]/((float) FIXEDP(1)));
      win=1;
    }
    if (smooth < min) {
      min = smooth;
      minat = i;
    }
  }
  printf("Best:\n");
  printf ("[%10u]: %f %f\n", minat,
	  Vsmooth[minat]/((float) FIXEDP(1)), Usmooth[minat]/((float) FIXEDP(1)));

  word v = v1 + minat*v2;
  signed_word u = u1 + minat*u2;
  printf ("x = %s", u>0? "": "-");
  print_word(v);
  printf ("/");
  print_word(abs_word(u));
  printf ("\n");  

  if (u>0)
    assert(((_BitInt(512))x*u)%P == v);
  else
    assert(((_BitInt(512))x*-u)%P == P-v);
  
  if (win) {
    printf ("Successfull smooth decomposition\n");
    return 0;
  }
  else {
    printf ("Failed to find smooth decomposition\n");
    return 1;
  }
}
