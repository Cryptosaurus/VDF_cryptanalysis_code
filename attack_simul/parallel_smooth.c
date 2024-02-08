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

typedef unsigned _BitInt(256+32) word;

#define PSIZE 256
// #define PSIZE 128

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

typedef struct {
  uint32_t B;
  uint32_t N;
  struct {
    uint32_t factor;
    uint32_t prime;
    uint32_t mod_inv;
    uint32_t lg;
  } data[0];
} prime_data;

// Fixed point arithmetic
#define FIXEDP(x) ((int32_t)((x)*(1<<16)))

uint32_t modinv(uint32_t x, uint32_t q) {
  // Extended Euclid
  int64_t u1 = 0, v1 = q;
  int64_t u2 = 1, v2 = x;

  while (v2 > 1) {
    uint64_t div = v1/v2;
    int64_t ut = u2, vt = v2;
    u2 = u1 - div*u2;
    v2 = v1 - div*v2;
    u1 = ut; v1 = vt;
  }

  assert(u2 <= q);
  assert(-u2 < q);
  return u2>0? u2: u2+q;
}

prime_data *precompute(uint32_t B, word P) {
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
	N++;
	assert(!GET(j));
	SET(j);
      }
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (uint32_t i=0; i<primes->N; i++) {
    uint32_t f = primes->data[i].factor;
    primes->data[i].mod_inv = f - modinv(P%f, f);
    uint32_t p = primes->data[i].prime;
    primes->data[i].lg = FIXEDP(log2f(p));
  }

  free(sieve);
  return primes;
}

int main (int argc, char *argv[]) {
#if PSIZE==256
  // 256-bit parameters
  uint32_t B = (uint32_t)-1; // Smooth part          ($B$ in the paper)
  uint64_t C = 1ULL<<26;     // Number of candidates ($R$ in the paper)
  int32_t BB = FIXEDP(45);   // Medium-factor size   ($B'$ in the paper)
  
  // 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
  word P = ((((word)0x4000000000000000LL << 128) + 0x224698fc094cf91bLL) << 64) +0x992d30ed00000001LL;
  int32_t lg = FIXEDP(254);
#else
  uint32_t B = 1ULL<<20;
  uint64_t C = 1ULL<<16;
  int32_t BB = FIXEDP(40);

  // 0x30000003000000010000000000000001
  word P = (((word)0x3000000300000001)<<64) + 1;
  int32_t lg = FIXEDP(125.584962506096);
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
    primes = precompute(B, P);
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

  int32_t *smoothness = calloc(C, sizeof(smoothness[0]));
  assert(smoothness);
  smoothness[0] = FIXEDP(log_word(x));
  for (uint64_t i=1; i<C; i++) {
    smoothness[i] = lg+FIXEDP(log2f(i));
  }


  /* Main loop */
  uint64_t steps = 0;
  
#pragma omp parallel for schedule(dynamic) reduction (+:steps)
  for (uint32_t i=0; i<primes->N; i++) {
    uint32_t p = primes->data[i].factor;
    uint32_t r = x % p;
    for (uint64_t j=((uint64_t)r*primes->data[i].mod_inv)%p; j<C; j+=p) {
      __atomic_add_fetch(&smoothness[j], -primes->data[i].lg, __ATOMIC_RELAXED);
      steps++;
    }
  }

  printf ("Steps: 2^%f\n", log2(steps));
  int win = 0;
  
  int32_t min = smoothness[0], minat = 0;
  for (uint64_t i=0; i<C; i++) {
    if (smoothness[i] < BB) {
      printf ("[%10lu]: %f\n", i, smoothness[i]/((float) FIXEDP(1)));
      win=1;
    }
    if (smoothness[i] < min) {
      min = smoothness[i];
      minat = i;
    }
  }
  printf("Best:\n");
  printf ("[%10u]: %f\n", minat, smoothness[minat]/((float) FIXEDP(1)));

  printf ("y = ");
  print_word(x);
  printf (" + %i*", minat);
  print_word(P);
  printf (" is almost-smooth\n");  
  
  if (win) {
    printf ("Successfull smooth decomposition\n");
    return 0;
  }
  else {
    printf ("Failed to find smooth decomposition\n");
    return 1;
  }
}
