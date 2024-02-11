//toy.h:
//toy Post-Quantum Public-Key Cryptosystem
#define TK_K 3
#define TK_N 4
#define TK_Q 97
#define invq 14
#define NEG(X) (TK_Q - (X))
#define w 22
#define sqrt_w 33
#define beta_q ((1 << 32) % q)
#define MONT_MUL(x, y) \
  (x = ((uint64_t)x * y) >> 32, x = (x >> 16) - (((x & 0xFFFF) * invq & 0xFFFF) * q >> 16))

static void toy_fill_small(short *buf, int n);
static void toy_polmul_naive(short* dst, const short* a, const short* b, int add);
static void toy_mulmv(short *dst, const short *mat, const short* vec);
static void toy_dot(short* dst, const short* v1, const short* v2);
static void toy_add(short* dst, const short* v1, const short* v2, int count, int v2_neg);
void toy_gen(short* A, short* t, short* s);
void toy_enc(const short *A, const short* t, int plain, short *u, short* v);
int toy_dec(const short* s, const short* u, const short* v);
//*******************************************
void reference_ntt(int *data, int forrward);
void permute_bitreverse( const short *data,  short *x);
void copy_array(const short *x,  short *data);
void ntt( short *data,int forwarrd,int anti_cyclic);

