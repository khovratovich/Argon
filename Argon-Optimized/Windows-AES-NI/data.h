
//#define KAT
//#define KAT_INTERNAL
#define _MEASURE
//#define _MEASUREINT 
#define THREADING
#define RANDOMIZE //make slice selection pseudo-random


#define MAX_THREADS 32
#define MAX_SLICE 32
#define MIN_SLICE 1
#define MAX_OUTLEN 0xFFFFFFFF
#define MIN_OUTLEN 4
#define MIN_MEMORY 1
#define MAX_MEMORY 0xFFFFFF
#define MIN_TIME 2
#define LENGTH_SIZE 4
#define MIN_PASSWORD 0
#define MAX_PASSWORD 0xFFFFFFFF
#define MAX_AD 0xFFFFFFFF
#define MAX_SALT  0xFFFFFFFF
#define MIN_SALT  8
#define MIN_SECRET  0
#define MAX_SECRET 16
#define MIN_AD  0
#define INPUT_SIZE (INPUT_BLOCKS*12)
#define INPUT_BLOCKS 32
#define CACHE_SIZE 128
#define MAX_CACHE 128
#define BATCH_SIZE 16   //should be not larger than 32
#define GROUP_SIZE 32


#define AES_ROUNDS 5


struct int128{
	uint64_t i0, i1;
	int128(uint64_t y0 = 0, uint64_t y1 = 0){ i0 = y0; i1 = y1; };
	int128& operator^=(const int128 &r){ i0 ^= r.i0; i1 ^= r.i1; return *this; }
	int128& operator=(const int128 &r){ i0 = r.i0; i1 = r.i1; return *this; }
	unsigned char operator[](unsigned i)
	{
		if (i<8)
			return (i0 >> (8 * i)) & 0xff;
		else if (i<16)
			return (i1 >> (8 * (i - 8))) & 0xff;
		return 0;
	}
	int128 operator^(const int128 &r){ return int128(i0 ^ r.i0, i1^r.i1); }
};

extern  void AES_reduced_batch_intr(__m128i* batch, uint32_t batch_size);
extern int ArgonFast64Ext(uint8_t *out, uint32_t outlen, const uint8_t *pwd, uint32_t pwdlen, const uint8_t *salt, uint32_t saltlen, const uint8_t *secret,
uint8_t secretlen, const void *ad, uint32_t adlen, uint32_t t_cost, uint32_t m_cost, uint8_t parallel_degree);
extern int PHS(void *out, size_t outlen, const void *in, size_t inlen, const void *salt, size_t saltlen,
	uint32_t t_cost, size_t m_cost, uint32_t thread_n);