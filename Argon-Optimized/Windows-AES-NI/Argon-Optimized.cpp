// Argon-Optimized.cpp
// Optimized code for Windows x86/x64 architectures. Produced with Visual Studio 2013 Express
// Code written by Dmitry Khovratovich, khovratovich@gmail.com, and Yann Le Corre
// Requires C++11
//To fileprint intermediate variables, uncomment #define KAT
//To enable measurement of SubGroups and ShuffleSlices, uncomment #define _MEASUREINT
//Threading is always enabled


#include "stdio.h"
#include "stdint.h"
#include <time.h>
#include <vector>
#include <thread>
#include <random>
#include <cstring>
using namespace std;

#include "wmmintrin.h"
#if defined(_MSC_VER)
#include "intrin.h"
#else
#include <x86intrin.h>
#endif
#include "emmintrin.h"

#pragma intrinsic(_mm_set_epi64x)  

#include "data.h"
#include "blake2.h"





void dumpState(FILE *fp, __m128i *state, size_t state_size)
{
	for(size_t i = 0; i < state_size; ++i)
	{
		uint64_t u0 = (uint64_t)_mm_extract_epi64(state[i], 0);
		uint64_t u1 = (uint64_t)_mm_extract_epi64(state[i], 1);
#if defined(_MSC_VER)
		fprintf(fp, "Block %3.3lu: H: %.16llx L: %.16llx \n", i, u1, u0);
#else
		fprintf(fp, "Block %3.3lu: H: %" PRIx64 " L: %" PRIx64 "\n", i, u1, u0);
#endif
		//fprintf(fp,"Block %3.3lu: H: %.16lx L: %.16lx\n", i, u1, u0);
	}
}



void ShuffleSlicesThr(__m128i* state, uint32_t slice_length, unsigned slices)
{
	for(uint32_t s=0; s<slices; ++s) //Loop on slices
	{
		uint32_t j = slice_length-1;
		for(uint32_t i=0; i<slice_length; ++i)
		{
			uint32_t index1 = s*(slice_length)+i;
			__m128i v1 = state[index1];
			uint64_t q1 = _mm_extract_epi64(v1, 0);
			uint32_t index2 = s*(slice_length)+j;
			__m128i v2 = state[index2];
			uint64_t q2 = _mm_extract_epi64(v2, 0);
			uint64_t value = q1 + q2;
			j = value% slice_length;
			state[index1]=v2;
			state[index2] = v1;
		}
	}
}







void SubGroupsThrExtended(__m128i* state, uint32_t group_n, uint32_t distance, uint64_t* Input, uint32_t shift) //SubGroupsThr prepended by Initial round
{

	__m128i groups[MAX_CACHE][GROUP_SIZE];
	uint32_t cached = 0;
	cached = (group_n<CACHE_SIZE) ? group_n : CACHE_SIZE;
	for (uint32_t i = 0; i<group_n; i += cached)
	{
		//Storing group inputs
		if (group_n< cached + i)
			cached = (group_n - i);
		for (unsigned l = 0; l<GROUP_SIZE; ++l)
		{
			for (uint32_t j = 0; j<cached; ++j)
			{
				groups[j][l] = _mm_set_epi64x((uint64_t)(i + j + shift + l*distance), Input[(i + j + shift + l*distance) % 4]); //I[l]||((i+j)*GROUP_SIZE+l)
			}
		}
		AES_reduced_batch_intr(groups[0], GROUP_SIZE*cached);
		for (unsigned j = 0; j<cached; ++j)
		{
			//Computing X_i:
			__m128i X[16];
			X[0] = _mm_xor_si128(groups[j][3], _mm_xor_si128(groups[j][7], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][15], _mm_xor_si128(groups[j][19], _mm_xor_si128(groups[j][23], _mm_xor_si128(groups[j][27], groups[j][31])))))));
			X[1] = _mm_xor_si128(groups[j][1], _mm_xor_si128(groups[j][3], _mm_xor_si128(groups[j][9], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][17], _mm_xor_si128(groups[j][19], _mm_xor_si128(groups[j][25], groups[j][27])))))));
			X[2] = _mm_xor_si128(groups[j][0], _mm_xor_si128(groups[j][2], _mm_xor_si128(groups[j][4], _mm_xor_si128(groups[j][6], _mm_xor_si128(groups[j][16], _mm_xor_si128(groups[j][18], _mm_xor_si128(groups[j][20], groups[j][22])))))));
			X[3] = _mm_xor_si128(groups[j][1], _mm_xor_si128(groups[j][3], _mm_xor_si128(groups[j][5], _mm_xor_si128(groups[j][7], _mm_xor_si128(groups[j][9], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][13], groups[j][15])))))));
			X[4] = _mm_xor_si128(groups[j][6], _mm_xor_si128(groups[j][7], _mm_xor_si128(groups[j][14], _mm_xor_si128(groups[j][15], _mm_xor_si128(groups[j][22], _mm_xor_si128(groups[j][23], _mm_xor_si128(groups[j][30], groups[j][31])))))));
			X[5] = _mm_xor_si128(groups[j][10], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][14], _mm_xor_si128(groups[j][15], _mm_xor_si128(groups[j][26], _mm_xor_si128(groups[j][27], _mm_xor_si128(groups[j][30], groups[j][31])))))));
			X[6] = _mm_xor_si128(groups[j][16], _mm_xor_si128(groups[j][17], _mm_xor_si128(groups[j][20], _mm_xor_si128(groups[j][21], _mm_xor_si128(groups[j][24], _mm_xor_si128(groups[j][25], _mm_xor_si128(groups[j][28], groups[j][29])))))));
			X[7] = _mm_xor_si128(groups[j][12], _mm_xor_si128(groups[j][13], _mm_xor_si128(groups[j][14], _mm_xor_si128(groups[j][15], _mm_xor_si128(groups[j][28], _mm_xor_si128(groups[j][29], _mm_xor_si128(groups[j][30], groups[j][31])))))));
			X[8] = _mm_xor_si128(groups[j][4], _mm_xor_si128(groups[j][5], _mm_xor_si128(groups[j][6], _mm_xor_si128(groups[j][7], _mm_xor_si128(groups[j][12], _mm_xor_si128(groups[j][13], _mm_xor_si128(groups[j][14], groups[j][15])))))));
			X[9] = _mm_xor_si128(groups[j][16], _mm_xor_si128(groups[j][17], _mm_xor_si128(groups[j][18], _mm_xor_si128(groups[j][19], _mm_xor_si128(groups[j][20], _mm_xor_si128(groups[j][21], _mm_xor_si128(groups[j][22], groups[j][23])))))));
			X[10] = _mm_xor_si128(groups[j][1], _mm_xor_si128(groups[j][5], _mm_xor_si128(groups[j][9], _mm_xor_si128(groups[j][13], _mm_xor_si128(groups[j][17], _mm_xor_si128(groups[j][21], _mm_xor_si128(groups[j][25], groups[j][29])))))));
			X[11] = _mm_xor_si128(groups[j][2], _mm_xor_si128(groups[j][6], _mm_xor_si128(groups[j][10], _mm_xor_si128(groups[j][14], _mm_xor_si128(groups[j][18], _mm_xor_si128(groups[j][22], _mm_xor_si128(groups[j][26], groups[j][30])))))));
			X[12] = _mm_xor_si128(groups[j][4], _mm_xor_si128(groups[j][5], _mm_xor_si128(groups[j][6], _mm_xor_si128(groups[j][7], _mm_xor_si128(groups[j][20], _mm_xor_si128(groups[j][21], _mm_xor_si128(groups[j][22], groups[j][23])))))));
			X[13] = _mm_xor_si128(groups[j][8], _mm_xor_si128(groups[j][9], _mm_xor_si128(groups[j][10], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][24], _mm_xor_si128(groups[j][25], _mm_xor_si128(groups[j][26], groups[j][27])))))));
			X[14] = _mm_xor_si128(groups[j][0], _mm_xor_si128(groups[j][1], _mm_xor_si128(groups[j][2], _mm_xor_si128(groups[j][3], _mm_xor_si128(groups[j][8], _mm_xor_si128(groups[j][9], _mm_xor_si128(groups[j][10], groups[j][11])))))));
			X[15] = _mm_xor_si128(groups[j][0], _mm_xor_si128(groups[j][4], _mm_xor_si128(groups[j][8], _mm_xor_si128(groups[j][12], _mm_xor_si128(groups[j][16], _mm_xor_si128(groups[j][20], _mm_xor_si128(groups[j][24], groups[j][28])))))));





			for (unsigned k = 0; k<16; k += BATCH_SIZE)
			{
				AES_reduced_batch_intr(X + k, BATCH_SIZE);//Computing F's
			}
			for (unsigned k = 0; k<16; ++k)
			{
				groups[j][2 * k] = _mm_xor_si128(groups[j][2 * k], X[k]); //XORs
				groups[j][2 * k + 1] = _mm_xor_si128(groups[j][2 * k + 1], X[k]);
			}
			for (unsigned k = 0; k<32; k += BATCH_SIZE)
			{
				AES_reduced_batch_intr(groups[j] + k, BATCH_SIZE);//Computing F's
			}
		}//end of group computation
		for (unsigned l = 0; l<32; ++l)
		{
			for (unsigned j = 0; j<cached; ++j)
			{
				state[i + j + l*distance] = groups[j][l];
			}
		}
	}
}


void SubGroupsThr(__m128i* state, size_t group_n, size_t distance) //elements of the group are 'distance' blocks from each other in memory
{

	__m128i groups[MAX_CACHE][GROUP_SIZE];
	size_t cached=0;
	cached = (group_n<CACHE_SIZE) ? group_n : CACHE_SIZE;
	for (size_t i = 0; i<group_n; i += cached)
	{
		//Storing group inputs
		if(group_n< cached+i)
			cached = (group_n-i);
		for(unsigned l=0; l<GROUP_SIZE; ++l)
		{
			for (size_t j = 0; j<cached; ++j)
			{
				groups[j][l]= state[i+j+l*distance];
			}
		}
		for(unsigned j=0; j<cached; ++j)
		{
			//Computing X_i:
			__m128i X[16];
			X[ 0] = _mm_xor_si128(groups[j][ 3], _mm_xor_si128(groups[j][ 7], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][15], _mm_xor_si128(groups[j][19], _mm_xor_si128(groups[j][23], _mm_xor_si128(groups[j][27], groups[j][31])))))));
			X[ 1] = _mm_xor_si128(groups[j][ 1], _mm_xor_si128(groups[j][ 3], _mm_xor_si128(groups[j][ 9], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][17], _mm_xor_si128(groups[j][19], _mm_xor_si128(groups[j][25], groups[j][27])))))));
			X[ 2] = _mm_xor_si128(groups[j][ 0], _mm_xor_si128(groups[j][ 2], _mm_xor_si128(groups[j][ 4], _mm_xor_si128(groups[j][ 6], _mm_xor_si128(groups[j][16], _mm_xor_si128(groups[j][18], _mm_xor_si128(groups[j][20], groups[j][22])))))));
			X[ 3] = _mm_xor_si128(groups[j][ 1], _mm_xor_si128(groups[j][ 3], _mm_xor_si128(groups[j][ 5], _mm_xor_si128(groups[j][ 7], _mm_xor_si128(groups[j][ 9], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][13], groups[j][15])))))));
			X[ 4] = _mm_xor_si128(groups[j][ 6], _mm_xor_si128(groups[j][ 7], _mm_xor_si128(groups[j][14], _mm_xor_si128(groups[j][15], _mm_xor_si128(groups[j][22], _mm_xor_si128(groups[j][23], _mm_xor_si128(groups[j][30], groups[j][31])))))));
			X[ 5] = _mm_xor_si128(groups[j][10], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][14], _mm_xor_si128(groups[j][15], _mm_xor_si128(groups[j][26], _mm_xor_si128(groups[j][27], _mm_xor_si128(groups[j][30], groups[j][31])))))));
			X[ 6] = _mm_xor_si128(groups[j][16], _mm_xor_si128(groups[j][17], _mm_xor_si128(groups[j][20], _mm_xor_si128(groups[j][21], _mm_xor_si128(groups[j][24], _mm_xor_si128(groups[j][25], _mm_xor_si128(groups[j][28], groups[j][29])))))));
			X[ 7] = _mm_xor_si128(groups[j][12], _mm_xor_si128(groups[j][13], _mm_xor_si128(groups[j][14], _mm_xor_si128(groups[j][15], _mm_xor_si128(groups[j][28], _mm_xor_si128(groups[j][29], _mm_xor_si128(groups[j][30], groups[j][31])))))));
			X[ 8] = _mm_xor_si128(groups[j][ 4], _mm_xor_si128(groups[j][ 5], _mm_xor_si128(groups[j][ 6], _mm_xor_si128(groups[j][ 7], _mm_xor_si128(groups[j][12], _mm_xor_si128(groups[j][13], _mm_xor_si128(groups[j][14], groups[j][15])))))));
			X[ 9] = _mm_xor_si128(groups[j][16], _mm_xor_si128(groups[j][17], _mm_xor_si128(groups[j][18], _mm_xor_si128(groups[j][19], _mm_xor_si128(groups[j][20], _mm_xor_si128(groups[j][21], _mm_xor_si128(groups[j][22], groups[j][23])))))));
			X[10] = _mm_xor_si128(groups[j][ 1], _mm_xor_si128(groups[j][ 5], _mm_xor_si128(groups[j][ 9], _mm_xor_si128(groups[j][13], _mm_xor_si128(groups[j][17], _mm_xor_si128(groups[j][21], _mm_xor_si128(groups[j][25], groups[j][29])))))));
			X[11] = _mm_xor_si128(groups[j][ 2], _mm_xor_si128(groups[j][ 6], _mm_xor_si128(groups[j][10], _mm_xor_si128(groups[j][14], _mm_xor_si128(groups[j][18], _mm_xor_si128(groups[j][22], _mm_xor_si128(groups[j][26], groups[j][30])))))));
			X[12] = _mm_xor_si128(groups[j][ 4], _mm_xor_si128(groups[j][ 5], _mm_xor_si128(groups[j][ 6], _mm_xor_si128(groups[j][ 7], _mm_xor_si128(groups[j][20], _mm_xor_si128(groups[j][21], _mm_xor_si128(groups[j][22], groups[j][23])))))));
			X[13] = _mm_xor_si128(groups[j][ 8], _mm_xor_si128(groups[j][ 9], _mm_xor_si128(groups[j][10], _mm_xor_si128(groups[j][11], _mm_xor_si128(groups[j][24], _mm_xor_si128(groups[j][25], _mm_xor_si128(groups[j][26], groups[j][27])))))));
			X[14] = _mm_xor_si128(groups[j][ 0], _mm_xor_si128(groups[j][ 1], _mm_xor_si128(groups[j][ 2], _mm_xor_si128(groups[j][ 3], _mm_xor_si128(groups[j][ 8], _mm_xor_si128(groups[j][ 9], _mm_xor_si128(groups[j][10], groups[j][11])))))));
			X[15] = _mm_xor_si128(groups[j][ 0], _mm_xor_si128(groups[j][ 4], _mm_xor_si128(groups[j][ 8], _mm_xor_si128(groups[j][12], _mm_xor_si128(groups[j][16], _mm_xor_si128(groups[j][20], _mm_xor_si128(groups[j][24], groups[j][28])))))));



			

			for(unsigned k=0; k<16;k+=BATCH_SIZE)
			{
				AES_reduced_batch_intr(X + k, BATCH_SIZE);//Computing F's
			}
			for(unsigned k=0; k<16;++k)
			{
				groups[j][2*k] =  _mm_xor_si128(groups[j][2*k],X[k]); //XORs
				groups[j][2*k+1] = _mm_xor_si128(groups[j][2*k+1],X[k]); 
			}
			for(unsigned k=0; k<32;k+=BATCH_SIZE)
			{
				AES_reduced_batch_intr(groups[j] + k, BATCH_SIZE);//Computing F's
			}
		}//end of group computation
		for(unsigned l=0; l<32; ++l)
		{
			for(unsigned j=0; j<cached; ++j)
			{
				state[i+j+l*distance]=groups[j][l];
			}
		}
	}
}


void GenPermutation32(uint8_t out[32])
{
	if (out == NULL)
		return;
	//Assume that memory for out is allocated
	vector<bool> filled(32, false);

	std::mt19937 gen;  //Starting PRNG
	gen.seed(time(0));
	for (unsigned i = 0; i < 32; ++i)
	{
		uint8_t value = gen();
		uint8_t index = 0;
		while (value != 0)//Go over unassigned values for the permutation
		{
			index++;
			if (index == 32)
				index = 0;
			if (!filled[index])//Skip allocated values
				value--;		
		}
		filled[index] = true;
		out[i] = index;
	}

}


int ArgonFast64Ext(uint8_t *out, uint32_t outlen, const uint8_t *pwd, uint32_t pwdlen, const uint8_t *salt, uint32_t saltlen, const uint8_t *secret,
	uint8_t secretlen, const void *ad, uint32_t adlen, uint32_t t_cost, uint32_t m_cost, uint8_t parallel_degree)
{
	__m128i* state;
#ifdef KAT
	FILE* fp = fopen("kat-opt.log", "a+");

	fprintf(fp, "=======================================\n");
	fprintf(fp, "Iterations: %d, Memory: %d KBytes, Parallelism: %d slices, Tag length: %d bytes\n", t_cost, m_cost, parallel_degree, outlen);

#endif
#ifdef _MEASUREINT
	uint64_t  i1, i2, i3, i4, i5, i6, d1, d2;
	uint32_t ui1, ui2, ui3, ui4, ui6;
	i1 = __rdtscp(&ui1);
#endif
	//0. Restricting parameters
	if (outlen>MAX_OUTLEN)
		outlen = MAX_OUTLEN;
	if (outlen < MIN_OUTLEN)
		return -1;  //Tag too short

	if (pwdlen> MAX_PASSWORD)
		pwdlen = MAX_PASSWORD;
	if (pwdlen < MIN_PASSWORD)
		return -2; //Password too short

	if (saltlen < MIN_SALT)
		return -3; //Salt too short
	if (saltlen> MAX_SALT)
		saltlen = MAX_SALT;

	if (secretlen> MAX_SECRET)
		secretlen = MAX_SECRET;
	if (secretlen < MIN_SECRET)
		return -4; //Secret too short

	if (adlen> MAX_AD)
		adlen = MAX_AD;
	if (adlen < MIN_AD)
		return -5; //Associated data too short

	//minumum m_cost =1
	if (m_cost<MIN_MEMORY)
		m_cost = MIN_MEMORY;
	if (m_cost>MAX_MEMORY)
		m_cost = MAX_MEMORY;

	//minimum t_cost =3
	if (t_cost<MIN_TIME)
		t_cost = MIN_TIME;

	if (parallel_degree<MIN_SLICE)
		parallel_degree = MIN_SLICE;
	if (parallel_degree>MAX_SLICE)
		parallel_degree = MAX_SLICE;


#ifdef KAT
	fprintf(fp, "Password: ");
	for (unsigned i = 0; i<pwdlen; ++i)
		fprintf(fp, "%2.2x ", ((unsigned char*)pwd)[i]);
	fprintf(fp, "\n");
	fprintf(fp, "Salt: ");
	for (unsigned i = 0; i<saltlen; ++i)
		fprintf(fp, "%2.2x ", ((unsigned char*)salt)[i]);
	fprintf(fp, "\n");

#endif

	//Initial hashing

	uint64_t input_hash[4];

	blake2b_state BlakeHash;
	blake2b_init(&BlakeHash, 32);

	blake2b_update(&BlakeHash, (const uint8_t*)&pwdlen, sizeof(pwdlen));
	blake2b_update(&BlakeHash, (const uint8_t*)pwd, pwdlen);
	blake2b_update(&BlakeHash, (const uint8_t*)&saltlen, sizeof(saltlen));
	blake2b_update(&BlakeHash, (const uint8_t*)salt, saltlen);
	blake2b_update(&BlakeHash, (const uint8_t*)&secretlen, sizeof(secretlen));
	blake2b_update(&BlakeHash, (const uint8_t*)secret, secretlen);
	blake2b_update(&BlakeHash, (const uint8_t*)&adlen, sizeof(adlen));
	blake2b_update(&BlakeHash, (const uint8_t*)ad, adlen);
	blake2b_update(&BlakeHash, (const uint8_t*)&parallel_degree, sizeof(parallel_degree));
	blake2b_update(&BlakeHash, (const uint8_t*)&m_cost, sizeof(m_cost));
	blake2b_update(&BlakeHash, (const uint8_t*)&t_cost, sizeof(t_cost));
	blake2b_update(&BlakeHash, (const uint8_t*)&outlen, sizeof(outlen));

	blake2b_final(&BlakeHash, (uint8_t*)input_hash, 32);


#if defined(KAT) && defined(KAT_INTERNAL)
	fprintf(fp, "Input hash:\n");
	for (unsigned i = 0; i<4; ++i)
	{
		fprintf(fp, "%.16llx ", input_hash[i]);
	}
	fprintf(fp, "\n");

#endif


	//2. Filling blocks
	size_t state_size = m_cost * 64;
#ifdef _MEASUREINT
	i3 = __rdtscp(&ui3);
	d1 = (i3 - i1) / (state_size);
	printf("First phase %2.2f cpb\n", (float)d1 / (16));
#endif
	state = new __m128i[state_size];
	if (state == NULL)
		return 1;
	//printf("Memory allocated: %lu MBytes, %d threads\n", (state_size*sizeof(__m128i)) / (1 << 20), parallel_degree);

	size_t width = state_size / GROUP_SIZE;

	

#ifdef _MEASUREINT
	i2 = __rdtscp(&ui2);
	d1 = (i2 - i3) / (state_size);
	printf("Filling phase %2.2f cpb\n", (float)d1 / (16));
#endif

	


	//4.1 First call to SubGroups: 
	vector<thread> Threads;

	//Each thread processes a certain number of groups
	size_t distance = state_size / GROUP_SIZE;
	size_t groups_per_thread = distance / (parallel_degree);//Area size of each thread.
	for (unsigned i = 0; i<parallel_degree - 1; ++i)
		Threads.push_back(thread(SubGroupsThrExtended, state + i*groups_per_thread, groups_per_thread, distance, input_hash, i*groups_per_thread));
	//Last thread with possibly larger area
	Threads.push_back(thread(SubGroupsThrExtended, state + (parallel_degree - 1)*groups_per_thread, distance - (parallel_degree - 1)*groups_per_thread, 
		distance, input_hash, (parallel_degree - 1)*groups_per_thread));
	for (auto& th : Threads)
		th.join();
	Threads.clear();
#ifdef _MEASUREINT
	i5 = __rdtscp(&ui3);
	d1 = (i5 - i4) / (state_size);
	printf("SubGroups %2.2f cpb\n", (float)d1 / (16));
#endif
#if defined(KAT) && defined(KAT_INTERNAL)
	fprintf(fp, "Round %d SubGroups:\nBlocks:\n", 1);
	dumpState(fp, state, state_size);
#endif



//4.2 Other rounds
memset(input_hash, 0, 4*sizeof(uint64_t));
for (unsigned l = 0; l <t_cost; ++l)
{
#ifdef RANDOMIZE
	uint8_t perm[32];
	GenPermutation32(perm);
	size_t slice_length = state_size / GROUP_SIZE;
	for(unsigned i=0; i<32; ++i)
	{
		Threads.push_back(thread(ShuffleSlicesThr, state + i*slice_length, slice_length, 1));
		if (((i + 1)%thread_n == 0) || (i+1 == 32))
		{
			for (auto& th : Threads)
				th.join();
			Threads.clear();
		}
	}
#else
	uint32_t slice_length = (state_size / parallel_degree);
	for (unsigned i = 0; i<parallel_degree - 1; ++i)
		Threads.push_back(thread(ShuffleSlicesThr, state + i*slice_length, slice_length,1));
	//Last thread
	Threads.push_back(thread(ShuffleSlicesThr, state + (parallel_degree - 1)*slice_length, slice_length,1));
	//Wrapping up
	for (auto& th : Threads)
		th.join();
	Threads.clear();
#endif


#ifdef _MEASUREINT
	i4 = __rdtscp(&ui4);
	d1 = (i4 - i5) / (state_size);
	printf("ShuffleSlices %2.2f cpb\n", (float)d1 / (16));
#endif
#if defined(KAT) && defined(KAT_INTERNAL)
	fprintf(fp, "\nShuffleSlices:\nBlocks:\n");
	dumpState(fp, state, state_size);
#endif

	
		//Each thread processes a certain number of groups
		size_t distance = state_size / GROUP_SIZE;
		size_t groups_per_thread = distance / (parallel_degree);//Area size of each thread.
		for (unsigned i = 0; i<parallel_degree - 1; ++i)
			Threads.push_back(thread(SubGroupsThr, state + i*groups_per_thread, groups_per_thread, distance));
		//Last thread with possibly larger area
		Threads.push_back(thread(SubGroupsThr, state + (parallel_degree - 1)*groups_per_thread, 
			distance - (parallel_degree - 1)*groups_per_thread, distance));
		for (auto& th : Threads)
			th.join();
		Threads.clear();


#ifdef _MEASUREINT
		i5 = __rdtscp(&ui3);
		d1 = (i5 - i4) / (state_size);
		printf("SubGroups %2.2f cpb\n", (float)d1 / (16));
#endif
#if defined(KAT) && defined(KAT_INTERNAL)
		fprintf(fp, "\nRound %d SubGroups:\nBlocks:\n", l + 1);
		dumpState(fp, state, state_size);
		fprintf(fp, "\n");
#endif

	}

	//5.Finalization

	__m128i a1_intr = _mm_set_epi64x(0,0);
	__m128i a2_intr= _mm_set_epi64x(0, 0);
	for (uint32_t i = 0; i< state_size/2; ++i)
	{
		a1_intr = _mm_xor_si128(a1_intr, state[i]);
		a2_intr = _mm_xor_si128(a2_intr, state[i+state_size / 2]);
	}
	memset(state, 0, state_size + sizeof(__m128i));
	int128 a1(_mm_extract_epi64(a1_intr, 0), _mm_extract_epi64(a1_intr, 1));
	int128 a2(_mm_extract_epi64(a2_intr, 0), _mm_extract_epi64(a2_intr, 1)); 
	uint8_t tag_buffer[32];

	blake2b_init(&BlakeHash, 32);
	blake2b_update(&BlakeHash, (const uint8_t*)&(a1.i0), 8);
	blake2b_update(&BlakeHash, (const uint8_t*)&(a1.i1), 8);
	blake2b_update(&BlakeHash, (const uint8_t*)&(a2.i0), 8);
	blake2b_update(&BlakeHash, (const uint8_t*)&(a2.i1), 8);

	while (outlen > 32)
	{
		blake2b_final(&BlakeHash, tag_buffer, 32);
		memcpy(out, tag_buffer, 32);
		out += 32;
		outlen -= 32;
	}
	blake2b_final(&BlakeHash, tag_buffer, outlen);
	memcpy(out, tag_buffer, outlen);
	memset(tag_buffer, 0, 32);

#ifdef _MEASUREINT
	i6 = __rdtscp(&ui6);
	d1 = (i6 - i5) / (state_size);
	printf("Final phase %2.2f cpb\n", (float)d1 / (16));
	d2 = (i6 - i1) / (state_size);
	printf("Total %2.2f cpb\n", (float)d2 / (16));
#endif
#ifdef KAT
	fprintf(fp, "Tag: ");
	for (unsigned i = 0; i<outlen; ++i)
		fprintf(fp, "%2.2x ", ((unsigned char*)out)[i]);
	fprintf(fp, "\n");
	fclose(fp);
#endif

	delete[] state;
	return 0;
}



int PHS(void *out, size_t outlen, const void *in, size_t inlen, const void *salt, size_t saltlen, 
	uint32_t t_cost, size_t m_cost, uint32_t thread_n=1)
{
	return ArgonFast64Ext((uint8_t*)out, outlen, (const uint8_t*)in, inlen, (const uint8_t*)salt, saltlen, NULL, 0, NULL,0,t_cost, m_cost, thread_n);
}
