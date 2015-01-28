#include "stdio.h"

#include "wmmintrin.h"
#include <immintrin.h> 
#include <intrin.h>
#include <stdint.h>
#include <time.h> 

#include <string>
using namespace std;

#include "blake2.h"
#include "data.h"









void SubGroups(int128* state, unsigned width)
{
	for(unsigned i=0; i<width/32; i++)
	{
		//Computing X_i:
		int128 X[16];
		int128 Input[32];
		for (uint8_t j = 0; j < 32; ++j)
			Input[j] = state[i + j*(width / 32)];

		X[ 0] =  Input[ 3]^ Input[ 7]^ Input[11]^ Input[15]^ Input[19]^ Input[23]^ Input[27]^ Input[31];
		X[ 1] =  Input[ 1]^ Input[ 3]^ Input[ 9]^ Input[11]^ Input[17]^ Input[19]^ Input[25]^ Input[27];
		X[ 2] =  Input[ 0]^ Input[ 2]^ Input[ 4]^ Input[ 6]^ Input[16]^ Input[18]^ Input[20]^ Input[22];
		X[ 3] =  Input[ 1]^ Input[ 3]^ Input[ 5]^ Input[ 7]^ Input[ 9]^ Input[11]^ Input[13]^ Input[15];
		X[ 4] =  Input[ 6]^ Input[ 7]^ Input[14]^ Input[15]^ Input[22]^ Input[23]^ Input[30]^ Input[31];
		X[ 5] =  Input[10]^ Input[11]^ Input[14]^ Input[15]^ Input[26]^ Input[27]^ Input[30]^ Input[31];
		X[ 6] =  Input[16]^ Input[17]^ Input[20]^ Input[21]^ Input[24]^ Input[25]^ Input[28]^ Input[29];
		X[ 7] =  Input[12]^ Input[13]^ Input[14]^ Input[15]^ Input[28]^ Input[29]^ Input[30]^ Input[31];
		X[ 8] =  Input[ 4]^ Input[ 5]^ Input[ 6]^ Input[ 7]^ Input[12]^ Input[13]^ Input[14]^ Input[15];
		X[ 9] =  Input[16]^ Input[17]^ Input[18]^ Input[19]^ Input[20]^ Input[21]^ Input[22]^ Input[23];
		X[10] =  Input[ 1]^ Input[ 5]^ Input[ 9]^ Input[13]^ Input[17]^ Input[21]^ Input[25]^ Input[29];
		X[11] =  Input[ 2]^ Input[ 6]^ Input[10]^ Input[14]^ Input[18]^ Input[22]^ Input[26]^ Input[30];
		X[12] =  Input[ 4]^ Input[ 5]^ Input[ 6]^ Input[ 7]^ Input[20]^ Input[21]^ Input[22]^ Input[23];
		X[13] =  Input[ 8]^ Input[ 9]^ Input[10]^ Input[11]^ Input[24]^ Input[25]^ Input[26]^ Input[27];
		X[14] =  Input[ 0]^ Input[ 1]^ Input[ 2]^ Input[ 3]^ Input[ 8]^ Input[ 9]^ Input[10]^ Input[11];
		X[15] =  Input[ 0]^ Input[ 4]^ Input[ 8]^ Input[12]^ Input[16]^ Input[20]^ Input[24]^ Input[28];
		
		

		for(unsigned j=0; j<16;++j)
		{
			AES_reduced(X[j]);//Computing F's
			state[i + 2 * j*(width / 32)] ^= X[j]; //XORs
			state[i + (2 * j + 1)*(width / 32)] ^= X[j];
			AES_reduced(state[i + 2 * j*(width / 32)]);
			AES_reduced(state[i + (2 * j + 1)*(width / 32)]);
		}
	}
}

void ShuffleSlices(int128* state, uint32_t  area, uint32_t slice_n)
{
	uint32_t slice_length = area / slice_n;
	for (uint8_t s = 0; s<slice_n; ++s) //Loop on slices
	{
		uint32_t j = slice_length - 1;
		for (uint32_t i = 0; i<slice_length; ++i)
		{
			//j <- S[j]+ S[i]
			//Swap(S[i],S[j])
			uint32_t index1 = s*slice_length + i;
			uint32_t index2 = s*slice_length + j;
			uint64_t value = state[index2].i0 + state[index1].i0;
			j = value % slice_length;
			swap(state[index1],state[index2]);
		}
	}
}

int ArgonRef(uint8_t *out, uint32_t outlen, const uint8_t *pwd, uint32_t pwdlen, const uint8_t *salt, uint32_t saltlen, const uint8_t *secret,
	uint8_t secretlen, const void *ad, uint32_t adlen, uint32_t t_cost, uint32_t m_cost, uint8_t parallel_degree)
{
	Init();  //Initializing Galois field multiplication table.
	int128* state;  //Array A of blocks

	//0. Restricting parameters
	if(outlen>MAX_OUTLEN)
		outlen=MAX_OUTLEN;
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
	if(m_cost<MIN_MEMORY)
		m_cost = MIN_MEMORY;
	if(m_cost>MAX_MEMORY)
		m_cost = MAX_MEMORY;
	
	//minimum t_cost =3
	if(t_cost<MIN_TIME)
		t_cost = MIN_TIME;

	if (parallel_degree<MIN_SLICE)
		parallel_degree = MIN_SLICE;
	if (parallel_degree>MAX_SLICE)
		parallel_degree = MAX_SLICE;
	
	
#ifdef KAT
	FILE* fp=fopen("kat.log","a+");
	fprintf(fp,"=======================================\n");
	fprintf(fp,"Iterations: %d, Memory: %d KBytes, Parallelism: %d slices, Tag length: %d bytes\n", t_cost, m_cost,parallel_degree, outlen);
	fprintf(fp,"Password: ");
	for(unsigned i=0; i<pwdlen; ++i)
		fprintf(fp,"%2.2x ",((uint8_t*)pwd)[i]);
	fprintf(fp,"\n");
	fprintf(fp,"Salt: ");
	for(unsigned i=0; i<saltlen; ++i)
		fprintf(fp,"%2.2x ",((uint8_t*)salt)[i]);
	fprintf(fp,"\n");
		
#endif
	//Initial hashing

	uint64_t input_hash[4];

	blake2b_state BlakeHash;
	blake2b_init(&BlakeHash,32);

	blake2b_update(&BlakeHash, (const uint8_t*)&pwdlen, sizeof(pwdlen));
	blake2b_update(&BlakeHash, (const uint8_t*)pwd,pwdlen);
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

#if defined KATINT
	fprintf(fp, "Input hash:\n");
	for (unsigned i = 0; i<4; ++i)
	{
		fprintf(fp, "%.16llx ", input_hash[i]);
	}
	fprintf(fp, "\n");

#endif

	//2. Filling blocks
	uint32_t state_size = m_cost*64;
	state = new int128[state_size];
	if(state==NULL)
		return 1;
	printf("Memory allocated: %d KBytes\n",state_size*sizeof(int128)/(1<<10));

	for(unsigned i=0; i<state_size; ++i)
	{
		//Input part
		state[i].i0 = input_hash[i%4];
		//Counter
		state[i].i1 ^= ((uint64_t)i);
	}
	memset(input_hash, 0, 4*sizeof(uint64_t));
	blake2b_init(&BlakeHash, 32);

	//3. Initial transformation
	for(unsigned i=0; i<state_size; ++i)
	{
		AES_reduced(state[i]);
	}



	//4. Rounds: 
	for(unsigned l=0; l <t_cost; ++l)
	{
		SubGroups(state,state_size);
#ifdef KATINT
	fprintf(fp,"Round %d SubGroups:\nBlocks:\n",l+1);
	for(unsigned i=0; i<state_size; ++i)
	{
		fprintf(fp,"Block %3.3d: H: %.16llx L: %.16llx",i,state[i].i1,state[i].i0);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
		
#endif
		ShuffleSlices(state,state_size,parallel_degree);
			#ifdef KATINT
	fprintf(fp,"ShuffleSlices:\nBlocks:\n");
	for(unsigned i=0; i<state_size; ++i)
	{
		fprintf(fp,"Block %3.3d: H: %.16llx L: %.16llx",i,state[i].i1,state[i].i0);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
		
#endif
	}

	//5.Finalization
	SubGroups(state,state_size);
#ifdef KATINT
	fprintf(fp,"Last round: SubGroups:\nBlocks:\n");
	for(unsigned i=0; i<state_size; ++i)
	{
		fprintf(fp,"Block %3.3d: H: %.16llx L: %.16llx",i,state[i].i1,state[i].i0);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
		
#endif


	int128 a1(0,0);
	int128 a2(0,0);
	for(unsigned i=0; i< state_size/2; ++i)
	{
		a1 ^= state[i];
		a2 ^= state[i+state_size/2];
		state[i] = int128(0,0);
		state[i+state_size/2] = int128(0,0);
	}

	uint8_t tag_buffer[32];

	blake2b_init(&BlakeHash, 32);
	blake2b_update(&BlakeHash, (const uint8_t*)&(a1.i0), 8);
	blake2b_update(&BlakeHash, (const uint8_t*)&(a1.i1), 8); 
	blake2b_update(&BlakeHash, (const uint8_t*)&(a2.i0), 8);
	blake2b_update(&BlakeHash, (const uint8_t*)&(a2.i1), 8);

	uint8_t* out_flex = out;
	uint32_t outlen_flex = outlen;
	while (outlen_flex > 16)
	{
		blake2b_final(&BlakeHash, tag_buffer, 32);
		memcpy(out_flex, tag_buffer, 16);
		out_flex += 16;
		outlen_flex -= 16;
		blake2b_init(&BlakeHash, 32);
		blake2b_update(&BlakeHash, tag_buffer, 32);
	}
	blake2b_final(&BlakeHash, tag_buffer, outlen_flex);
	memcpy(out, tag_buffer, outlen_flex);
	memset(tag_buffer, 0, 32);
#ifdef KAT
	fprintf(fp,"Tag: ");
	for(unsigned i=0; i<outlen; ++i)
		fprintf(fp,"%2.2x ",((uint8_t*)out)[i]);
	fprintf(fp,"\n");
	fclose(fp);
#endif KAT
	

	delete state;
	return 0;
}






int PHS(void *out, size_t outlen, const void *in, size_t inlen, const void *salt, uint32_t  saltlen, 
	uint32_t t_cost, uint32_t m_cost)
{
	return ArgonRef((uint8_t*)out, outlen, (const uint8_t*)in, inlen, (const uint8_t*)salt, saltlen, NULL, 0, NULL,0, t_cost, m_cost,1);
}
