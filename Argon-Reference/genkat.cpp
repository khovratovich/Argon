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

void GenKat()
{
	unsigned char out[128];
	unsigned char zero_array[256];
	memset(zero_array, 0, 256);
	unsigned char one_array[256];
	memset(one_array, 1, 256);
	unsigned t_cost = 3;
	//unsigned m_cost = 2;
#ifdef KAT
	remove("kat.log");
#endif
	for (unsigned m_cost = 1; m_cost <= 10; m_cost *= 10)
	{

		for (unsigned p_len = 0; p_len < 256; p_len += 128)
		{
			for (unsigned s_len = 8; s_len <= 24; s_len += 16)
			{
				for (unsigned thr = 1; thr <= 4; ++thr)
				{
					for (uint16_t outlen = 8; outlen <= 128; outlen *= 4)
					{

#ifdef _MEASURE
						uint64_t  i2, i3, d2;
						uint32_t ui2, ui3;
#endif


						
#ifdef _MEASURE
						clock_t start = clock();
						i2 = __rdtscp(&ui2);
#endif
						ArgonRef(out, outlen, zero_array, p_len, one_array, s_len, NULL, 0, NULL, 0, t_cost, m_cost, thr);

#ifdef _MEASURE
						i3 = __rdtscp(&ui3);
						clock_t finish = clock();

						d2 = (i3 - i2) / (m_cost);
						float mcycles = (float)(i3 - i2) / (1 << 20);
						printf("Argon:  %d iterations %2.2f cpb %2.2f Mcycles\n", t_cost, (float)d2 / 1000, mcycles);

						printf("Tag: ");
						for (unsigned i = 0; i < outlen; ++i)
							printf("%2.2x ", ((unsigned char*)out)[i]);
						printf("\n");

						float run_time = ((float)finish - start) / (CLOCKS_PER_SEC);
						printf("%2.4f seconds\n", run_time);
#endif
					}
				}
			}
		}
	}
}

int main(int argc, char* argv[])
{
	GenKat();
}