#include "types.h"

unsigned char CompMap[256];

static bool InitMap()
	{
	for (unsigned i = 0; i < 256; ++i)
		CompMap[i] = i;

	CompMap[(int) 'a'] = 't';
	CompMap[(int) 'c'] = 'g';
	CompMap[(int) 'g'] = 'c';
	CompMap[(int) 't'] = 'a';

	CompMap[(int) 'A'] = 'T';
	CompMap[(int) 'C'] = 'G';
	CompMap[(int) 'G'] = 'C';
	CompMap[(int) 'T'] = 'A';

	return true;
	}

static bool InitDone = InitMap();

void RevComp(char *S, unsigned L)
	{
	unsigned char *s = (unsigned char *) S;
	unsigned char *t = (unsigned char *) (S + L - 1);
	while (s < t)
		{	
		unsigned char c = *s;
		*s++ = CompMap[*t];
		*t-- = CompMap[c];
		}
	if (s == t)
		*s = CompMap[*s];
	}
