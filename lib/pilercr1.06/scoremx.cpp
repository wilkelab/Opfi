#include "multaln.h"

const float NUC_OPEN = -400;
const float NUC_EXTEND = -25;
const float NUC_SP_CENTER = 2*NUC_EXTEND;

float g_scoreGapExtend = NUC_EXTEND;
float g_scoreGapOpen = NUC_OPEN;

#define v(x)	((float) x + NUC_SP_CENTER)
#define ROW(A, C, G, T) \
	{ v(A), v(C), v(G), v(T) },

float NUC_SP[32][32] =
	{
//         A        C        G        T
ROW(     100,     -50,     -50,     -50) // A

ROW(     -50,     100,     -50,     -50) // C

ROW(     -50,     -50,     100,     -50) // G

ROW(     -50,     -50,     -50,     100) // T
	};

PTR_SCOREMATRIX g_ptrScoreMatrix = &NUC_SP;

void LogMx()
	{
	for (int i = 0; i < 32; ++i)
		{
		for (int j = 0; j < 32; ++j)
			{
			Log(" %5.1f", NUC_SP[i][j]);
			}
		Log("\n");
		}
	}
