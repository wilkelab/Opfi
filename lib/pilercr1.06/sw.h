#ifndef SW_H
#define SW_H

#include <string>

typedef float score_t;
const score_t MINUS_INF = -999999;

#undef DPM
#undef DPD
#undef DPI

#undef DPMb
#undef DPDb
#undef DPIb

// Macros to simulate 2D matrices
#define DPM(i, j)	DPM_[(j)*PCA + (i)]
#define DPD(i, j)	DPD_[(j)*PCA + (i)]
#define DPI(i, j)	DPI_[(j)*PCA + (i)]

#define DPMb(i, j)	DPMb_[BandIndex((i), (j), LA, LB, r)]
#define DPDb(i, j)	DPDb_[BandIndex((i), (j), LA, LB, r)]
#define DPIb(i, j)	DPIb_[BandIndex((i), (j), LA, LB, r)]

//#define SUBST(cA, cB)		\
//	((*g_ptrScoreMatrix)[CharToLetter[(unsigned char) (cA)]][CharToLetter[(unsigned char) (cB)]])

static inline float SUBST(unsigned char cA, unsigned char cB)
	{
	unsigned LetA = CharToLetter[(unsigned char) cA];
	unsigned LetB = CharToLetter[(unsigned char) cB];
	assert(LetA < 32);
	assert(LetB < 32);
	return (*g_ptrScoreMatrix)[LetA][LetB];
	}

extern float g_scoreGapExtend;
extern float g_scoreGapOpen;

#define GAPOPEN		g_scoreGapOpen
#define GAPEXTEND	g_scoreGapExtend

#define EQ(a, b)	((a) == (b))

score_t SWSimple(const char *A_, unsigned LA, const char *B_,
  unsigned LB, unsigned *ptrStartA, unsigned *ptrStartB, std::string &Path);
class Seq;
score_t SW(const Seq &A, const Seq &B, unsigned *ptrStartA, unsigned *ptrStartB,
  std::string &Path);

#endif	// SW_H
