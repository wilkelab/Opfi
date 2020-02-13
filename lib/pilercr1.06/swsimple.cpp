#include "pilercr.h"
#include "multaln.h"
#include "sw.h"

// Simple version of Smith-Waterman.
// Could certainly be made smaller & faster if desired.

#define	TRACE	0

#define ScoreToStr	LocalScoreToStr

void SWTraceBack(const unsigned char *A, unsigned LA, const unsigned char *B,
  unsigned LB, unsigned Topi, unsigned Topj, score_t TopScore, const score_t *DPM_,
  const score_t *DPD_, const score_t *DPI_, std::string &Path, unsigned *ptrStartA,
  unsigned *ptrStartB);

static const char *ScoreToStr(score_t s)
	{
	static char str[16];
	if (s <= MINUS_INF/2)
		return "     *";
#if	INTEGER_DP
	sprintf(str, "%6d", s);
#else
	sprintf(str, "%6.0f", s);
#endif
	return str;
	}

#if	TRACE
static void LogDP(const score_t *DPM_, const unsigned char *A, 
  const unsigned char *B, unsigned PCA, unsigned PCB)
	{
	Log("        ");
	for (unsigned j = 0; j < PCB; ++j)
		{
		char c = ' ';
		if (j > 0)
			c = B[j - 1];
		Log(" %4u:%c", j, c);
		}
	Log("\n");
	for (unsigned i = 0; i < PCA; ++i)
		{
		char c = ' ';
		if (i > 0)
			c = A[i - 1];
		Log("%4u:%c  ", i, c);
		for (unsigned j = 0; j < PCB; ++j)
			Log(" %s", ScoreToStr(DPM(i, j)));
		Log("\n");
		}
	}
#endif	// TRACE

score_t SW(const Seq &A, const Seq &B, unsigned *ptrStartA, unsigned *ptrStartB,
  std::string &Path)
	{
	unsigned LA = A.Length();
	unsigned LB = B.Length();

	char *sA = all(char, LA + 1);
	char *sB = all(char, LB + 1);

	A.ToString(sA, LA + 1);
	B.ToString(sB, LB + 1);

	score_t Score = SWSimple(sA, LA, sB, LB, ptrStartA, ptrStartB, Path);

	freemem(sA);
	freemem(sB);

	return Score;
	}

score_t SWSimple(const char *A_, unsigned LA, const char *B_,
  unsigned LB, unsigned *ptrStartA, unsigned *ptrStartB, std::string &Path)
	{
	assert(LB > 0 && LA > 0);

	const unsigned char *A = (const unsigned char *) A_;
	const unsigned char *B = (const unsigned char *) B_;

	Path.clear();

	const unsigned PCA = LA + 1;
	const unsigned PCB = LB + 1;

// Allocate DP matrices
	const unsigned LM = PCA*PCB;
	score_t *DPM_ = all(score_t, LM);
	score_t *DPD_ = all(score_t, LM);
	score_t *DPI_ = all(score_t, LM);

// Boundaries
	for (unsigned j = 0; j < PCB; ++j)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(0, j) = 0;
		DPD(0, j) = 0;
		DPI(0, j) = 0;
		}

	for (unsigned i = 0; i < PCA; ++i)
		{
		DPM(i, 0) = 0;
		DPD(i, 0) = 0;
		DPI(i, 0) = 0;
		}

// ============
// Main DP loop
// ============
	score_t TopScore = -1;
	unsigned Topi = 0;
	unsigned Topj = 0;
	for (unsigned i = 1; i < PCA; ++i)
		{
		const unsigned char cA = A[i-1];

		for (unsigned j = 1; j < PCB; ++j)
			{
			const unsigned char cB = B[j-1];

			{
		// Match M=LetterA+LetterB
			score_t scoreLL = SUBST(cA, cB);

			score_t scoreMM = DPM(i-1, j-1);
			score_t scoreDM = DPD(i-1, j-1);
			score_t scoreIM = DPI(i-1, j-1);

			score_t scoreBest;
			if (scoreMM >= scoreDM && scoreMM >= scoreIM)
				scoreBest = scoreMM;
			else if (scoreDM >= scoreMM && scoreDM >= scoreIM)
				scoreBest = scoreDM;
			else 
				{
				assert(scoreIM >= scoreMM && scoreIM >= scoreDM);
				scoreBest = scoreIM;
				}
			score_t s = scoreBest + scoreLL;
			if (s > 0)
				{
				DPM(i, j) = s;
				if (s > TopScore)
					{
					TopScore = s;
					Topi = i;
					Topj = j;
					}
				}
			else
				DPM(i, j) = 0;
			}

			{
		// Delete D=LetterA+GapB
			score_t scoreMD = DPM(i-1, j) + GAPOPEN;
			score_t scoreDD = DPD(i-1, j) + GAPEXTEND;

			score_t scoreBest;
			if (scoreMD >= scoreDD)
				scoreBest = scoreMD;
			else
				{
				assert(scoreDD >= scoreMD);
				scoreBest = scoreDD;
				}
			if (scoreBest < 0)
				scoreBest = 0;
			DPD(i, j) = scoreBest;
			}

		// Insert I=GapA+LetterB
			{
			score_t scoreMI = DPM(i, j-1) + GAPOPEN;
			score_t scoreII = DPI(i, j-1) + GAPEXTEND;

			score_t scoreBest;
			if (scoreMI >= scoreII)
				scoreBest = scoreMI;
			else 
				{
				assert(scoreII > scoreMI);
				scoreBest = scoreII;
				}
			if (scoreBest < 0)
				scoreBest = 0;
			DPI(i, j) = scoreBest;
			}
			}
		}

#if TRACE
	Log("DPM:\n");
	LogDP(DPM_, A, B, PCA, PCB);
	Log("DPD:\n");
	LogDP(DPD_, A, B, PCA, PCB);
	Log("DPI:\n");
	LogDP(DPI_, A, B, PCA, PCB);
#endif

	if (TopScore < 0)
		{
		*ptrStartA = 0;
		*ptrStartB = 0;
		return 0;
		}

#if	TRACE
	Log("TopScore=%d Topi=%u Topj=%u\n", TopScore, Topi, Topj);
#endif

	SWTraceBack(A, LA, B, LB, Topi, Topj, TopScore, DPM_, DPD_, DPI_, Path,
	  ptrStartA, ptrStartB);

	freemem(DPM_);
	freemem(DPD_);
	freemem(DPI_);

	return TopScore;
	}
