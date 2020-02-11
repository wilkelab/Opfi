#include "pilercr.h"
#include "multaln.h"
#include "sw.h"
#include <math.h>
#include <algorithm>

#define TRACE	0

void SWTraceBack(const unsigned char *A, unsigned LA, const unsigned char *B,
  unsigned LB, unsigned Topi, unsigned Topj, score_t TopScore, const score_t *DPM_,
  const score_t *DPD_, const score_t *DPI_, std::string &Path, unsigned *ptrStartA,
  unsigned *ptrStartB)
	{
#if	TRACE
	Log("\n");
	Log("SWTraceBack LengthA=%u LengthB=%u\n", LA, LB);
#endif
	assert(LB > 0 && LA > 0);

	Path.clear();

	unsigned i = Topi;
	unsigned j = Topj;

	const unsigned PCA = LA + 1;

	char cEdgeType = 'M';
#if	TRACE
	Log("TraceBack\n");
	Log("  i    j  Edge\n");
	Log("---  ---  ----\n");
#endif

	for (;;)
		{
		if (i == 0 || j == 0)
			goto Done;

		char cPrevEdgeType = cEdgeType;
#if	TRACE
		unsigned Prev_i = i;
		unsigned Prev_j = j;
#endif

/***
i is prefix length in A
j is prefix length in B
Current cell is DP<cEdgeType>(i,j).
Determine the predecessor cell.
***/
		switch (cEdgeType)
			{
		case 'M':
			{
			assert(i > 0);
			assert(j > 0);
			const unsigned char cA = A[i-1];
			const unsigned char cB = B[j-1];

			const score_t Score = DPM(i, j);
			const score_t scoreMatch = SUBST(cA, cB);

			score_t scoreMM = MINUS_INF;
			score_t scoreDM = MINUS_INF;
			score_t scoreIM = MINUS_INF;
			if (i > 0 && j > 0)
				scoreMM = DPM(i-1, j-1) + scoreMatch;
			if (i > 1)
				scoreDM = DPD(i-1, j-1) + scoreMatch;
			if (j > 1)
				scoreIM = DPI(i-1, j-1) + scoreMatch;

			if (EQ(scoreMM, Score))
				cEdgeType = 'M';
			else if (EQ(scoreDM, Score))
				cEdgeType = 'D';
			else if (EQ(scoreIM, Score))
				cEdgeType = 'I';
			else if (EQ(0, Score))
				goto Done;
			else
				Quit("TraceBack: failed to match M");

			--i;
			--j;
			break;
			}

		case 'D':
			{
			assert(i > 0);
			const score_t Score = DPD(i, j);

			score_t scoreMD = MINUS_INF;
			score_t scoreDD = MINUS_INF;
			if (i > 1)
				{
				scoreMD = DPM(i-1, j) + GAPOPEN;
				scoreDD = DPD(i-1, j) + GAPEXTEND;
				}

			if (EQ(Score, scoreMD))
				cEdgeType = 'M';
			else if (EQ(Score, scoreDD))
				cEdgeType = 'D';
			else
				Quit("TraceBack: failed to match D");

			--i;
			break;
			}

		case 'I':
			{
			assert(j > 0);
			const score_t Score = DPI(i, j);

			score_t scoreMI = MINUS_INF;
			score_t scoreII = MINUS_INF;
			if (j > 1)
				{
				scoreMI = DPM(i, j-1) + GAPOPEN;
				scoreII = DPI(i, j-1) + GAPEXTEND;
				}

			if (EQ(Score, scoreMI))
				cEdgeType = 'M';
			else if (EQ(Score, scoreII))
				cEdgeType = 'I';
			else
				Quit("TraceBack: failed to match I");

			--j;
			break;
			}

		default:
			assert(false);
			}

		Path.push_back(cPrevEdgeType);

#if TRACE
		{
		char cA = (cPrevEdgeType == 'M' || cPrevEdgeType == 'D') ? A[Prev_i-1] : '-';
		char cB = (cPrevEdgeType == 'M' || cPrevEdgeType == 'I') ? B[Prev_j-1] : '-';
		Log("%3d  %3d  %4c  %c%c\n", Prev_i-1, Prev_j-1, cPrevEdgeType, cA, cB);
		}
#endif
		}

Done:
#if	TRACE
	Log("StartA=%d StartB=%d\n", i, j);
#endif
	*ptrStartA = i;
	*ptrStartB = j;
	std::reverse(Path.begin(), Path.end());

#if	TRACE
	{
	score_t sp = ScorePathLocal(A, LA, B, LB, *ptrStartA, *ptrStartB, Path.c_str());
	if (!EQ(sp, TopScore))
		{
		Log("Path=%s\n", Path.c_str());
		Quit("TraceBackSW: scoreMax=%d != ScorePath=%d", TopScore, sp);
		}
	}
#endif
	}
