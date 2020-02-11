#include "pilercr.h"
#include <algorithm>

#define TRACE	0

extern void Align_Recursion(char *A, int Alen, char *B, int Blen,
  Trapezoid *b, int current, int comp, int MinLen, double MaxDiff, int Traplen);
extern int TSORT(const void *l, const void *r);
extern int StSORT(const void *l, const void *r);
extern int FnSORT(const void *l, const void *r);

// Sort by Lo
bool Cmp_Lo(const DPHit &Hit1, const DPHit &Hit2)
	{
	if (Hit1.aLo < Hit2.aLo)
		return true;
	return false;
	}

// Sort by Hi
bool Cmp_Hi(const DPHit &Hit1, const DPHit &Hit2)
	{
	if (Hit1.aHi < Hit2.aHi)
		return true;
	return false;
	}

// shared with aligntraps.cpp
Trapezoid **Tarray = NULL;
int *Covered;
int SegMax = -1;

static int fseg;
static int TarMax = -1;

void AlignTraps(char *A, int Alen, char *B, int Blen, Trapezoid *Traplist,
  int Traplen, int comp, const DPParams &DP)
{ 
	Tarray = NULL;
	Covered = 0;
	SegMax = -1;
	fseg = 0;
	TarMax = -1;

	Progress("Initializing trapezoid alignment");
  const int MinLen = DP.g_MinHitLength;
  const double MaxDiff = 1.0 - DP.MinId;

  Trapezoid *b;
  int i;


  if (Traplen >= TarMax)
    { TarMax = (int) (1.2*Traplen + 500);
      Tarray = (Trapezoid **)
               ckrealloc(Tarray,(sizeof(Trapezoid *) + sizeof(int))*TarMax,"Trapezoid array");
      Covered = (int *) (Tarray + TarMax);
    }

// Arbitrary number to avoid thrashing
	g_Hits.reserve(10000);

  i = 0;
  b = Traplist;
  for (i = 0; i < Traplen; i++)
    { Tarray[i] = b;
      Covered[i] = 0;
	  if (b == 0)
		  {
		  Traplen = i;
		  Warning("Internal problem: in AlignTraps, b=0, truncating list");
		  break;
		  }
	  b = b->next;
    }


  qsort(Tarray,Traplen,sizeof(Trapezoid *),TSORT);


#ifdef SHOW_TRAP
  { int i;
    Trapezoid *t;


    for (i = 0; i < Traplen; i++)             
      { t = Tarray[i];
        printf("  [%d,%d] x [%d,%d]\n",t->bot,t->top,t->lft,t->rgt);
      }
  }
#endif


#ifdef REPORT_DPREACH
  Al_depth = 0;
#endif
  g_HitCount = 0;
  fseg = g_HitCount;
  ProgressStart("Aligning trapezoids");
  int Step = Traplen/50;
  if (Step == 0)
	  Step = 1;
  for (i = 0; i < Traplen; i++)
	  {
	  if (i%Step == 0)
		  ProgressStep(i, Traplen);

		if (! Covered[i])
			{ b = Tarray[i];
				if (b->top - b->bot < k)
					continue;
#if	TRACE
			Log("Align trap %d,%d - %d,%d\n", b->lft, b->top, b->rgt, b->top);
#endif
			if (b->lft < g_MinDiag)
				continue;
			Align_Recursion(A,Alen,B,Blen,b,i,comp,MinLen,MaxDiff,Traplen);
			}
	  }
	ProgressDone();

  /* Remove lower scoring segments that begin or end at
       the same point as a higher scoring segment.       */


  Progress("Removing redundant hits");
  if (g_HitCount > fseg)
    { int i, j;


//      qsort(g_Hits+fseg,g_HitCount-fseg,sizeof(DPHit),StSORT);
	std::vector<DPHit>::iterator pFrom = g_Hits.begin();
	std::vector<DPHit>::iterator pTo = g_Hits.begin();

	pFrom += fseg;
	pTo -= fseg;
	std::sort(pFrom, pTo, Cmp_Lo);

      for (i = fseg; i < g_HitCount; i = j)
        { for (j = i+1; j < g_HitCount; j++)
            { if (g_Hits[j].aLo != g_Hits[i].aLo) break;
              if (g_Hits[j].bLo != g_Hits[i].bLo) break;
              if (g_Hits[j].score > g_Hits[i].score)
                { g_Hits[i].score = -1; i = j; }
              else
                g_Hits[j].score = -1;
            }
        }


	std::sort(pFrom, pTo, Cmp_Hi);
//      qsort(g_Hits+fseg,g_HitCount-fseg,sizeof(DPHit),FnSORT);
      for (i = fseg; i < g_HitCount; i = j)
        { for (j = i+1; j < g_HitCount; j++)
            { if (g_Hits[j].aHi != g_Hits[i].aHi) break;
              if (g_Hits[j].bHi != g_Hits[i].bHi) break;
              if (g_Hits[j].score > g_Hits[i].score)
                { g_Hits[i].score = -1; i = j; }
              else
                g_Hits[j].score = -1;
            }
        }


      for (i = fseg; i < g_HitCount; i++)
        if (g_Hits[i].score >= 0)
          g_Hits[fseg++] = g_Hits[i];
      g_HitCount = fseg;
    }

#ifdef REPORT_SIZES
  printf("\n  %9d segments\n",g_HitCount);
  fflush(stdout);
#endif

  free(Tarray);
  Tarray = 0;

  g_HitCount = (int) g_Hits.size();
}

int SumDPLengths(const DPHit *DPHits, int HitCount)
	{
	int Sum = 0;
	for (int i = 0; i < HitCount; ++i)
		{
		const DPHit &Hit = DPHits[i];
		const int Length = Hit.aHi - Hit.aLo;
		if (Length < 0)
			Quit("SumDPLengths, Length < 0");
		Sum += Length;
		}
	return Sum;
	}
