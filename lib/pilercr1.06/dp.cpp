#include "pilercr.h"

#undef  TEST_DPREACH
#undef  TEST_TRAPTRIM
#undef  STAT_DPREACH
#undef  SHOW_DIAGS

#undef  REPORT_DPREACH
#undef	REPORT_DPREACH

#define MAXIGAP  5
#define DIFFCOST 3
#define SAMECOST 1

extern Trapezoid **Tarray;
extern int *Covered;
extern int SegMax;

//extern DPHit *SegSols;
//extern int NumSegs;

static int    BLOCKCOST  = DIFFCOST*MAXIGAP;
static int    MATCHCOST  = DIFFCOST+SAMECOST;
static double RMATCHCOST = DIFFCOST+1.;

/*** FORWARD AND REVERSE D.P. EXTENSION ROUTINES ***/
/*      Called at the mid-point of trapezoid -- mid X [lo,hi], the extension
        is computed to an end point and the lowest and highest diagonals
        are recorded.  These are returned in a partially filled DPHit
        record, that will be merged with that returned for extension in the
        opposite direction.
*/


static int  VecMax = -1;
static int *Vec1 = NULL;
static int *Vec2;

DPHit *TraceForwardPath(char *A, int Alen, char *B, int Blen,
                                       int mid, int lo, int hi)
{ static DPHit rez;
  int *V, odd;
  int  mxv, mxl, mxr, mxi, mxj;
  int  i, j;
#ifdef STAT_DPREACH
  long dparea;
#endif


  /* Set basis from (mid,lo) .. (mid,hi) */


  if (lo < 0) lo = 0;
  if (hi > Blen) hi = Blen;


  if ((hi-lo)+MAXIGAP >= VecMax)
    { VecMax = (int) (1.2*((hi-lo) + MAXIGAP + 1) + 10000); 
      Vec1   = (int *) ckrealloc(Vec1,2*VecMax*sizeof(int),"Vector D.P. arrays");
      Vec2   = Vec1 + VecMax;
    }


  V   = Vec1 - lo;
  odd = 1;


  for (j = lo; j <= hi; j++)
    V[j] = 0;


  hi += MAXIGAP;
  if (hi > Blen) hi = Blen;


  for (; j <= hi; j++)
    V[j] = V[j-1] - DIFFCOST;


  mxv = 0;
  mxr = mid - lo;
  mxl = mid - hi;
  mxi = mid;
  mxj = lo;


#ifdef STAT_DPREACH
  dparea = (hi-lo+1);
#endif
 
  /* Advance to next row */


    int maxi = mid + g_MaxDPHitLen - 1;
	if (maxi >= Alen)
		maxi = Alen;
  for (i = mid; lo <= hi && i < maxi; i++)

//  for (i = mid; lo <= hi && i < Alen; i++)
    { int  c, v;
      int *W;


      W = V;
      if (odd)
        V = Vec2 - lo;
      else
        V = Vec1 - lo;
      odd = 1-odd;


#ifdef TEST_DPREACH
      Log("\n%d [%d,%d] x = %d\n      ",i,lo,hi,mxv);
      for (j = lo; j <= hi; j++)
        Log("  %c ",B[j]);
      Log("\n ");
      for (j = lo; j <= hi; j++)
        Log(" %3d",W[j]);
      Log("\n%c",A[i]);
#endif


      v = W[lo];
      c = V[lo] = v - DIFFCOST;
#ifdef TEST_DPREACH
      Log(" %3d",c);
#endif
      for (j = lo+1; j <= hi; j++)
        { int r, t;


          t = c;
          c = v;
          v = W[j];
          if (A[i] == B[j-1] && CharToLetter[(unsigned char) (A[i])] >= 0) c += MATCHCOST;


          r = c;
          if (v > r) r = v;
          if (t > r) r = t;


          V[j] = c = r - DIFFCOST;
          if (c >= mxv)
            { mxv = c;
              mxi = i+1;
              mxj = j;
            }
#ifdef TEST_DPREACH
          Log(" %3d",c);
#endif
        }


      if (j <= Blen)
        { int r;


          if (A[i] == B[j-1] && CharToLetter[(unsigned char) (A[i])] >= 0) v += MATCHCOST;


          r = v;
          if (c > r) r = c;


          V[j] = v = r - DIFFCOST;
          if (v > mxv)
            { mxv = v;
              mxi = i+1;
              mxj = j;
            }
#ifdef TEST_DPREACH
          Log(" %3d",v);
#endif
     
          for (j++; j <= Blen; j++)
            { v -= DIFFCOST;
              if (v < mxv - BLOCKCOST) break;
              V[j] = v;
#ifdef TEST_DPREACH
              Log(" %3d",v);
#endif
            }
        }
#ifdef TEST_DPREACH
      Log("\n");
#endif


      hi = j-1;


      while (lo <= hi && V[lo] < mxv - BLOCKCOST)
        lo += 1;
      while (lo <= hi && V[hi] < mxv - BLOCKCOST)
        hi -= 1;


      if ((hi-lo)+2 > VecMax)
        { VecMax = (int) (1.2*((hi-lo) + 2) + 10000); 
          Vec1   = (int *) ckrealloc(Vec1,2*VecMax*sizeof(int),"Vector D.P. arrays");
          Vec2   = Vec1 + VecMax;
        }


      if ((i+1) - lo > mxr)
        mxr = (i+1) - lo;
      if ((i+1) - hi < mxl)
        mxl = (i+1) - hi;


#ifdef STAT_DPREACH
      dparea += (hi-lo+1);
#endif
    }


#ifdef STAT_DPREACH
   Log("  DP_Area = %ld  Peak is %d @ (%d,%d) in [%d,%d]\n",
          dparea,mxv,mxi,mxj,mxl,mxr);
#endif


  rez.aHi = mxj;
  rez.bHi = mxi;
  rez.ldiag = mxl;
  rez.hdiag = mxr;
  rez.score = mxv;
  return (&rez);
}


DPHit *TraceReversePath(char *A, int Alen, char *B, int Blen,
                                       int top, int lo, int hi, int bot,
                                       int xfactor)
{ static DPHit rez;
  int *V, odd;
  int  mxv, mxl, mxr, mxi, mxj;
  int  i, j;
#ifdef STAT_DPREACH
  long dparea;
#endif


  /* Set basis from (top,lo) .. (top,hi) */


  if (lo < 0)    lo = 0;
  if (hi > Blen) hi = Blen;


  if ((hi-lo)+MAXIGAP >= VecMax)
    { VecMax = (int) (1.2*((hi-lo) + MAXIGAP + 1) + 10000); 
      Vec1   = (int *) ckrealloc(Vec1,2*VecMax*sizeof(int),"Vector D.P. arrays");
      Vec2   = Vec1 + VecMax;
    }


  V   = Vec1 + ((VecMax-1) - hi);
  odd = 1;


  for (j = hi; j >= lo; j--)
    V[j] = 0;


  lo -= MAXIGAP;
  if (lo < 0) lo = 0;


  for (; j >= lo; j--)
    V[j] = V[j+1] - DIFFCOST;


  mxv = 0;
  mxr = top - lo;
  mxl = top - hi;
  mxi = top;
  mxj = lo;


#ifdef STAT_DPREACH
  dparea = (hi-lo+1);
#endif
 
  /* Advance to next row */


  if (top-1 <= bot) xfactor = BLOCKCOST;


    int mini = top - g_MaxDPHitLen;
	if (mini < 0)
		mini = 0;

  for (i = top-1; lo <= hi && i >= mini; i--)
//  for (i = top-1; lo <= hi && i >= 0; i--)
    { int  c, v;
      int *W;


      W = V;
      if (odd)
        V = Vec2 + ((VecMax-1) - hi);
      else
        V = Vec1 + ((VecMax-1) - hi);
      odd = 1-odd;


#ifdef TEST_DPREACH
      Log("\n%d [%d,%d] x = %d\n          ",i,lo,hi,mxv);
      for (j = hi; j >= lo; j--)
        Log("  %c ",B[j-1]);
      Log("\n ");
      for (j = hi; j >= lo; j--)
        Log(" %3d",W[j]);
      Log("\n%c",A[i]);
#endif


      v = W[hi];
      c = V[hi] = v - DIFFCOST;
#ifdef TEST_DPREACH
      Log(" %3d",c);
#endif
      for (j = hi-1; j >= lo; j--)
        { int r, t;


          t = c;
          c = v;
          v = W[j];
          if (A[i] == B[j] && CharToLetter[(unsigned char) (A[i])] >= 0) c += MATCHCOST;


          r = c;
          if (v > r) r = v;
          if (t > r) r = t;


          V[j] = c = r - DIFFCOST;
          if (c >= mxv)
            { mxv = c;
              mxi = i;
              mxj = j;
            }
#ifdef TEST_DPREACH
          Log(" %3d",c);
#endif
        }


      if (j >= 0)
        { int r;


          if (A[i] == B[j] && CharToLetter[(unsigned char) (A[i])] >= 0) v += MATCHCOST;


          r = v;
          if (c > r) r = c;


          V[j] = v = r - DIFFCOST;
          if (v > mxv)
            { mxv = v;
              mxi = i;
              mxj = j;
            }
#ifdef TEST_DPREACH
          Log(" %3d",v);
#endif
     
          for (j--; j >= 0; j--)
            { v -= DIFFCOST;
              if (v < mxv - xfactor) break;
              V[j] = v;
#ifdef TEST_DPREACH
              Log(" %3d",v);
#endif
            }
        }
#ifdef TEST_DPREACH
      Log("\n");
#endif


      lo = j+1;


      while (lo <= hi && V[lo] < mxv - xfactor)
        lo += 1;
      while (lo <= hi && V[hi] < mxv - xfactor)
        hi -= 1;


      if (i == bot) xfactor = BLOCKCOST;


      if ((hi-lo)+2 > VecMax)
        { VecMax = (int) (1.2*((hi-lo) + 2) + 10000); 
          Vec1   = (int *) ckrealloc(Vec1,2*VecMax*sizeof(int),"Vector D.P. arrays");
          Vec2   = Vec1 + VecMax;
        }


      if (i-lo > mxr)
        mxr = i-lo;
      if (i-hi < mxl)
        mxl = i-hi;


#ifdef STAT_DPREACH
      dparea += (hi-lo+1);
#endif
    }


#ifdef STAT_DPREACH
   Log("  DP_Area = %ld  Peak is %d @ (%d,%d) in [%d,%d]\n",
          dparea,mxv,mxi,mxj,mxl,mxr);
#endif


  rez.aLo = mxj;
  rez.bLo = mxi;
  rez.ldiag = mxl;
  rez.hdiag = mxr;
  rez.score = mxv;
  return (&rez);
}

/*** FINDING ALIGNMENTS WITHIN A TRAPEZOIDAL ZONE ***/

int TSORT(const void *l, const void *r)
{ Trapezoid *x, *y;
  x = *((Trapezoid **) l);
  y = *((Trapezoid **) r);
  return (x->bot - y->bot);
}

int StSORT(const void *l, const void *r)
{ DPHit *x, *y;
  x = (DPHit *) l;
  y = (DPHit *) r;
  if (x->aLo < y->aLo)
    return (-1);
  else if (x->aLo > y->aLo)
    return (1);
  else
    return (x->bLo - y->bLo);
}

int FnSORT(const void *l, const void *r)
{ DPHit *x, *y;
  x = (DPHit *) l;
  y = (DPHit *) r;
  if (x->aHi < y->aHi)
    return (-1);
  else if (x->aHi > y->aHi)
    return (1);
  else
    return (x->bHi - y->bHi);
}


#ifdef REPORT_DPREACH
static int  Al_depth;
#endif


void Align_Recursion(char *A, int Alen, char *B, int Blen,
                            Trapezoid *b, int current, int comp,
                            int MinLen, double MaxDiff, int Traplen)
{ int j, mid, indel;
  float pcnt;
  DPHit *hend, *lend;
  Trapezoid ltrp, htrp;


  mid = (b->bot + b->top) / 2;


#ifdef REPORT_DPREACH
  Log(" [%d,%d]x[%d,%d] = %d (Depth = %d)\n",
         b->bot,b->top,b->lft,b->rgt,b->top - b->bot + 1,Al_depth);
#endif


  lend = TraceForwardPath(B,Blen,A,Alen,mid,mid-b->rgt,mid-b->lft);


  { int x;


    x = 0;
    do
      { x += 1;
        hend = TraceReversePath(B,Blen,A,Alen,
                                lend->bHi,lend->aHi,lend->aHi,
                                mid+MAXIGAP,BLOCKCOST+2*x*DIFFCOST);
      }
    while (hend->bLo > mid + x*MAXIGAP && hend->score < lend->score);
  }


  hend->aHi = lend->aHi;
  hend->bHi = lend->bHi;


#ifdef REPORT_DPREACH
  Log("  Got [%d,%d]x[%d,%d] ([%d,%d]) at score = %d\n",
         hend->bLo,hend->bHi,hend->ldiag,hend->hdiag,hend->aLo,hend->aHi,hend->score);
#endif


  ltrp = htrp = *b;
  ltrp.top = hend->bLo - MAXIGAP;
  htrp.bot = hend->bHi + MAXIGAP;


  if (hend->bHi - hend->bLo >= MinLen &&
      hend->aHi - hend->aLo >= MinLen   )


    { indel = abs( (hend->aLo - hend->bLo)
                 - (hend->aHi - hend->bHi) );
      pcnt = (float)((1/RMATCHCOST)
           - (hend->score - indel)
           / (RMATCHCOST*(hend->bHi - hend->bLo)));


      if (pcnt <= MaxDiff)
    
        { hend->error = pcnt;
    
          for (j = current+1; j < Traplen; j++)
            { Trapezoid *t;
              int   ta, tb, ua, ub; 
        
              t = Tarray[j];
              if (t->bot >= hend->bHi) break;
        
              tb = t->top - t->bot + 1;
              ta = t->rgt - t->lft + 1;
              if (t->lft < hend->ldiag)
                ua = hend->ldiag;
              else
                ua = t->lft;
              if (t->rgt > hend->hdiag)
                ub = hend->hdiag;
              else
                ub = t->rgt;
        
              if (ua > ub) continue;
        
              ua = ub - ua + 1;
              if (t->top > hend->bHi)
                ub = hend->bHi - t->bot + 1;
              else
                ub = tb;
        
              if (((1.*ua)/ta)*((1.*ub)/tb) > .99)
                Covered[j] = 1;
            }
        
          //if (NumSegs >= SegMax)
          //  { SegMax = (int)(1.2*NumSegs + 500);
          //    SegSols = (DPHit *) ckrealloc(SegSols,
          //                                          sizeof(DPHit)*SegMax,
          //                                          "Segment Alignment array");
          //  }
        
          { int d;
        
            d = hend->ldiag;  /*  Oops, diags to this point are b-a, not a-b. */
            hend->ldiag = - (hend->hdiag);
            hend->hdiag = - d;
            if (comp)
              { hend->bLo = Blen - hend->bLo;
                hend->bHi = Blen - hend->bHi;
                hend->ldiag = Blen + hend->ldiag;
                hend->hdiag = Blen + hend->hdiag;
              }
          }

		if (AcceptHit(*hend))
	        g_Hits.push_back(*hend);
		else
			{
			if (fDiscardedHits)
				{
				WriteDPHit(fDiscardedHits, *hend, false);
				++g_DiscardedHitCount;
				}
			}
        
#ifdef REPORT_DPREACH
          Log("  Hit from (%d,%d) to (%d,%d) within [%d,%d] score %d\n",
                 hend->aLo,hend->bLo,hend->aHi,hend->bHi,
                 hend->ldiag,hend->hdiag,hend->score);
#endif
        }
    }


#ifdef REPORT_DPREACH
  Al_depth += 1;
#endif
  if (ltrp.top - ltrp.bot > MinLen && ltrp.top < b->top - MAXIGAP)
    Align_Recursion(A,Alen,B,Blen,&ltrp,current,comp,MinLen,MaxDiff,Traplen);
  if (htrp.top - htrp.bot > MinLen)
    Align_Recursion(A,Alen,B,Blen,&htrp,current,comp,MinLen,MaxDiff,Traplen);
#ifdef REPORT_DPREACH
  Al_depth -= 1;
#endif
}
