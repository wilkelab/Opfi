#ifndef MULTALN_H
#define MULTALN_H

#ifdef _MSC_VER
#pragma warning(disable: 4996)	// deprecated functions
#define _CRT_SECURE_NO_DEPRECATE	1
#endif

#include <vector>
#include <memory.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include "types.h"
#include "params.h"
#include "myassert.h"

#ifndef _MSC_VER
#define stricmp strcasecmp
#define strnicmp strncasecmp
#define	_snprintf snprintf
#define _fsopen(name, mode, share)	fopen((name), (mode))
#endif

enum
	{
	BIT_MM = 0x00,
	BIT_DM = 0x01,
	BIT_IM = 0x02,
	BIT_xM = 0x03,

	BIT_DD = 0x00,
	BIT_MD = 0x04,
	//  ID not allowed
	BIT_xD = 0x04,

	BIT_II = 0x00,
	BIT_MI = 0x08,
	//  DI not allowed
	BIT_xI = 0x08,
	};

// NX=Nucleotide alphabet
enum NX
	{
	NX_A,
	NX_C,
	NX_G,
	NX_T,
	NX_U = NX_T,
	NX_N,
	NX_GAP
	};

const size_t MAX_ALPHA = 4;
const size_t MAX_ALPHA_EX = 6;
const size_t g_AlphaSize = 4;

typedef float BASETYPE;
typedef BASETYPE FCOUNT;
typedef BASETYPE SCORE;

typedef float SCOREMATRIX[32][32];
typedef SCOREMATRIX *PTR_SCOREMATRIX;
extern PTR_SCOREMATRIX g_ptrScoreMatrix;

static inline bool BTEq2(BASETYPE b1, BASETYPE b2)
	{
	double diff = fabs(b1 - b2);
	if (diff < 0.0001)
		return true;
	double sum = fabs(b1) + fabs(b2);
	return diff/sum < 0.005;
	}

static inline bool BTEq(double b1, double b2)
	{
	return BTEq2((BASETYPE) b1, (BASETYPE) b2);
	}

static inline bool ScoreEq(SCORE s1, SCORE s2)
	{
	return BTEq(s1, s2);
	}

extern SCORE g_scoreGapOpen;
extern SCORE g_scoreGapExtend;

const unsigned uInsane = UINT_MAX;
const double dInsane = -9e-9;
const BASETYPE BTInsane = (BASETYPE) -9e-9;
const SCORE MINUS_INFINITY = (SCORE) -1e37;

#include "utils.h"
#include "seq.h"
#include "seqvect.h"
#include "msa.h"
#include "tree.h"
#include "distfunc.h"
#include "pwpath.h"
#include "estring.h"
#include "profile.h"
#include "distcalc.h"

extern unsigned g_CharToLetter[];
extern unsigned g_CharToLetterEx[];

extern char g_LetterToChar[];
extern char g_LetterExToChar[];

extern char g_UnalignChar[];
extern char g_AlignChar[];

extern bool g_IsWildcardChar[];
extern bool g_IsResidueChar[];

#define CharToLetter(c)		(g_CharToLetter[(unsigned char) (c)])
#define CharToLetterEx(c)	(g_CharToLetterEx[(unsigned char) (c)])

#define LetterToChar(u)		(g_LetterToChar[u])
#define LetterExToChar(u)	(g_LetterExToChar[u])

#define IsResidueChar(c)	(g_IsResidueChar[(unsigned char) (c)])
#define IsGapChar(c)		('-' == (c) || '.' == (c))
#define IsWildcardChar(c)	(g_IsWildcardChar[(unsigned char) (c)])

#define AlignChar(c)		(g_AlignChar[(unsigned char) (c)])
#define UnalignChar(c)		(g_UnalignChar[(unsigned char) (c)])

void GetGuideTree(const SeqVect &Seqs, Tree &GuideTree);
void KmerDist(const SeqVect &Seqs, DistFunc &DF);
void ProgressiveAlign(const SeqVect &Seqs, const Tree &GuideTree, MSA &Aln);
SCORE GlobalAlign(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path);
void BitTraceBack(char **TraceBack, unsigned uLengthA, unsigned uLengthB,
  char LastEdge, PWPath &Path);
SCORE AlignProfiles(
  const ProfPos *PA, unsigned uLengthA,
  const ProfPos *PB, unsigned uLengthB,
  PWPath &Path, ProfPos **ptrPout, unsigned *ptruLengthOut);
void AlignTwoProfsGivenPath(const PWPath &Path,
  const ProfPos *PA, unsigned uPrefixLengthA,
  const ProfPos *PB, unsigned uPrefixLengthB,
  ProfPos **ptrPOut, unsigned *ptruLengthOut);
void MakeRootMSA(const SeqVect &v, const Tree &GuideTree, ProgNode Nodes[],
  MSA &Aln);
void MultipleAlign(SeqVect &Seqs, MSA &Aln);
void UPGMA(const DistCalc &DC, Tree &tree);
void GetConsSeq(const MSA &Aln, double MinCons, int *ptrStartCol,
  int *ptrEndCol, Seq &ConsSeq);
void GetConsSymbols(const MSA &Aln, Seq &ConsSymbols);

#endif	// MULTALN_H
