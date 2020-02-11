#include "pilercr.h"

/***
Generate default aligner parameters given:
	minimum hit length.
	minimum fractional identity.
	sequence lengths.
	maximum memory.
***/

const int DEFAULT_LENGTH = 400;
const double DEFAULT_MIN_ID = 0.94;
const double RAM_FRACT = 0.80;

/***
For minimum word length, choose k=4 arbitrarily.
For max, k=16 definitely won't work with 32-bit ints
because 4^16 = 2^32 = 4GB.
k=14 might be OK, but would have to look carefully at
boundary cases, which I haven't done.
k=13 is definitely safe, so set this as upper bound.
***/
const int MIN_WORD_LENGTH = 4;
const int MAX_WORD_LENGTH = 13;
const int MAX_AVG_INDEX_LIST_LENGTH = 10;
const int TUBE_OFFSET_DELTA = 32;

static double FilterMemRequired(int SeqLengthQ, const FilterParams &FP)
	{
	const double Index = 2*sizeof(INDEX_ENTRY)*pow4d(FP.WordSize);

	const int TubeWidth = FP.TubeOffset + FP.SeedDiffs;
	const double MaxActiveTubes = (SeqLengthQ + TubeWidth - 1)/FP.TubeOffset + 1;
	const double Tubes = MaxActiveTubes*sizeof(TubeState);
	return Index + Tubes;
	}

double AvgIndexListLength(int SeqLengthQ, const FilterParams &FP)
	{
	return SeqLengthQ / pow4d(FP.WordSize);
	}

double TotalMemRequired(int SeqLengthQ, const FilterParams &FP)
	{
	const double Filter = FilterMemRequired(SeqLengthQ, FP);
	const double Seq = SeqLengthQ;
	return Filter + Seq;
	}

static void CalcParams(int g_MinHitLength, double MinId, int SeqLengthQ,
  int g_Diameter, double MaxMem, FilterParams *ptrFP, DPParams *ptrDP)
	{
	if (MinId < 0 || MinId > 1.0)
		Quit("CalcParams: bad MinId=%g", MinId);
	if (g_MinHitLength <= 4)
		Quit("CalcParams: bad g_MinHitLength=%d", g_MinHitLength);

// Lower bound on word length k by requiring manageable index.
// Given kmer occurs once every 4^k positions.
// Hence average number of index entries is i = N/(4^k) for random
// string of length N.
// Require i <= I, then k > log_4(N/i).
	const double dSeqLengthA = (double) SeqLengthQ;
	const int MinWordSize = (int) (log4(g_Diameter) - log4(MAX_AVG_INDEX_LIST_LENGTH) + 0.5);

// First choice is that filter criteria are same as DP criteria,
// but this may not be possible.
	int SeedLength = g_MinHitLength;
	int SeedDiffs = (int) (g_MinHitLength*(1.0 - MinId) + 0.5);

// Find filter valid filter parameters,
// starting from preferred case.
	int WordSize = -1;
	for (;;)
		{
		int MinWords = -1;
		for (WordSize = MAX_WORD_LENGTH; WordSize >= MinWordSize; --WordSize)
			{
			ptrFP->WordSize = WordSize;
			ptrFP->SeedLength = SeedLength;
			ptrFP->SeedDiffs = SeedDiffs;
			ptrFP->TubeOffset = ptrFP->SeedDiffs + TUBE_OFFSET_DELTA;

			double Mem = TotalMemRequired(SeqLengthQ, *ptrFP);
			if (MaxMem > 0 && Mem > MaxMem)
				{
				Log("Parameters seedlength=%3d k=%2d seeddiffs=%2d, mem=%.0f Mb > maxmem=%.0f Mb\n",
				  ptrFP->SeedLength,
				  ptrFP->WordSize,
				  ptrFP->SeedDiffs,
				  Mem/1e6,
				  MaxMem/1e6);
				MinWords = -1;
				continue;
				}

			MinWords = MinWordsPerFilterHit(SeedLength, WordSize, SeedDiffs);
			if (MinWords <= 0)
				{
				Log("Parameters seedlength=%3d k=%2d seeddiffs=%2d, B=%d\n",
				  ptrFP->SeedLength,
				  ptrFP->WordSize,
				  ptrFP->SeedDiffs,
				  MinWords);
				MinWords = -1;
				continue;
				}

			const double Len = AvgIndexListLength(g_Diameter, *ptrFP);
			if (Len > MAX_AVG_INDEX_LIST_LENGTH)
				{
				Log("Parameters n=%d k=%d e=%d, B=%d avgixlen=%g > max = %d\n",
				  ptrFP->SeedLength,
				  ptrFP->WordSize,
				  ptrFP->SeedDiffs,
				  MinWords,
				  Len,
				  MAX_AVG_INDEX_LIST_LENGTH);
				MinWords = -1;
				continue;
				}
			break;
			}
		if (MinWords > 0)
			break;

	// Failed to find filter parameters, try
	// fewer errors and shorter seed.
		if (SeedLength >= g_MinHitLength/4)
			{
			SeedLength -= 4;
			continue;
			}

		if (SeedDiffs > 0)
			{
			--SeedDiffs;
			continue;
			}
		
		Quit("Failed to find filter parameters");
		}

	ptrDP->g_MinHitLength = g_MinHitLength;
	ptrDP->MinId = MinId;
	}

// Alignment parameters are specified by length and pctid.
void GetParams(int SeqLengthQ, int g_Diameter, int Length, double MinId,
  FilterParams *ptrFP, DPParams *ptrDP)
	{
	const char *strMaxMem = ValueOpt("maxmem");
	const double MaxMem = (0 == strMaxMem) ? GetRAMSize()*RAM_FRACT : atof(strMaxMem);

	CalcParams(Length, MinId, SeqLengthQ, g_Diameter, MaxMem, ptrFP, ptrDP);
	}
