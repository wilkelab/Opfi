#include "multaln.h"

#define TRACE 0

#define MIN(x, y)	(((x) < (y)) ? (x) : (y))
#define MAX(x, y)	(((x) > (y)) ? (x) : (y))

const unsigned TUPLE_COUNT = 4*4*4*4;
static unsigned char Count1[TUPLE_COUNT];
static unsigned char Count2[TUPLE_COUNT];

// Nucleotide groups according to MAFFT (sextet5)
// 0 =  A
// 1 =  C
// 2 =  G
// 3 =  T
// 4 =  other

static unsigned ResidueGroup[] =
	{
	0,		// NX_A,
	1,		// NX_C,
	2,		// NX_G,
	3,		// NX_T/U
	0,		// NX_N,
	0,		// NX_R,
	0,		// NX_Y,
	0,		// NX_GAP
	};
static unsigned uResidueGroupCount = sizeof(ResidueGroup)/sizeof(ResidueGroup[0]);

static char *TupleToStr(int t)
	{
	static char Letter[4] = { 'A', 'C', 'G', 'T' };
	static char s[7];
	int t1, t2, t3, t4;

	t1 = t%4;
	t2 = (t/4)%4;
	t3 = (t/(4*4))%4;
	t4 = (t/(4*4*4))%4;

	s[3] = Letter[t1];
	s[2] = Letter[t2];
	s[1] = Letter[t3];
	s[0] = Letter[t4];
	return s;
	}

static unsigned GetTuple(const unsigned uLetters[], unsigned n)
	{
	assert(uLetters[n] < uResidueGroupCount);
	assert(uLetters[n+1] < uResidueGroupCount);
	assert(uLetters[n+2] < uResidueGroupCount);
	assert(uLetters[n+3] < uResidueGroupCount);

	unsigned u1 = ResidueGroup[uLetters[n]];
	unsigned u2 = ResidueGroup[uLetters[n+1]];
	unsigned u3 = ResidueGroup[uLetters[n+2]];
	unsigned u4 = ResidueGroup[uLetters[n+3]];

	return u4 + u3*4 + u2*4*4 + u1*4*4*4;
	}

static void CountTuples(const unsigned L[], unsigned uTupleCount, unsigned char Count[])
	{
	memset(Count, 0, TUPLE_COUNT*sizeof(unsigned char));
	for (unsigned n = 0; n < uTupleCount; ++n)
		{
		const unsigned uTuple = GetTuple(L, n);
		++(Count[uTuple]);
		}
	}

static void ListCount(const unsigned char Count[])
	{
	for (unsigned n = 0; n < TUPLE_COUNT; ++n)
		{
		if (0 == Count[n])
			continue;
		Log("%s  %u\n", TupleToStr(n), Count[n]);
		}
	}

void KmerDist(const SeqVect &Seqs, DistFunc &DF)
	{
	const unsigned uSeqCount = Seqs.Length();

	DF.SetCount(uSeqCount);
	if (0 == uSeqCount)
		return;

// Initialize distance matrix to zero
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		DF.SetDist(uSeq1, uSeq1, 0);
		DF.SetId(uSeq1, Seqs[uSeq1]->GetId());
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			DF.SetDist(uSeq1, uSeq2, 0);
		}

// Convert to letters
	unsigned **Letters = new unsigned *[uSeqCount];
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq &s = *(Seqs[uSeqIndex]);
		const unsigned uSeqLength = s.Length();
		unsigned *L = new unsigned[uSeqLength];
		Letters[uSeqIndex] = L;
		for (unsigned n = 0; n < uSeqLength; ++n)
			{
			char c = s[n];
			L[n] = CharToLetterEx(c);
			if (L[n] >= 4)
				L[n] = 3;
			}
		}

	unsigned **uCommonTupleCount = new unsigned *[uSeqCount];
	for (unsigned n = 0; n < uSeqCount; ++n)
		{
		uCommonTupleCount[n] = new unsigned[uSeqCount];
		memset(uCommonTupleCount[n], 0, uSeqCount*sizeof(unsigned));
		}

	const unsigned uPairCount = (uSeqCount*(uSeqCount + 1))/2;
	unsigned uCount = 0;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &seq1 = *(Seqs[uSeq1]);
		const unsigned uSeqLength1 = seq1.Length();
		if (uSeqLength1 < 5)
			continue;

		const unsigned uTupleCount = uSeqLength1 - 5;
		const unsigned *L = Letters[uSeq1];
		CountTuples(L, uTupleCount, Count1);
#if	TRACE
		{
		Log("Seq1=%d\n", uSeq1);
		Log("Groups:\n");
		for (unsigned n = 0; n < uSeqLength1; ++n)
			Log("%u", ResidueGroup[L[n]]);
		Log("\n");

		Log("Tuples:\n");
		ListCount(Count1);
		}
#endif

		for (unsigned uSeq2 = 0; uSeq2 <= uSeq1; ++uSeq2)
			{
			++uCount;
			Seq &seq2 = *(Seqs[uSeq2]);
			const unsigned uSeqLength2 = seq2.Length();
			if (uSeqLength2 < 4)
				{
				if (uSeq1 == uSeq2)
					DF.SetDist(uSeq1, uSeq2, 0);
				else
					DF.SetDist(uSeq1, uSeq2, 1);
				continue;
				}

		// First pass through seq 2 to count tuples
			const unsigned uTupleCount = uSeqLength2 <= 5 ? 0 : uSeqLength2 - 5;
			const unsigned *L = Letters[uSeq2];
			CountTuples(L, uTupleCount, Count2);
#if	TRACE
			Log("Seq2=%d Counts=\n", uSeq2);
			ListCount(Count2);
#endif

		// Second pass to accumulate sum of shared tuples
			unsigned uSum = 0;
			for (unsigned n = 0; n < uTupleCount; ++n)
				{
				const unsigned uTuple = GetTuple(L, n);
				uSum += MIN(Count1[uTuple], Count2[uTuple]);

			// This is a hack to make sure each unique tuple counted only once.
				Count2[uTuple] = 0;
				}
#if	TRACE
			{
			Seq &s1 = *(Seqs[uSeq1]);
			Seq &s2 = *(Seqs[uSeq2]);
			const char *pName1 = s1.GetName();
			const char *pName2 = s2.GetName();
			Log("Common count %s(%d) - %s(%d) =%u\n",
			  pName1, uSeq1, pName2, uSeq2, uSum);
			}
#endif
			uCommonTupleCount[uSeq1][uSeq2] = uSum;
			uCommonTupleCount[uSeq2][uSeq1] = uSum;
			}
		}

	uCount = 0;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &s1 = *(Seqs[uSeq1]);
		const char *pName1 = s1.GetName();

		double dCommonTupleCount11 = uCommonTupleCount[uSeq1][uSeq1];
		if (0 == dCommonTupleCount11)
			dCommonTupleCount11 = 1;

		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			{
			++uCount;

			double dCommonTupleCount22 = uCommonTupleCount[uSeq2][uSeq2];
			if (0 == dCommonTupleCount22)
				dCommonTupleCount22 = 1;

			const double dDist1 = (dCommonTupleCount11 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount11;
			const double dDist2 = (dCommonTupleCount22 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount22;

			const double dMinDist = MIN(dDist1, dDist2);
			DF.SetDist(uSeq1, uSeq2, (float) dMinDist);
#if	TRACE
			Log("%u,%d Common11=%g Common22=%g Common12=%u d1=%g d2=%g\n",
			  uSeq1, uSeq2,
			  dCommonTupleCount11, dCommonTupleCount22,
			  uCommonTupleCount[uSeq1][uSeq2],
			  dDist1, dDist2);
#endif
			}
		}

	for (unsigned n = 0; n < uSeqCount; ++n)
		{
		delete[] uCommonTupleCount[n];
		delete[] Letters[n];
		}
	delete[] uCommonTupleCount;
	delete[] Letters;
	}
//
//void Test()
//	{
//	DistFunc DF;
//	Seq *s1 = new Seq;
//	Seq *s2 = new Seq;
//	s1->FromString("GCATCGCCCGCCAGCAATGGCGGGCGCGGATTGAAAC", "s1");
//	s2->FromString("GCATCGCCCCTCGGCAACGAGGGGCGCGGATTGAAAC", "s2");
//	SeqVect Seqs;
//	Seqs.push_back(s1);
//	Seqs.push_back(s2);
//	KmerDist(Seqs, DF);
//	Log("Dist=%.1f\n", DF.GetDist(0, 1));
//	Quit("Test");
//	}
