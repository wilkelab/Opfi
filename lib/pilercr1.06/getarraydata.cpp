#include "pilercr.h"

static bool PileByPos(int PileIndex1, int PileIndex2)
	{
	return GetPile(PileIndex1).Lo < GetPile(PileIndex2).Lo;
	}

bool GetArrayData(const IntVec &ArrayPileIndexes, ArrayData &AD)
	{
	if (ArrayPileIndexes.size() == 0)
		Quit("GetArrayData, size=0");

#if	TRACE
	Log("\n");
	Log("Array\n");
	Log("=====\n");
#endif

	for_CIntVec(ArrayPileIndexes, p)
		AD.PileIndexes.push_back(*p);

// Sort by position
	std::sort(AD.PileIndexes.begin(), AD.PileIndexes.end(), PileByPos);

	SeqVect Seqs;
	int Lo = -1;
	int Hi = -1;
	for_CIntVec(ArrayPileIndexes, p)
		{
		int PileIndex = *p;
		Seq *s = new Seq;
		PileToSeq(g_Piles, PileIndex, *s, &Lo, &Hi);
		if (Seqs.size() == 0)
			AD.Lo = Lo;
		Seqs.push_back(s);
		}

	MSA Aln;
	MultipleAlign(Seqs, Aln);

	int StartCol;
	int EndCol;
	double MinCons = g_DraftMinCons;
	const int SeqCount = (int) Aln.GetSeqCount();
	if (SeqCount >= g_MinOneDiff)
		{
		double AltMinCons = ((double) (SeqCount - 1) / (double) SeqCount) - 0.1;
		if (AltMinCons < MinCons)
			MinCons = AltMinCons;
		}

	GetConsSeq(Aln, MinCons, &StartCol, &EndCol, AD.ConsSeq);
	if ((int) AD.ConsSeq.Length() < g_DraftMinRepeatLength)
		return false;

	char Label[128];
	sprintf(Label, "Cons%d", AD.Id);
	AD.ConsSeq.SetName(Label);
	AD.ConsSeq.SetId(0);
	AD.RepeatLength = AD.ConsSeq.Length();
	AD.SpacerLength = -1;

	if (g_LogAlns)
		{
		Log("\n");
		Log("Draft alignment for array %d:\n", AD.Id);
		Aln.LogMe();
		Log("\n");
		Log("Consensus: ");
		AD.ConsSeq.LogMeSeqOnly();
		Log("\n");
		Log("\n");
		}

#if	TRACE
	Aln.LogMe();
	AD.ConsSeq.LogMe();
	Log("\n");
#endif

	return true;
	}
