#include "pilercr.h"
#include "multaln.h"

#define TRACE		0

void PileToSeq(const std::vector<PileData> &Piles, int PileIndex, Seq &s,
  int *ptrPosLo, int *ptrPosHi)
	{
	const int PileCount = (int) Piles.size();
	assert(PileIndex >= 0 && PileIndex < PileCount);

	const PileData &Pile = Piles[PileIndex];
	const ContigData &Contig = GlobalPosToContig(Pile.Lo);

	int Lo = Pile.Lo - g_HitPadding;
	int Hi = Pile.Hi + g_HitPadding;

	if (Lo < Contig.From)
		Lo = Contig.From;
	int ContigTo = Contig.From + Contig.Length - 1;
	if (Hi >= ContigTo)
		Hi = ContigTo;

	int Length = Hi - Lo + 1;
	assert(Length > 0);

	int n = (int) strlen(Contig.Label) + 64;
	char *Name = all(char, n);
	sprintf(Name, "%s:%d-%d", Contig.Label, Lo, Hi);
	s.SetName(Name);
	freemem(Name);

	s.reserve(Length);
	for (int Pos = Lo; Pos <= Hi; ++Pos)
		s.push_back(g_SeqQ[Pos]);

	*ptrPosLo = Lo;
	*ptrPosHi = Hi;
	}

static void OutputArrayHeaderSim()
	{
	Out("\n\n");
	Out("Array          Sequence    Position      Length  # Copies  Repeat  Spacer  +  Consensus\n");
	Out("=====  ================  ==========  ==========  ========  ======  ======  =  =========\n");
	}

static void OutputArrayHeaderPos()
	{
	Out("\n");
	Out("Array          Sequence    Position      Length  # Copies  Repeat  Spacer    Distance  Consensus\n");
	Out("=====  ================  ==========  ==========  ========  ======  ======  ==========  =========\n");
	}

static void OutputArraySummaryLineSim(const ArrayAln &AA)
	{
	const ContigData &Contig = GlobalPosToContig(AA.Pos);
	const char *Label = Contig.Label;
	int ArrayLength = GetArrayLength(AA);
	int RepeatCount = (int) AA.Repeats.size();
	int AvgRepeatLength = GetAvgRepeatLength(AA);
	int AvgSpacerLength = GetAvgSpacerLength(AA);
	char Strand = (AA.ClusterRevComp ? '-' : '+');

	char *TmpLabel;
	int LocalLo = GlobalToLocal(AA.Pos, &TmpLabel);
	Out("%5d  %16.16s  %10d  %10d  %8d  %6d  %6d  %c  ",
	  AA.Id, Label, LocalLo + 1, ArrayLength, RepeatCount, AvgRepeatLength,
	  AvgSpacerLength, Strand);

	const Seq &s = AA.AlignedConsSeq;
	int Len = s.Length();

	for (int i = 0; i < Len; ++i)
		Out("%c", s.GetChar(i));
	Out("\n");
	}

static void OutputArraySummaryLinePos(const ArrayAln &AA, int Distance)
	{
	const ContigData &Contig = GlobalPosToContig(AA.Pos);
	const char *Label = Contig.Label;
	int ArrayLength = GetArrayLength(AA);
	int RepeatCount = (int) AA.Repeats.size();
	int AvgRepeatLength = GetAvgRepeatLength(AA);
	int AvgSpacerLength = GetAvgSpacerLength(AA);

	char *TmpLabel;
	int LocalLo = GlobalToLocal(AA.Pos, &TmpLabel);
	Out("%5d  %16.16s  %10d  %10d  %8d  %6d  %6d  ",
	  AA.Id, Label, LocalLo + 1, ArrayLength, RepeatCount, AvgRepeatLength, AvgSpacerLength);

	if (Distance == -1)
		Out("            ");
	else
		Out("%10d  ", Distance);

	const Seq &s = AA.ConsSeq;
	int Len = s.Length();

	for (int i = 0; i < Len; ++i)
		Out("%c", s.GetChar(i));
	Out("\n");
	}

static void OutputArrayFooter(const Seq &ConsSymbols)
	{
	for (int i = 0; i < (5+2+16+2+10+2+10+8+2+6+2+6+2+2+3); i++)
		Out(" ");
	int Len = ConsSymbols.Length();
	for (int i = 0; i < Len; ++i)
		Out("%c", ConsSymbols.GetChar(i));
	Out("\n");
	}

static bool ByPos(const ArrayAln *A1, const ArrayAln *A2)
	{
	return A1->Pos < A2->Pos;
	}

void OutArrays(std::vector<ArrayAln *> &AAVec)
	{
	IntVecVec Clusters;
	ClusterConsRC(AAVec, Clusters);
	const size_t ClusterCount = Clusters.size();

	Progress("Creating report");

// Detail report
	Out("\n");
	Out("\n");
	Out("DETAIL REPORT\n");
	Out("\n");
	const size_t ArrayCount = AAVec.size();
	for (size_t ArrayIndex = 0; ArrayIndex < ArrayCount; ++ArrayIndex)
		{
		const ArrayAln &AA = *(AAVec[ArrayIndex]);
		OutDetailAA(AA);
		}

// Summary reports
	Out("\n");
	Out("\n");
	Out("SUMMARY BY SIMILARITY\n");
	Out("\n");
	OutputArrayHeaderSim();

	for (size_t ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		const IntVec &Cluster = Clusters[ClusterIndex];
		Seq ConsSymbols;
		AlignCluster(AAVec, Cluster, ConsSymbols);
		const size_t ArrayCount = Cluster.size();
		for (size_t i = 0; i < ArrayCount; ++i)
			{
			unsigned ArrayIndex = Cluster[i];
			ArrayAln &AA = *(AAVec[ArrayIndex]);
			if (ArrayCount == 1)
				AA.ClusterRevComp = false;
			OutputArraySummaryLineSim(AA);
			}
		if (ArrayCount > 1)
			OutputArrayFooter(ConsSymbols);
		Out("\n");
		}

	Out("\n");
	Out("\n");
	Out("SUMMARY BY POSITION\n");
	Out("\n");

	std::sort(AAVec.begin(), AAVec.end(), ByPos);

	int LastContigIndex = -1;
	int LastPos = -1;
	for (size_t ArrayIndex = 0; ArrayIndex < ArrayCount; ++ArrayIndex)
		{
		const ArrayAln &AA = *(AAVec[ArrayIndex]);
		int Pos = AA.Pos;
		int ContigIndex = GlobalPosToContigIndex(Pos);
		if (ContigIndex != LastContigIndex)
			{
			LastPos = -1;
			const ContigData &Contig = GlobalPosToContig(Pos);
			Out("\n");
			Out("\n");
			Out(">");
			Out(Contig.Label);
			Out("\n");
			OutputArrayHeaderPos();
			LastContigIndex = Contig.Index;
			}

		int Distance = -1;
		if (ArrayIndex > 0)
			Distance = DistanceAA(*(AAVec[ArrayIndex-1]), AA);

		OutputArraySummaryLinePos(AA, Distance);
		}
	}
