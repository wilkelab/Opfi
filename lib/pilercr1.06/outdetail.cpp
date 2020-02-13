#include "pilercr.h"
#include "sw.h"

static void OutputDetailHeader(const ArrayData &AD, unsigned ColCount)
	{
	char *Label;
	GlobalToLocal(AD.Lo, &Label);
	Out("\n\n");
	Out("Array %d\n", AD.Id);
	
	if (g_ShowHits)
		{
		Out("Root: ");
		LogHit(AD.PileIndexes.front());
		Out("\n");
		}

	Out(">%s\n", Label);
	Out("\n");

	Out("       Pos  Repeat     %%id  Spacer  Left flank    Repeat");
	int Blanks = (int) ColCount - 6;
	for (int i = 0; i < Blanks; ++i)
		Out(" ");
	Out("    Spacer\n");
	Out("==========  ======  ======  ======  ==========    ");
	for (unsigned i = 0; i < ColCount; ++i)
		Out("=");
	Out("    ======\n");
	}

static void OutputDetailFooter(const ArrayData &AD, unsigned ColCount)
	{
	Out("==========  ======  ======  ======  ==========    ");
	for (unsigned i = 0; i < ColCount; ++i)
		Out("=");
	Out("\n");
	}

int PathALength(const char *Path)
	{
	int Length = 0;
	while (char c = *Path++)
		if (c == 'M' || c == 'D')
			++Length;
	return Length;
	}

double GetPctId(const char *A, const char *B, const char *Path)
	{
	int APos = 0;
	int BPos = 0;
	int Same = 0;
	int LA = 0;
	int LB = 0;
	while (char c = *Path++)
		{
		switch (c)
			{
		case 'M':
			if (toupper(A[APos++]) == toupper(B[BPos++]))
				++Same;
			++LA;
			++LB;
			break;

		case 'D':
			++APos;
			++LA;
			break;

		case 'I':
			++BPos;
			++LB;
			break;
			}
		}

	double MaxLen = std::max(LA, LB);
	if (MaxLen == 0)
		return 0;

	return (double) Same * 100.0 / MaxLen;
	}

void OutputArrayDetail(ArrayData &AD)
	{
	const IntVec &PileIndexes = AD.PileIndexes;

	std::vector<RowData> Rows;
	size_t RowCount = PileIndexes.size();
	Rows.reserve(RowCount);

	const unsigned ConsSeqLen = AD.ConsSeq.Length();
	char *ConsSeqStr = all(char, ConsSeqLen+1);
	AD.ConsSeq.ToString(ConsSeqStr, ConsSeqLen+1);

	SeqVect HitSeqs;
	int RowIndex = 0;
	for_CIntVec(PileIndexes, p)
		{
		RowData &Row = Rows[RowIndex];

		int PileIndex = *p;
		assert(PileIndex >= 0 && PileIndex < g_PileCount);
		Seq *PileSeq = new Seq;
		int Lo;
		int Hi;
		PileToSeq(g_Piles, PileIndex, *PileSeq, &Lo, &Hi);

		unsigned PileSeqLen = PileSeq->Length();
		char *PileSeqStr = all(char, PileSeqLen+1);
		PileSeq->ToString(PileSeqStr, PileSeqLen+1);

		std::string Path;
		unsigned PileOffset;
		unsigned ConsOffset;
		SWSimple(PileSeqStr, PileSeqLen, ConsSeqStr, ConsSeqLen, &PileOffset,
		  &ConsOffset, Path);

		Row.PctId = GetPctId(PileSeqStr+PileOffset, ConsSeqStr+ConsOffset, Path.c_str());

		const PileData &Pile = g_Piles[PileIndex];
		Row.RepeatLo = Lo + PileOffset;
		if (RowIndex == 0)
			AD.Lo = Row.RepeatLo;

		int RepeatLength = PathALength(Path.c_str());
		Row.RepeatLength = RepeatLength;

		Seq *HitSeq = new Seq;
		const int ToPos = PileOffset + RepeatLength - 1;
		for (int Pos = PileOffset; Pos <= ToPos; ++Pos)
			HitSeq->push_back(PileSeqStr[Pos]);
		HitSeq->SetName("hit");
		HitSeq->SetId(0);

		HitSeqs.push_back(HitSeq);

		++RowIndex;
		}

	freemem(ConsSeqStr);
	ConsSeqStr = 0;

	Seq *ConsSeq = new Seq;
	ConsSeq->Copy(AD.ConsSeq);
	HitSeqs.push_back(ConsSeq);
	MSA Aln;
	MultipleAlign(HitSeqs, Aln);

	const unsigned ColCount = Aln.GetColCount();

	OutputDetailHeader(AD, ColCount);

	int SumSpacerLengths = 0;
	for (size_t RowIndex = 0; RowIndex < RowCount; ++RowIndex)
		{
		const RowData &Row = Rows[RowIndex];
		char *Label;
		int LocalRepeatLo = GlobalToLocal(Row.RepeatLo, &Label);
		Out("%10d  %6d  %6.1f%%", LocalRepeatLo + 1, Row.RepeatLength, Row.PctId);

		int NextLo = -1;
		int Hi = -1;
		int SpacerLength = -1;

	// Spacer
		if (RowIndex + 1 == RowCount)
			Out("         ");
		else
			{
			NextLo = Rows[RowIndex+1].RepeatLo;
			Hi = Row.RepeatLo + Row.RepeatLength - 1;
			SpacerLength = NextLo - Hi;
			Out("  %5d  ", SpacerLength);
			SumSpacerLengths += SpacerLength;
			}

	// Left flank
		const ContigData &Contig = GlobalPosToContig(Row.RepeatLo);
		int ContigFrom = Contig.From;
		int ContigTo = ContigFrom + Contig.Length - 1;
		int LeftFlankStartPos = Row.RepeatLo - g_FlankSize;
		int LeftFlankBlanks = 0;
		if (LeftFlankStartPos < ContigFrom)
			{
			LeftFlankStartPos = ContigFrom;
			LeftFlankBlanks = g_FlankSize - (Row.RepeatLo - ContigFrom);
			}
		for (int i = 0; i < LeftFlankBlanks; ++i)
			Out(" ");
		for (int Pos = LeftFlankStartPos; Pos < Row.RepeatLo; ++Pos)
			Out("%c", g_SeqQ[Pos]);

		Out("    ");

	// Repeat
		for (unsigned Col = 0; Col < ColCount; ++Col)
			Out("%c", Aln.GetChar((unsigned) RowIndex, Col));

		Out("    ");

	// Spacer
		if (RowIndex + 1 == RowCount)
			{
			int RightFlankPos = Row.RepeatLo + Row.RepeatLength;
			int ToPos = RightFlankPos + 9;
			if (ToPos >= ContigTo)
				ToPos = ContigTo;
			for (int Pos = RightFlankPos; Pos <= ToPos; ++Pos)
				Out("%c", g_SeqQ[Pos]);
			}
		else
			{
			for (int Pos = Hi + 1; Pos < NextLo; ++Pos)
				Out("%c", g_SeqQ[Pos]);
			}

		if (g_ShowHits)
			{
			Out("  ");
			int PileIndex = AD.PileIndexes[RowIndex];
			LogPile(PileIndex);
			}

		Out("\n");
		}

	AD.SpacerLength = (int) ((double) SumSpacerLengths / (double) (RowCount - 1) + 0.5);
	OutputDetailFooter(AD, ColCount);

	Out("%10d  %6d          %6d                ",
	  (int) AD.PileIndexes.size(), AD.RepeatLength, AD.SpacerLength);

	for (unsigned Col = 0; Col < ColCount; ++Col)
		Out("%c", Aln.GetChar((unsigned) RowCount, Col));
	Out("\n");

	Seq ConsSymbols;
	GetConsSymbols(Aln, ConsSymbols);

	for (int i = 0; i < (10+2+6+4+28); ++i)
		Out(" ");

	assert(ConsSymbols.Length() == ColCount);
	for (unsigned Col = 0; Col < ColCount; ++Col)
		Out("%c", ConsSymbols.GetChar(Col));
	Out("\n");
	}
