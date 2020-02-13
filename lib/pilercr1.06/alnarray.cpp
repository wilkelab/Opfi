#include "pilercr.h"
#include "sw.h"

static void AllocAA(ArrayAln &AA, int RepeatCount)
	{
	for (int i = 0; i < RepeatCount; ++i)
		{
		AA.LeftFlanks.push_back(*new std::string);
		AA.Repeats.push_back(*new std::string);
		AA.Spacers.push_back(*new std::string);
		}
	}

static ArrayAln *CopyAA(const ArrayAln &AA)
	{
	ArrayAln *NewAA = new ArrayAln;
	*NewAA = AA;
	return NewAA;
	}

static bool EqAA(const ArrayAln &AA1, const ArrayAln &AA2)
	{
	if (AA1.Id != AA2.Id)
		return false;

	if (AA1.Pos != AA2.Pos)
		return false;

	const size_t RepeatCount = AA1.Repeats.size();
	if (RepeatCount != AA2.Repeats.size())
		return false;

	for (unsigned RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		if (AA1.Repeats != AA2.Repeats)
			return false;
		if (AA1.LeftFlanks != AA2.LeftFlanks)
			return false;
		if (AA1.Spacers != AA2.Spacers)
			return false;
		}
	return true;
	}

static size_t UngappedLength(const std::string &s)
	{
	size_t Length =  0;
	for (size_t i = 0; i < s.size(); ++i)
		if (s[i] != '-')
			++Length;
	return Length;
	}

int GetArrayLength(const ArrayAln &AA)
	{
	size_t Length = 0;
	const size_t RepeatCount = AA.Repeats.size();
	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		size_t RepeatLength = UngappedLength(AA.Repeats[RepeatIndex]);
		Length += RepeatLength;
		if (RepeatIndex != RepeatCount - 1)
			{
			size_t SpacerLength = AA.Spacers[RepeatIndex].size();
			Length += SpacerLength;
			}
		}
	return (int) Length;
	}

void AlignArray(const ArrayData &AD, ArrayAln &AA)
	{
	const IntVec &PileIndexes = AD.PileIndexes;

	std::vector<RowData> Rows;
	size_t RowCount = PileIndexes.size();
	Rows.resize(RowCount);

	AllocAA(AA, (int) RowCount);

	AA.Id = AD.Id;

	const unsigned ConsSeqLen = AD.ConsSeq.Length();
	char *ConsSeqStr = all(char, ConsSeqLen+1);
	AD.ConsSeq.ToString(ConsSeqStr, ConsSeqLen+1);

	SeqVect HitSeqs;
	int RowIndex = 0;
	for_CIntVec(PileIndexes, p)
		{
		RowData &Row = Rows[RowIndex++];

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
		int RepeatLength = PathALength(Path.c_str());
		Row.RepeatLength = RepeatLength;

		Seq *HitSeq = new Seq;
		const int ToPos = PileOffset + RepeatLength - 1;
		for (int Pos = PileOffset; Pos <= ToPos; ++Pos)
			HitSeq->push_back(PileSeqStr[Pos]);
		HitSeq->SetName("hit");
		HitSeq->SetId(0);

		HitSeqs.push_back(HitSeq);
		}

	freemem(ConsSeqStr);
	ConsSeqStr = 0;

	Seq *ConsSeq = new Seq;
	ConsSeq->Copy(AD.ConsSeq);
	HitSeqs.push_back(ConsSeq);
	MSA Aln;
	MultipleAlign(HitSeqs, Aln);

	if (g_LogAlns)
		{
		Out("Refined alignment array %d:\n", AD.Id);
		Aln.LogMe();
		Out("\n");
		}

	const unsigned ColCount = Aln.GetColCount();

	int SumSpacerLengths = 0;
	for (size_t RowIndex = 0; RowIndex < RowCount; ++RowIndex)
		{
		const RowData &Row = Rows[RowIndex];
		char *Label;
		int LocalRepeatLo = GlobalToLocal(Row.RepeatLo, &Label);

		int NextLo = -1;
		int Hi = -1;
		int SpacerLength = -1;
		if (RowIndex + 1 != RowCount)
			{
			NextLo = Rows[RowIndex+1].RepeatLo;
			Hi = Row.RepeatLo + Row.RepeatLength - 1;
			SpacerLength = NextLo - Hi;
			if (SpacerLength < 4)
				{
				AA.Repeats.clear();
				AA.Spacers.clear();
				AA.LeftFlanks.clear();
				return;
				}
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

		std::string &LeftFlank = AA.LeftFlanks[RowIndex];
		for (int i = 0; i < LeftFlankBlanks; ++i)
			LeftFlank += ' ';
		for (int Pos = LeftFlankStartPos; Pos < Row.RepeatLo; ++Pos)
			LeftFlank += g_SeqQ[Pos];

	// Repeat
		std::string &Repeat = AA.Repeats[RowIndex];
		for (unsigned Col = 0; Col < ColCount; ++Col)
			Repeat += Aln.GetChar((unsigned) RowIndex, Col);

	// Spacer
		std::string &Spacer = AA.Spacers[RowIndex];
		if (RowIndex + 1 == RowCount)
			{
			int RightFlankPos = Row.RepeatLo + Row.RepeatLength;
			int ToPos = RightFlankPos + 9;
			if (ToPos >= ContigTo)
				ToPos = ContigTo;
			for (int Pos = RightFlankPos; Pos <= ToPos; ++Pos)
				Spacer += g_SeqQ[Pos];
			}
		else
			{
			for (int Pos = Hi + 1; Pos < NextLo; ++Pos)
				Spacer += g_SeqQ[Pos];
			}
		}

	AA.Pos = Rows[0].RepeatLo;
	FixBounds(AA);

	GetConsSeqAA(AA, AA.AlignedConsSeq);
	AA.ConsSeq.Copy(AA.AlignedConsSeq);
	StripBlanksAndGaps(AA.ConsSeq);
	}
