#include "pilercr.h"
#include "sw.h"

int UnalignedLength(const std::string &s)
	{
	size_t L = s.size();
	int N = 0;
	for (size_t i = 0; i < L; ++i)
		if (s[i] != '-')
			++N;
	return N;
	}

static bool GloballyAlignable(const char *Seq1, unsigned L1, const char *Seq2,
  unsigned L2)
	{
	std::string Path;
	unsigned Start1;
	unsigned Start2;
	score_t Score = SWSimple(Seq1, L1, Seq2, L2, &Start1, &Start2, Path);

	if (Score <= 0)
		return false;

	int AlnLen = (int) Path.size();
	if (GetRatio(AlnLen, (int) L1) < 0.95)
		return false;

	return true;
	}

// is i = nj for n approximately integer in range 1 .. 4?
static int Multiple(int i, int j)
	{
	if (GetRatio(i, 1*j) >= 0.95)
		return 1;

	if (GetRatio(i, 2*j) >= 0.95)
		return 2;

	if (GetRatio(i, 3*j) >= 0.95)
		return 3;

	if (GetRatio(i, 4*j) >= 0.95)
		return 4;

	return -1;
	}

int DistanceAA(const ArrayAln &AA1, const ArrayAln &AA2)
	{
	int End1 = AA1.Pos + GetArrayLength(AA1) + GetAvgSpacerLength(AA1) - 1;
	int Pos2 = AA2.Pos;

	int ContigIndex1 = GlobalPosToContigIndex(End1);
	int ContigIndex2 = GlobalPosToContigIndex(Pos2);

	if (ContigIndex1 != ContigIndex2)
		return -1;

	int ci1 = GlobalPosToContigIndex(AA1.Pos);
	if (ContigIndex1 != ci1)
		Warning("ContigIndex1=%u GlobalPosToContigIndex(AA1.Pos=%d) = %d",
		  ContigIndex1, AA1.Pos, ci1);

	return Pos2 - End1 - 1;
	}

bool MergeAAs(ArrayAln &AA1, const ArrayAln &AA2)
	{
	int RepeatLength1 = GetAvgRepeatLength(AA1);
	int RepeatLength2 = GetAvgRepeatLength(AA2);

	if (RepeatLength1 != RepeatLength2)
		return false;

	int SpacerLength1 = GetAvgSpacerLength(AA1);
	int SpacerLength2 = GetAvgSpacerLength(AA2);

	if (GetRatio(SpacerLength1, SpacerLength2) < 0.95)
		return false;

// Is distance between arrays a multiple of repeat+spacer length?
	int Distance = DistanceAA(AA1, AA2);
	int RepeatPlusSpacerLength = RepeatLength1 + SpacerLength1;
	int MissedRepeatCount = Multiple(Distance, RepeatPlusSpacerLength);
	if (MissedRepeatCount < 0)
		return false;

	int ConsSeqLength1 = AA1.ConsSeq.Length();
	int ConsSeqLength2 = AA2.ConsSeq.Length();

	char *Seq1 = all(char, ConsSeqLength1 + 1);
	char *Seq2 = all(char, ConsSeqLength2 + 1);

	AA1.ConsSeq.ToString(Seq1, ConsSeqLength1 + 1);
	AA2.ConsSeq.ToString(Seq2, ConsSeqLength2 + 1);

	if (!GloballyAlignable(Seq1, ConsSeqLength1, Seq2, ConsSeqLength2))
		return false;

	IntVec RepeatStartPos;
	IntVec RepeatLengths;

	int LastSpacerStartPos = AA1.Pos + GetArrayLength(AA1);
	int EndPosOfLastRepeat = LastSpacerStartPos - 1;
	StrVec Repeats;
	for (int i = 0; i < MissedRepeatCount; ++i)
		{
		int GuessedRepeatPos = EndPosOfLastRepeat + SpacerLength1;
		int StartPos = GuessedRepeatPos - 8;
		int EndPos = GuessedRepeatPos + RepeatLength1 + 8;
		int SeqLength = EndPos - StartPos + 1;
		char *Seq = all(char, SeqLength);
		for (int i = 0; i < SeqLength; ++i)
			Seq[i] = g_SeqQ[StartPos + i];

		char *ConsSeq = all(char, ConsSeqLength1 + 1);
		AA1.ConsSeq.ToString(ConsSeq, ConsSeqLength1 + 1);

		std::string Path;
		unsigned StartSeq;
		unsigned StartCons;
		score_t Score = SWSimple(Seq, SeqLength, ConsSeq, ConsSeqLength1, 
		  &StartSeq, &StartCons, Path);

		int AlnLength = (int) Path.size();
		if (GetRatio(AlnLength, ConsSeqLength1) < 0.95)
			return false;

		std::string &Repeat = *new std::string;
		for (int i = 0; i < (int) StartCons; ++i)
			Repeat.push_back('-');

		int SeqPos = StartPos + StartSeq;
		int RepeatPos = SeqPos;
		RepeatStartPos.push_back(RepeatPos);
		int RepeatLength = 0;
		for (int i = 0; i < AlnLength; ++i)
			{
			if (Path[i] == 'M' || Path[i] == 'D')
				{
				char c = g_SeqQ[SeqPos];
				Repeat.push_back(c);
				++SeqPos;
				++RepeatLength;
				}
			else
				Repeat.push_back('-');
			}
		if (!RepeatLengths.empty() && RepeatLength != RepeatLengths[0])
			return false;
		RepeatLengths.push_back(RepeatLength);
		EndPosOfLastRepeat = RepeatPos + RepeatLength - 1;

		while (Repeat.size() < (size_t) ConsSeqLength1)
			Repeat.push_back('-');

		Repeats.push_back(Repeat);
		}

	for (int i = 0; i < MissedRepeatCount; ++i)
		AA1.Repeats.push_back(Repeats[i]);

	for (int i = 0; i < MissedRepeatCount; ++i)
		{
		std::string &LeftFlank = *new std::string;
		int LeftFlankStartPos = RepeatStartPos[i] - g_FlankSize;
		for (int i = 0; i < g_FlankSize; ++i)
			{
			char c = g_SeqQ[LeftFlankStartPos + i];
			LeftFlank.push_back(c);
			}
		AA1.LeftFlanks.push_back(LeftFlank);
		}

	int OldCount = (int) AA1.Spacers.size();

	std::string &Spacer = *new std::string;
	int RepeatStartPos0 = RepeatStartPos[0];
	int LastSpacerEndPos = RepeatStartPos0 - 1;
	assert(LastSpacerEndPos > LastSpacerStartPos);
	int Length = LastSpacerEndPos - LastSpacerStartPos + 1;
	for (int i = 0; i < Length; ++i)
		{
		char c = g_SeqQ[LastSpacerStartPos + i];
		Spacer.push_back(c);
		}
	AA1.Spacers[OldCount-1] = Spacer;

	for (int i = 0; i < MissedRepeatCount; ++i)
		{
		std::string &Spacer = *new std::string;
		int SpacerStartPos = RepeatStartPos[i] + UnalignedLength(AA1.Repeats[OldCount+i]);
		int SpacerEndPos;
		if (i == MissedRepeatCount - 1)
			SpacerEndPos = AA2.Pos - 1;
		else
			SpacerEndPos = RepeatStartPos[i+1] - 1;
		assert(SpacerEndPos > SpacerStartPos);
		int Length = SpacerEndPos - SpacerStartPos + 1;
		for (int i = 0; i < Length; ++i)
			{
			char c = g_SeqQ[SpacerStartPos + i];
			Spacer.push_back(c);
			}
		AA1.Spacers.push_back(Spacer);
		}

	size_t RepeatCount2 = AA2.Repeats.size();
	size_t RepeatLength = 0;
	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount2; ++RepeatIndex)
		{
		AA1.LeftFlanks.push_back(AA2.LeftFlanks[RepeatIndex]);
		std::string Repeat = AA2.Repeats[RepeatIndex];
		AA1.Repeats.push_back(Repeat);
		const std::string &Spacer = AA2.Spacers[RepeatIndex];
		size_t L = Spacer.size();
		AA1.Spacers.push_back(Spacer);
		}

	ValidateAA(AA1);
	Log("MergeAAs\n");
	return true;
	}
