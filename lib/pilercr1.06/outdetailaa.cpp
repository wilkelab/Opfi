#include "pilercr.h"

static void OutputDetailHeader(const ArrayAln &AA)
	{
	unsigned RepeatColCount = (unsigned) AA.Repeats[0].size();

	char *Label;
	GlobalToLocal(AA.Pos, &Label);
	Out("\n\n");
	Out("Array %d\n", AA.Id);
	
	Out(">%s\n", Label);
	Out("\n");

	Out("       Pos  Repeat     %%id  Spacer  Left flank    Repeat");
	int Blanks = (int) RepeatColCount - 6;
	for (int i = 0; i < Blanks; ++i)
		Out(" ");
	Out("    Spacer\n");
	Out("==========  ======  ======  ======  ==========    ");
	//   1234567890  123456  123456  123456  1234567890
	for (unsigned i = 0; i < RepeatColCount; ++i)
		Out("=");
	Out("    ======\n");
	}

int GetAvgSpacerLength(const ArrayAln &AA)
	{
	size_t RepeatCount = AA.Repeats.size();
	size_t Sum = 0;
	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount - 1; ++RepeatIndex)
		Sum += AA.Spacers[RepeatIndex].size();
	return (int) (Sum / (RepeatCount - 1));
	}

int GetAvgRepeatLength(const ArrayAln &AA)
	{
	size_t RepeatCount = AA.Repeats.size();
	size_t Sum = 0;
	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		Sum += AA.Repeats[RepeatIndex].size();
	return (int) (Sum / RepeatCount);
	}

static void OutputDetailFooter(const ArrayAln &AA)
	{
	unsigned RepeatColCount = (unsigned) AA.Repeats[0].size();

	Out("==========  ======  ======  ======  ==========    ");
	for (unsigned i = 0; i < RepeatColCount; ++i)
		Out("=");
	Out("\n");

	int RepeatCount = (int) AA.Repeats.size();
	int RepeatLength = (int) AA.Repeats[0].size();
	int SpacerLength = GetAvgSpacerLength(AA);

	Out("%10d  %6d          %6d                ",
	  RepeatCount, RepeatLength, SpacerLength);

	assert(AA.AlignedConsSeq.size() == RepeatColCount);
	for (unsigned i = 0; i < RepeatColCount; ++i)
		Out("%c", AA.AlignedConsSeq[i]);

	Out("\n");
	}

static double GetRepeatPctId(const ArrayAln &AA, size_t RepeatIndex)
	{
	assert(RepeatIndex < AA.Repeats.size());
	const std::string &Repeat = AA.Repeats[RepeatIndex];
	const Seq &ConsSeq = AA.AlignedConsSeq;
	size_t RepeatLength = Repeat.size();
	assert(ConsSeq.size() == RepeatLength);

	int Same = 0;
	for (size_t i = 0; i < RepeatLength; ++i)
		{
		if (toupper(ConsSeq[i]) == toupper(Repeat[i]))
			++Same;
		}

	return 100.0 * (double) Same / (double) RepeatLength;
	}

static void OutDetailRepeat(const ArrayAln &AA, size_t RepeatIndex, int &Pos)
	{
	const std::string &LeftFlank = AA.LeftFlanks[RepeatIndex];
	const std::string &Repeat = AA.Repeats[RepeatIndex];
	const std::string &Spacer = AA.Spacers[RepeatIndex];
	const Seq &ConsSeq = AA.AlignedConsSeq;

	size_t LeftFlankLength = LeftFlank.size();
	size_t RepeatLength = Repeat.size();
	size_t SpacerLength = Spacer.size();

	assert(LeftFlankLength == g_FlankSize);

	double PctId = GetRepeatPctId(AA, RepeatIndex);

	const ContigData &Contig = GlobalPosToContig(Pos);
	const char *Label = Contig.Label;
	char *TmpLabel;
	int LocalLo = GlobalToLocal(Pos, &TmpLabel);

	//          Pos  Repeat     %id  Spacer  Left flank    Repeat
	//   ==========  ======  ======  ======  ==========    ======
	//   1234567890  123456  123456  123456  1234567890
	Out("%10d  %6d  %6.1f  ",
	  LocalLo + 1,	// 1-based for user's benefit
	  RepeatLength,
	  PctId);

	const size_t RepeatCount = AA.Repeats.size();
	if (RepeatIndex == RepeatCount - 1)
		Out("        ");
	else
		Out("%6d  ", SpacerLength);

	for (size_t i = 0; i < LeftFlankLength; ++i)
		Out("%c", LeftFlank[i]);

	Out("    ");

	for (size_t i = 0; i < RepeatLength; ++i)
		{
		if (Repeat[i] == ConsSeq[i])
			{
			if (Repeat[i] == '-')
				Out("-");
			else
				Out(".");
			}
		else
			Out("%c", Repeat[i]);
		}

	Out("    ");

	for (size_t i = 0; i < SpacerLength; ++i)
		Out("%c", Spacer[i]);

	Out("\n");

	Pos += (int) (RepeatLength + SpacerLength);
	}

void OutDetailAA(const ArrayAln &AA)
	{
	OutputDetailHeader(AA);

	int Pos = AA.Pos;
	size_t RepeatCount = AA.Repeats.size();
	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		OutDetailRepeat(AA, RepeatIndex, Pos);

	OutputDetailFooter(AA);
	}
