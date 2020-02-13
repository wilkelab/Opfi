#include "pilercr.h"

#define	VALIDATE_ALL		1

#if	VALIDATE_ALL
#define ValidateAA_All	ValidateAA
#else
#define ValidateAA_All	/* empty */
#endif

static char LeftCol(ArrayAln &AA, size_t RepeatIndex)
	{
	std::string &LeftFlank = AA.LeftFlanks[RepeatIndex];
	size_t Length = LeftFlank.size();
	if (Length == 0)
		return 0;
	return LeftFlank[Length-1];
	}

static bool LeftColConserved(ArrayAln &AA)
	{
	const size_t RepeatCount = AA.LeftFlanks.size();
	assert(RepeatCount == AA.Repeats.size());
	assert(RepeatCount == AA.Spacers.size());

	char c0 = LeftCol(AA, 0);
	if (c0 == 0 || c0 == '-')
		return false;

	for (size_t RepeatIndex = 1; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		char c = LeftCol(AA, RepeatIndex);
		if (c != c0)
			return false;
		}

	return true;
	}

static char RightCol(ArrayAln &AA, size_t RepeatIndex)
	{
	std::string &Spacer = AA.Spacers[RepeatIndex];
	size_t Length = Spacer.size();
	if (Length == 0)
		return 0;
	return Spacer[0];
	}

static bool RightColConserved(ArrayAln &AA)
	{
	const size_t RepeatCount = AA.Spacers.size();
	assert(RepeatCount == AA.Repeats.size());
	assert(RepeatCount == AA.Spacers.size());

	char c0 = RightCol(AA, 0);
	if (c0 == 0)
		return false;

	for (size_t RepeatIndex = 1; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		char c = RightCol(AA, RepeatIndex);
		if (c != c0)
			return false;
		}

	return true;
	}

static bool FirstColConserved(const ArrayAln &AA)
	{
	const int *Counts = GetCountsAA(AA, 0);

	int MaxCount = 0;
	for (int i = 0; i < 4; ++i)
		if (Counts[i] > MaxCount)
			MaxCount = Counts[i];

	const size_t RepeatCount = AA.LeftFlanks.size();
	double Cons = (double) MaxCount / (double) RepeatCount;
	return RepeatCount - MaxCount <= 1 || Cons >= g_DraftMinColCons;
	}

static bool LastColConserved(const ArrayAln &AA)
	{
	int RepeatLength = (int) AA.Repeats[0].size();

	const int *Counts = GetCountsAA(AA, RepeatLength - 1);

	int MaxCount = 0;
	for (int i = 0; i < 4; ++i)
		if (Counts[i] > MaxCount)
			MaxCount = Counts[i];

	const size_t RepeatCount = AA.LeftFlanks.size();
	double Cons = (double) MaxCount / (double) RepeatCount;
	return RepeatCount - MaxCount <= 1 || Cons >= g_DraftMinColCons;
	}

static void MoveLeftColToRepeat(ArrayAln &AA)
	{
	--(AA.Pos);

	const size_t RepeatCount = AA.LeftFlanks.size();
	assert(RepeatCount == AA.Repeats.size());
	assert(RepeatCount == AA.Spacers.size());
	if (RepeatCount == 0)
		return;

	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		std::string &LeftFlank = AA.LeftFlanks[RepeatIndex];
		std::string &Repeat = AA.Repeats[RepeatIndex];
		std::string &Spacer = AA.Spacers[RepeatIndex];

		size_t LeftFlankLength = LeftFlank.size();
		size_t SpacerLength = Spacer.size();

		assert(LeftFlankLength > 0);
		char c = LeftFlank[LeftFlankLength-1];
		if (c == ' ')
			c = '-';

		Repeat = c + Repeat;
		if (RepeatIndex != RepeatCount - 1)
			Spacer = Spacer.substr(0, SpacerLength-1);

		if (RepeatIndex > 0)
			{
			std::string &PrevSpacer = AA.Spacers[RepeatIndex-1];
			size_t PrevSpacerLength = PrevSpacer.size();
			if (PrevSpacerLength >= (size_t) g_FlankSize)
				{
				size_t StartPos = PrevSpacerLength - g_FlankSize;
				LeftFlank = PrevSpacer.substr(StartPos, g_FlankSize);
				}
			else
				{
				LeftFlank.clear();
				size_t BlankCount = g_FlankSize - PrevSpacerLength;
				for (size_t i = 0; i < BlankCount; ++i)
					LeftFlank += ' ';
				LeftFlank += PrevSpacer;
				}
			}
		else
			{
			char *Label;
			int LocalPos = GlobalToLocal(AA.Pos, &Label);
			char c = ' ';
			if (LocalPos >= g_FlankSize)
				c = g_SeqQ[AA.Pos - g_FlankSize];
			LeftFlank = c + LeftFlank.substr(0, g_FlankSize-1);
			}
		}

	Log("MoveLeftColToRepeat\n");
	}

// Move left-most column in repeat matrix to left flank
static void MoveRepeatToLeftCol(ArrayAln &AA)
	{
	if (AA.Repeats[0][0] != '-')
		++(AA.Pos);

	const size_t RepeatCount = AA.LeftFlanks.size();
	assert(RepeatCount == AA.Repeats.size());
	assert(RepeatCount == AA.Spacers.size());
	if (RepeatCount == 0)
		return;

	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		std::string &LeftFlank = AA.LeftFlanks[RepeatIndex];
		std::string &Repeat = AA.Repeats[RepeatIndex];

		size_t RepeatLength = Repeat.size();
		size_t LeftFlankLength = LeftFlank.size();

		char c = Repeat[0];

		Repeat = Repeat.substr(1, RepeatLength - 1);
		if (c != '-')
			{
			LeftFlank = LeftFlank.substr(1, LeftFlankLength - 1);
			LeftFlank += c;
			}

		if (RepeatIndex > 0 && c != '-')
			{
			std::string &PrevSpacer = AA.Spacers[RepeatIndex - 1];
			PrevSpacer += c;
			}

		assert(Repeat.size() == RepeatLength - 1);
		assert(LeftFlank.size() == LeftFlankLength);
		}

	Log("MoveRepeatToLeftCol\n");
	}

// Move right-most column in repeat matrix to spacer
static void MoveRepeatToRightCol(ArrayAln &AA)
	{
	const size_t RepeatCount = AA.LeftFlanks.size();
	assert(RepeatCount == AA.Repeats.size());
	assert(RepeatCount == AA.Spacers.size());
	if (RepeatCount == 0)
		return;

	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		std::string &Repeat = AA.Repeats[RepeatIndex];
		std::string &Spacer = AA.Spacers[RepeatIndex];

		size_t RepeatLength = Repeat.size();

		char c = Repeat[RepeatLength - 1];

		Repeat = Repeat.substr(0, RepeatLength - 1);
		std::string NewSpacer;
		if (c != '-')
			NewSpacer = c;
		NewSpacer += Spacer;
		Spacer = NewSpacer;
		}

	Log("MoveRepeatToRightCol\n");
	}

static void MoveRightColToRepeat(ArrayAln &AA)
	{
	const size_t RepeatCount = AA.LeftFlanks.size();
	assert(RepeatCount == AA.Repeats.size());
	assert(RepeatCount == AA.Spacers.size());
	if (RepeatCount == 0)
		return;

	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		std::string &Spacer = AA.Spacers[RepeatIndex];
		std::string &Repeat = AA.Repeats[RepeatIndex];

		size_t SpacerLength = Spacer.size();
		size_t RepeatLength = Repeat.size();
		assert(SpacerLength > 0);

		char c = Spacer[0];
		std::string NewRepeat;
		NewRepeat += AA.Repeats[RepeatIndex];
		NewRepeat += c;
		Repeat = NewRepeat;

		Spacer = Spacer.substr(1, SpacerLength-1);

		assert(Spacer.size() == SpacerLength - 1);
		assert(Repeat.size() == RepeatLength + 1);
		}

	Log("MoveRightColToRepeat\n");
	}

static void FixLeftTerminalGap(ArrayAln &AA, size_t RepeatIndex)
	{
	std::string &Repeat = AA.Repeats[RepeatIndex];
	const size_t RepeatLength = Repeat.size();
	if (RepeatLength == 0)
		return;

	size_t GapLength = 0;
	for (size_t i = 0; i < RepeatLength && Repeat[i] == '-'; ++i)
		++GapLength;

	if (GapLength == 0)
		return;

	if (GapLength > 3)
		return;

	if (RepeatIndex == 0)
		return;

	std::string &PrevSpacer = AA.Spacers[RepeatIndex - 1];
	size_t PrevSpacerLength = PrevSpacer.size();
	if (PrevSpacerLength < (size_t) g_FlankSize)
		return;
	if (PrevSpacerLength < GapLength)
		return;
	if (PrevSpacerLength - GapLength < (size_t) g_FlankSize)
		return;

	const std::string &LeftFlank = AA.LeftFlanks[RepeatIndex];
	size_t LeftFlankLength = LeftFlank.size();
	if (LeftFlankLength <= GapLength)
		return;

	std::string NewRepeat;
	NewRepeat += LeftFlank.substr(LeftFlankLength - GapLength, GapLength);
	NewRepeat += Repeat.substr(GapLength, RepeatLength - GapLength);
	AA.Repeats[RepeatIndex] = NewRepeat;

	PrevSpacer = PrevSpacer.substr(0, PrevSpacerLength - GapLength);
	PrevSpacerLength = PrevSpacer.size();

	assert(PrevSpacerLength >= (size_t) g_FlankSize);

	std::string TmpStr = PrevSpacer.substr(PrevSpacerLength - g_FlankSize, g_FlankSize);
	AA.LeftFlanks[RepeatIndex] = TmpStr;

	Log("FixLeftTerminalGap %d\n", (int) RepeatIndex);
	}

static void FixRightTerminalGap(ArrayAln &AA, size_t RepeatIndex)
	{
	std::string &Repeat = AA.Repeats[RepeatIndex];
	char *Label;
	int LocalPos = GlobalToLocal(AA.Pos, &Label);
	const size_t RepeatLength = Repeat.size();
	if (RepeatLength == 0)
		return;

	size_t GapLength = 0;
	for (int i = (int) RepeatLength - 1; i >= 0 && Repeat[i] == '-'; --i)
		{
		++GapLength;
		}

	if (GapLength == 0)
		return;

	if (GapLength > 2)
		return;

	std::string &Spacer = AA.Spacers[RepeatIndex];
	size_t SpacerLength = Spacer.size();

	const size_t RepeatCount = AA.Repeats.size();
	if (SpacerLength < GapLength)
		{
		if (RepeatIndex != RepeatCount - 1)
			Quit("SpacerLength < GapLength");
		return;
		}

	Repeat = Repeat.substr(0, RepeatLength - GapLength);
	Repeat += Spacer.substr(0, GapLength);
	if (SpacerLength >= GapLength)
		Spacer = Spacer.substr(GapLength, SpacerLength - GapLength);
	else
		Spacer.clear();

	Log("FixRightTerminalGap %d\n", (int) RepeatIndex);
	}

static int GapCount(const std::string &s)
	{
	int Sum = 0;
	for (size_t i = 0; i < s.size(); ++i)
		if (s[i] == '-')
			++Sum;
	return Sum;
	}

static int AvgGapCount(const ArrayAln &AA)
	{
	const size_t RepeatCount = AA.Repeats.size();
	int Sum = 0;
	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		const std::string &Repeat = AA.Repeats[RepeatIndex];
		Sum += GapCount(Repeat);
		}
	return (int) (((double) Sum / (double) RepeatCount) + 0.5);
	}

static bool HasInitialPartialRepeat(const ArrayAln &AA)
	{
	int Avg = AvgGapCount(AA);
	int Count = GapCount(AA.Repeats[0]);
	return Count >= 3 && Count > Avg*2;
	}

static bool HasFinalPartialRepeat(const ArrayAln &AA)
	{
	int Avg = AvgGapCount(AA);
	size_t RepeatCount = AA.Repeats.size();
	int Count = GapCount(AA.Repeats[RepeatCount-1]);
	return Count >= 3 && Count > Avg*2;
	}

static void DeleteInitialPartialRepeat(ArrayAln &AA)
	{
	const std::string &FirstRepeat = AA.Repeats[0];
	const std::string &FirstSpacer = AA.Spacers[0];

	int Offset = 0;
	for (size_t i = 0; i < FirstRepeat.size(); ++i)
		if (FirstRepeat[i] != '-')
			++Offset;
	Offset += (int) FirstSpacer.size();
	AA.Pos += Offset;

	const size_t RepeatCount = AA.Repeats.size();
	for (int i = 0; i < (int) RepeatCount - 1; ++i)
		{
		AA.LeftFlanks[i] = AA.LeftFlanks[i+1];
		AA.Repeats[i] = AA.Repeats[i+1];
		AA.Spacers[i] = AA.Spacers[i+1];
		}

	AA.LeftFlanks.resize(RepeatCount-1);
	AA.Repeats.resize(RepeatCount-1);
	AA.Spacers.resize(RepeatCount-1);

	Log("DeleteInitialPartialRepeat\n");
	}

static void DeleteFinalPartialRepeat(ArrayAln &AA)
	{
	const std::string &LastRepeat = AA.Repeats[0];
	const std::string &LastSpacer = AA.Spacers[0];

	size_t RepeatCount = AA.Repeats.size();
	AA.LeftFlanks.resize(RepeatCount-1);
	AA.Repeats.resize(RepeatCount-1);
	AA.Spacers.resize(RepeatCount-1);

	Log("DeleteFinalPartialRepeat\n");
	}

void FixBounds(ArrayAln &AA)
	{
	ValidateAA(AA);

	if (HasInitialPartialRepeat(AA))
		{
		if (AA.Repeats.size() <= (size_t) g_DraftMinArraySize)
			return;
		DeleteInitialPartialRepeat(AA);
		}

	if (HasFinalPartialRepeat(AA))
		{
		if (AA.Repeats.size() <= (size_t) g_DraftMinArraySize)
			return;
		DeleteFinalPartialRepeat(AA);
		}

	size_t RepeatCount = AA.Repeats.size();
	assert(RepeatCount == AA.LeftFlanks.size());
	assert(RepeatCount == AA.Spacers.size());

	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		ValidateAA_All(AA);
		FixLeftTerminalGap(AA, RepeatIndex);
		ValidateAA_All(AA);
		FixRightTerminalGap(AA, RepeatIndex);
		ValidateAA_All(AA);
		}

	while (LeftColConserved(AA))
		{
		ValidateAA_All(AA);
		MoveLeftColToRepeat(AA);
		ValidateAA_All(AA);
		}
	while (RightColConserved(AA))
		{
		ValidateAA_All(AA);
		MoveRightColToRepeat(AA);
		ValidateAA_All(AA);
		}

	while (!FirstColConserved(AA))
		{
		ValidateAA_All(AA);
		MoveRepeatToLeftCol(AA);
		ValidateAA_All(AA);
		}
	while (!LastColConserved(AA))
		{
		ValidateAA_All(AA);
		MoveRepeatToRightCol(AA);
		ValidateAA_All(AA);
		}

	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		ValidateAA_All(AA);
		FixLeftTerminalGap(AA, RepeatIndex);
		ValidateAA_All(AA);
		FixRightTerminalGap(AA, RepeatIndex);
		ValidateAA_All(AA);
		}

	ValidateAA(AA);
	}
