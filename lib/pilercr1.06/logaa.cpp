#include "pilercr.h"

static void LogStdStr(const std::string &s)
	{
	for (size_t i = 0; i < s.size(); ++i)
		Log("%c", s[i]);
	}

void LogAA(const ArrayAln &AA)
	{
	const int RepeatCount = (int) AA.LeftFlanks.size();
	assert(RepeatCount == AA.Repeats.size());
	assert(RepeatCount == AA.Spacers.size());

	Log("AA Id=%d Pos=%d RepeatCount=%d\n", AA.Id, AA.Pos, RepeatCount);

	for (int RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		const std::string &LeftFlank = AA.LeftFlanks[RepeatIndex];
		const std::string &Repeat = AA.Repeats[RepeatIndex];
		const std::string &Spacer = AA.Spacers[RepeatIndex];

		LogStdStr(LeftFlank);
		Log("  ");
		LogStdStr(Repeat);
		Log("  ");
		LogStdStr(Spacer);
		Log("\n");
		}
	}
