#include "pilercr.h"

static double GetCons(const std::string &s1, const Seq &s2)
	{
	assert(s1.size() == s2.size());
	size_t L = s1.size();
	int Sum = 0;
	for (size_t i = 0; i < L; ++i)
		{
		if (toupper(s1[i]) == toupper(s2[i]))
			++Sum;
		}
	return (double) Sum / (double) L;
	}

bool AcceptedAA(const ArrayAln &AA)
	{
	const size_t ArraySize = AA.Repeats.size();
	if (ArraySize < (size_t) g_MinArraySize)
		{
		Log("Rejected size=%d, minarray=%d\n", (int) ArraySize, g_MinArraySize);
		LogAA(AA);
		return false;
		}

	int MinRepeatLength = (int) AA.Repeats[0].size();
	int MaxRepeatLength = 0;

	int AvgSpacerLength = GetAvgSpacerLength(AA);
	int AvgRepeatLength = GetAvgRepeatLength(AA);

	if (AvgSpacerLength < g_MinSpacerLength)
		{
		Log("Rejected, avg spacer %d < min %d\n",
		  AvgSpacerLength, g_MinSpacerLength);
		return false;
		}

	if (AvgSpacerLength > g_MaxSpacerLength)
		{
		Log("Rejected, avg spacer %d < max %d\n",
		  AvgSpacerLength, g_MaxSpacerLength);
		return false;
		}

	if (AvgRepeatLength < g_MinRepeatLength)
		{
		Log("Rejected min repeat length array min %d min allowed %d\n",
		  AvgRepeatLength, g_MinRepeatLength);
		LogAA(AA);
		return false;
		}

	if (AvgRepeatLength > g_MaxRepeatLength)
		{
		Log("Rejected Max repeat length array Max %d Max allowed %d\n",
		  AvgRepeatLength, g_MaxRepeatLength);
		LogAA(AA);
		return false;
		}

	int RunLength = 0;
	int MaxRunLength = 0;
	int RunStart = 0;
	int LastRunStart = 0;
	int RunEnd = 0;
	for (size_t RepeatIndex = 0; RepeatIndex < ArraySize; ++RepeatIndex)
		{
		const std::string &Repeat = AA.Repeats[RepeatIndex];

		double Cons = GetCons(Repeat, AA.AlignedConsSeq);
		if (Cons >= g_MinCons)
			{
			if (RunLength == 0)
				LastRunStart = (int) RepeatIndex;
			++RunLength;
			if (RunLength > MaxRunLength)
				{
				MaxRunLength = RunLength;
				RunStart = LastRunStart;
				RunEnd = (int) RepeatIndex;
				}
			}
		else
			RunLength = 0;
		}

	if (RunLength < g_MinArraySize)
		{
		Log("Rejected, min cons / size\n");
		LogAA(AA);
		return false;
		}

	double BestSpacerRatio = 1.0;
	for (int Base = 0; Base < (int) ArraySize - g_MinArraySize - 1; ++Base)
		{
		int MinSpacerLength = 999999;
		int MaxSpacerLength = 0;
		for (int RepeatIndex = Base; RepeatIndex < Base + g_MinArraySize; ++RepeatIndex)
			{
			assert(RepeatIndex < (int) ArraySize - 1);
			const std::string &Spacer = AA.Spacers[RepeatIndex];
			int SpacerLength = (int) Spacer.size();
			MinSpacerLength = std::min(MinSpacerLength, SpacerLength);
			MaxSpacerLength = std::max(MaxSpacerLength, SpacerLength);
			}
		double r = GetRatio(MinSpacerLength, MaxSpacerLength);
		if (r > BestSpacerRatio)
			BestSpacerRatio = r;
		}

	if (BestSpacerRatio < g_MinSpacerRatio)
		{
		Log("Rejected spacer ratio\n");
		LogAA(AA);
		return false;
		}

	return true;
	}
