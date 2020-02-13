#include "pilercr.h"

void StripBlanksAndGaps(std::string &s)
	{
	std::string tmp;
	for (size_t i = 0; i < s.size(); ++i)
		{
		char c = s[i];
		if (c != '-' && c != ' ')
			tmp += c;
		}
	s = tmp;
	}

void StripBlanksAndGaps(Seq &s)
	{
	Seq tmp;
	tmp.SetName(s.GetName());
	for (size_t i = 0; i < s.size(); ++i)
		{
		char c = s[i];
		if (c != '-' && c != ' ')
			tmp.push_back(c);
		}
	s = tmp;
	}

static void ValidateRep(const ArrayAln &AA, int RepeatIndex, int &Pos)
	{
	std::string Repeat = AA.Repeats[RepeatIndex];
	std::string LeftFlank = AA.LeftFlanks[RepeatIndex];
	const std::string &Spacer = AA.Spacers[RepeatIndex];

	StripBlanksAndGaps(Repeat);
	StripBlanksAndGaps(LeftFlank);

	size_t LeftFlankLength = LeftFlank.size();
	size_t RepeatLength = Repeat.size();
	size_t SpacerLength = Spacer.size();

	assert(Pos >= (int) LeftFlankLength);
	int LeftFlankStartPos = Pos - (int) LeftFlankLength;

	for (int i = 0; i < (int) LeftFlankLength; ++i)
		{
		char c1 = g_SeqQ[LeftFlankStartPos + i];
		char c2 = LeftFlank[i];
		if (toupper(c1) != toupper(c2))
			{
			Log("Repeat %d, flank[%d]=%c != Seq[%d]=%c LeftFlankStartPos=%d\n",
			  RepeatIndex, i, c1, LeftFlankStartPos + i, c2, LeftFlankStartPos);
			Log("Flank=");
			for (size_t j = 0; j < LeftFlankLength; ++j)
				Log("%c", LeftFlank[j]);
			Log("\n");
			Log("Seq[%d;%d]=", LeftFlankStartPos, LeftFlankLength);
			for (size_t j = 0; j < LeftFlankLength; ++j)
				Log("%c", g_SeqQ[LeftFlankStartPos + j]);
			Log("\n");
			LogAA(AA);
			Quit("ValidateRep failed");
			}
		}

	for (size_t i = 0; i < RepeatLength; ++i)
		{
		char c1 = g_SeqQ[Pos + i];
		char c2 = Repeat[i];
		if (toupper(c1) != toupper(c2))
			{
			Log("Repeat %d, Repeat[%d]=%c != Seq[%d]=%c\n",
			  RepeatIndex, i, c1, Pos + i, c2);
			Log("Repeat=");
			for (size_t j = 0; j < RepeatLength; ++j)
				Log("%c", Repeat[j]);
			Log("\n");
			Log("Seq[%d;%d]=", Pos, RepeatLength);
			for (size_t j = 0; j < RepeatLength; ++j)
				Log("%c", g_SeqQ[Pos + j]);
			Log("\n");
			LogAA(AA);
			Quit("ValidateRep failed");
			}
		}

	for (size_t i = 0; i < SpacerLength; ++i)
		{
		char c1 = g_SeqQ[Pos + RepeatLength + i];
		char c2 = Spacer[i];
		if (toupper(c1) != toupper(c2))
			{
			Log("Repeat %d, Spacer[%d]=%c != Seq[%d]=%c\n",
			  RepeatIndex, i, c1, Pos + RepeatLength + i, c2);
			Log("Spacer=");
			for (size_t j = 0; j < SpacerLength; ++j)
				Log("%c", Spacer[j]);
			Log("\n");
			Log("Seq[%d;%d]=", Pos + RepeatLength, SpacerLength);
			for (size_t j = 0; j < SpacerLength; ++j)
				Log("%c", g_SeqQ[Pos + RepeatLength + j]);
			Log("\n");
			LogAA(AA);
			Quit("ValidateRep failed");
			}
		}

	Pos += (int) (RepeatLength + SpacerLength);
	}

void ValidateAA(const ArrayAln &AA)
	{
	const size_t RepeatCount = AA.Spacers.size();
	assert(RepeatCount == AA.Repeats.size());
	assert(RepeatCount == AA.Spacers.size());

	int Pos = AA.Pos;
	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		const std::string &LeftFlank = AA.LeftFlanks[RepeatIndex];
		const std::string &Repeat = AA.Repeats[RepeatIndex];
		const std::string &Spacer = AA.Spacers[RepeatIndex];

		assert(Repeat.size() == AA.Repeats[0].size());
		assert(LeftFlank.size() == g_FlankSize);

		ValidateRep(AA, (int) RepeatIndex, Pos);
		}
	}
