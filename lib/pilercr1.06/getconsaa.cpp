#include "pilercr.h"

const int *GetCountsAA(const ArrayAln &AA, int ColIndex)
	{
	static int Counts[5];

	const size_t RepeatCount = AA.Repeats.size();

	for (int i = 0; i < 5; ++i)
		Counts[i] = 0;

	for (size_t RepeatIndex = 0; RepeatIndex < RepeatCount; ++RepeatIndex)
		{
		const std::string &Repeat = AA.Repeats[RepeatIndex];
		assert(ColIndex >= 0 && ColIndex < (int) Repeat.size());
		char c = Repeat[ColIndex];
		switch (toupper(c))
			{
		case 'A':
			++(Counts[0]);
			break;
		case 'C':
			++(Counts[1]);
			break;
		case 'G':
			++(Counts[2]);
			break;
		case 'T':
			++(Counts[3]);
			break;
		case '-':
			++(Counts[4]);
			break;
			}
		}
	return Counts;
	}

void GetConsSeqAA(const ArrayAln &AA, Seq &ConsSeq)
	{
	ConsSeq.clear();
	ConsSeq.SetName("cons");

	size_t ColCount = AA.Repeats[0].size();
	for (size_t ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		const int *Counts = GetCountsAA(AA, (int) ColIndex);
		int MaxCount = 0;
		int MaxLetter = 0;
		for (int i = 0; i < 5; ++i)
			{
			if (Counts[i] > MaxCount)
				{
				MaxCount = Counts[i];
				MaxLetter = i;
				}
			}

		char c = '?';
		switch (MaxLetter)
			{
		case 0:
			c = 'A';
			break;
		case 1:
			c = 'C';
			break;
		case 2:
			c = 'G';
			break;
		case 3:
			c = 'T';
			break;
		case 4:
			c = '-';
			break;
		default:
			Quit("Bad maxletter");
			}

		ConsSeq.push_back(c);
		}
	}
