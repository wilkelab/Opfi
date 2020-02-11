#include "pilercr.h"

static void WriteArray(FILE *f, const ArrayData &AD)
	{
	const int PileIndexCount = (int) AD.PileIndexes.size();
	for (int i = 0; i < PileIndexCount; ++i)
		{
		int PileIndex = AD.PileIndexes[i];
		assert(PileIndex >= 0 && PileIndex < g_PileCount);
		const PileData &Pile = g_Piles[PileIndex];

		char *Label;
		int Lo = GlobalToLocal(Pile.Lo, &Label);
		int Hi = Lo + Pile.Hi - Pile.Lo + 1;

	// GFF Fields are:
	// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
	//     0         1         2        3      4      5        6       7         8           9
		fprintf(f, "%s\tcrisper\tpile\t%d\t%d\t0\t+\t.\tPile %d ; Array %d\n",
		  Label,
		  Lo + 1,
		  Hi + 1,
		  PileIndex,
		  AD.Id);
		}
	}

void WriteArrays(FILE *f, const std::vector<ArrayData *> &ADVec)
	{
	const size_t ArrayCount = ADVec.size();
	for (size_t ArrayIndex = 0; ArrayIndex < ArrayCount; ++ArrayIndex)
		{
		const ArrayData &AD = *(ADVec[ArrayIndex]);
		WriteArray(f, AD);
		}
	}
