#include "pilercr.h"

void WritePiles()
	{
	const char *FileName = ValueOpt("piles");
	if (FileName == 0)
		return;
	FILE *f = OpenStdioFile(FileName, FILEIO_MODE_WriteOnly);

	ProgressStart("Writing piles");
	for (int PileIndex = 0; PileIndex < g_PileCount; ++PileIndex)
		{
		if (PileIndex%100 == 0)
			ProgressStep(PileIndex, g_PileCount);

		PileData &Pile = g_Piles[PileIndex];
		IntVec &HitIndexes = Pile.HitIndexes;
		for_CIntVec(HitIndexes, p)
			{
			int HitIndex = *p;
			const DPHit &Hit = GetHit(HitIndex);
			char Annot[128];
			sprintf(Annot, "Hit(%d) ; Pile(%d)", HitIndex, PileIndex);
			WriteDPHit(f, Hit, false, Annot);
			}
		}
	ProgressDone();

	fclose(f);
	}
