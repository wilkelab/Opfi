#include "pilercr.h"

static void ExpandPile(int PileIndex)
	{
	PileData &Pile = GetModifiablePile(PileIndex);

	int Pos = Pile.Lo;
	int ContigIndex = GlobalPosToContigIndex(Pos);
	assert(ContigIndex >= 0 && ContigIndex < ContigCountQ);
	const ContigData &Contig = ContigsQ[ContigIndex];
	if (Pos < Contig.From || Pos >= Contig.From + Contig.Length)
		return;

	int ContigTo = Contig.From + Contig.Length - 1;

	Pile.Lo -= g_PileExpandMargin;
	Pile.Hi += g_PileExpandMargin;

	if (Pile.Lo < Contig.From)
		Pile.Lo = Contig.From;

	if (Pile.Hi > ContigTo)
		Pile.Hi = ContigTo;
	}

void ExpandPiles()
	{
	if (g_PileExpandMargin == 0)
		return;

	for (int PileIndex = 0; PileIndex < g_PileCount; ++PileIndex)
		ExpandPile(PileIndex);
	}
