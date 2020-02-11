#include "pilercr.h"

void MakeContigMap(const ContigData *Contigs, int ContigCount, int **ptrMap)
	{
	*ptrMap = 0;
	if (ContigCount <= 0)
		return;
	const ContigData &LastContig = Contigs[ContigCount-1];
	const int SeqLength = LastContig.From + LastContig.Length - 1;
	const int BinCount = (SeqLength + CONTIG_MAP_BIN_SIZE - 1)/CONTIG_MAP_BIN_SIZE + 1;
	int *Map = all(int, BinCount);

// Initialize, enables correctness check
	for (int i = 0; i < BinCount; ++i)
		Map[i] = -1;

	for (int ContigIndex = 0; ContigIndex < ContigCount; ++ContigIndex)
		{
		const ContigData &Contig = Contigs[ContigIndex];

	// Contig required to start on bin boundary
		const int From = Contig.From;
		if (From%CONTIG_MAP_BIN_SIZE)
			Quit("MakeContigMap: Contig does not start on bin boundary");

		const int To = From + Contig.Length - 1;
		const int BinFrom = From/CONTIG_MAP_BIN_SIZE;
		const int BinTo = To/CONTIG_MAP_BIN_SIZE;

		for (int Bin = BinFrom; Bin <= BinTo; ++Bin)
			{
			if (-1 != Map[Bin])
				Quit("MakeContigMap logic error 1");
			Map[Bin] = ContigIndex;
			}
		}

	*ptrMap = Map;
	for (int ContigIndex = 0; ContigIndex < ContigCount; ++ContigIndex)
		{
		const ContigData &Contig = Contigs[ContigIndex];
		int ContigIndex2 = GlobalPosToContigIndex(Contig.From);
		if (ContigIndex2 != ContigIndex)
			{
			LogContigs(Contigs, ContigCount);
			Quit("MakeContigMap: start pos %d of contig %d mapped to contig %d",
			  Contig.From, ContigIndex, ContigIndex2);
			}

		int To = Contig.From + Contig.Length - 1;
		int ContigIndex3 = GlobalPosToContigIndex(To);
		if (ContigIndex2 != ContigIndex)
			{
			LogContigs(Contigs, ContigCount);
			Quit("MakeContigMap: end pos %d of contig %d mapped to contig index %d",
			  Contig.From, ContigIndex, ContigIndex3);
			}
		}
	}

const ContigData &GlobalPosToContig(int Pos)
	{
	int ContigIndex = GlobalPosToContigIndex(Pos);
	assert(ContigIndex >= 0 && ContigIndex < ContigCountQ);
	const ContigData &CD = ContigsQ[ContigIndex];
	if (Pos < CD.From || Pos >= CD.From + CD.Length)
		{
		LogContigs(ContigsQ, ContigCountQ);
		Quit("GlobalPosToContig(%d): returned contig %d, out of range",
		  Pos, ContigIndex);
		}
	return CD;
	}

int GlobalPosToContigIndex(int Pos)
	{
	int Bin = Pos/CONTIG_MAP_BIN_SIZE;
	const int BinCount = (SeqLengthQ + CONTIG_MAP_BIN_SIZE - 1)/CONTIG_MAP_BIN_SIZE;

// Hack to avoid crashes
	if (Bin < 0)
		Bin = 0;
	if (Bin >= BinCount)
		Bin = BinCount;

	if (ContigMapQ == 0)
		Quit("ContigMap = 0");

	int ContigIndex = -1;
// Awful hack...
	for (;;)
		{
		ContigIndex = ContigMapQ[Bin];
		if (ContigIndex >= 0)
			break;

		if (Bin == 0)
			Quit("GlobalPosToContigIndex");
		--Bin;
		}

	assert(ContigIndex >= 0 && ContigIndex < ContigCountQ);
	return ContigIndex;
	}

int GlobalToLocal(int Pos, char **ptrLabel)
	{
	int ContigIndex = GlobalPosToContigIndex(Pos);
	assert(ContigIndex >= 0 && ContigIndex < ContigCountQ);
	const ContigData &Contig = ContigsQ[ContigIndex];
	*ptrLabel = Contig.Label;
	return Pos - Contig.From;
	}
