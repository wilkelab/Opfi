#include "pilercr.h"

void WriteGFFRecord(FILE *f, int QueryFrom, int QueryTo, int TargetFrom,
  int TargetTo, bool Comp, char *Feature, int Score, char *Annot)
	{
	if (Comp)
		{
// Another lazy hack: empirically found off-by-one
// complemented hits have query pos off-by-one, so
// fix by this change. Should figure out what's really
// going on.
//		QueryFrom = SeqLengthQ - Hit.bHi - 1;
//		QueryTo = SeqLengthQ - Hit.bLo - 1;
		QueryFrom = SeqLengthQ - QueryTo;
		QueryTo = SeqLengthQ - QueryFrom;
		if (QueryFrom >= QueryTo)
			Quit("WriteGFFRecord: QueryFrom >= QueryTo (comp)");
		}
	else
		{
		if (QueryFrom >= QueryTo)
			Quit("WriteGFFRecord: QueryFrom >= QueryTo (not comp)");
		}

// DPHit coordinates sometimes over/underflow.
// This is a lazy hack to work around it, should really figure
// out what is going on.
	if (QueryFrom < 0)
		QueryFrom = 0;
	if (QueryTo >= SeqLengthQ)
		QueryTo = SeqLengthQ - 1;
	if (TargetFrom < 0)
		TargetFrom = 0;
	if (TargetTo >= SeqLengthQ)
		TargetTo = SeqLengthQ - 1;

// Take midpoint of segment -- lazy hack again, endpoints
// sometimes under / overflow
	const int TargetBin = (TargetFrom + TargetTo)/(2*CONTIG_MAP_BIN_SIZE);
	const int QueryBin = (QueryFrom + QueryTo)/(2*CONTIG_MAP_BIN_SIZE);
	const int BinCountT = (SeqLengthQ + CONTIG_MAP_BIN_SIZE - 1)/CONTIG_MAP_BIN_SIZE;
	const int BinCountQ = (SeqLengthQ + CONTIG_MAP_BIN_SIZE - 1)/CONTIG_MAP_BIN_SIZE;
	if (TargetBin < 0 || TargetBin >= BinCountT)
		{
		Warning("Target bin out of range");
		return;
		}
	if (QueryBin < 0 || QueryBin >= BinCountQ)
		{
		Warning("Query bin out of range");
		return;
		}

	if (ContigMapQ == 0)
		Quit("ContigMap = 0");

	const int TargetContigIndex = ContigMapQ[TargetBin];
	const int QueryContigIndex = ContigMapQ[QueryBin];

	if (TargetContigIndex < 0 || TargetContigIndex >= ContigCountQ)
		Quit("WriteGFFRecord: bad target contig index");

	if (QueryContigIndex < 0 || QueryContigIndex >= ContigCountQ)
		Quit("WriteGFFRecord: bad query contig index");

	const int TargetLength = TargetTo - TargetFrom + 1;
	const int QueryLength = QueryTo - QueryFrom + 1;

	if (TargetLength < 0 || QueryLength < 0)
		Warning("WriteGFFRecord: Length < 0");

	const ContigData &TargetContig = ContigsQ[TargetContigIndex];
	const ContigData &QueryContig = ContigsQ[QueryContigIndex];

	int TargetContigFrom = TargetFrom - TargetContig.From + 1;
	int QueryContigFrom = QueryFrom - QueryContig.From + 1;

	int TargetContigTo = TargetContigFrom + TargetLength - 1;
	int QueryContigTo = QueryContigFrom + QueryLength - 1;

	if (TargetContigFrom < 1)
		TargetContigFrom = 1;
	if (TargetContigTo > TargetContig.Length)
		TargetContigTo = TargetContig.Length;
	if (QueryContigFrom < 1)
		QueryContigFrom = 1;
	if (QueryContigTo > QueryContig.Length)
		QueryContigTo = QueryContig.Length;

	const char *TargetLabel = TargetContig.Label;
	const char *QueryLabel = QueryContig.Label;

	const char Strand = Comp ? '-' : '+';

// GFF Fields are:
// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
//     0         1         2        3      4      5        6       7         8           9
	fprintf(f, "%s\tcrisper\t%s\t%d\t%d\t%d\t%c\t.\tTarget %s %d %d",
	  QueryLabel,
	  Feature,
	  QueryContigFrom,
	  QueryContigTo,
	  Score,
	  Strand,
	  TargetLabel,
	  TargetContigFrom,
	  TargetContigTo);

	if (Annot != 0)
		fprintf(f, "; %s", Annot);

	fprintf(f, "\n");
	}
