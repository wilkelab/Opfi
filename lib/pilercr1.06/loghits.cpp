#include "pilercr.h"

static void LocalLogHit(int HitIndex)
	{
	const DPHit &Hit = GetHit(HitIndex);

	Seq *A = new Seq;
	Seq *B = new Seq;

	GetSubseq(Hit.aLo, Hit.aHi-1, *A);
	GetSubseq(Hit.bLo, Hit.bHi-1, *B);

	A->SetName("A");
	B->SetName("B");

	SeqVect Seqs;
	Seqs.push_back(A);
	Seqs.push_back(B);

	MSA Aln;
	MultipleAlign(Seqs, Aln);

	unsigned ColCount = Aln.GetColCount();

	Log("%5d  %7d  %7d  %6d  %5d  ", HitIndex, Hit.aLo, Hit.aHi, Hit.aHi - Hit.aLo + 1,
	  GetSpacerLength(Hit));
	for (unsigned Col = 0; Col < ColCount; ++Col)
		Log("%c", Aln.GetChar(0, Col));
	Log("\n");

	Log("%5d  %7d  %7d  %6d         ", HitIndex, Hit.bLo, Hit.bHi, Hit.bHi - Hit.bLo + 1);
	for (unsigned Col = 0; Col < ColCount; ++Col)
		Log("%c", Aln.GetChar(1, Col));
	Log("\n");
	Log("\n");
	}

void LogHits()
	{
	Log("\n");
	Log("  Hit       Lo       Hi  Length  Space\n");
	Log("=====  =======  =======  ======  =====\n");
	for (int HitIndex = 0; HitIndex < g_HitCount; ++HitIndex)
		LocalLogHit(HitIndex);
	}
