#include "multaln.h"

void GetGuideTree(const SeqVect &Seqs, Tree &GuideTree)
	{
	DistFunc DF;
	KmerDist(Seqs, DF);

	DistCalcDF DC;
	DC.Init(DF);

	UPGMA(DC, GuideTree);
	}
