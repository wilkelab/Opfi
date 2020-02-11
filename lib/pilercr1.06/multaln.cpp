#include "pilercr.h"
#include "multaln.h"

void MultipleAlign(SeqVect &Seqs, MSA &Aln)
	{
	assert(Seqs.GetSeqCount() > 0);
	if ((int) Seqs.GetSeqCount() > g_MaxMSASeqs)
		Quit("Array length %u, max is %d", Seqs.GetSeqCount(), g_MaxMSASeqs);

	if (Seqs.GetSeqCount() == 1)
		{
		Seq &s = Seqs.GetSeq(0);
		s.SetId(0);
		Aln.FromSeq(s);
		return;
		}

	unsigned uSeqCount = Seqs.GetSeqCount();
	MSA::SetIdCount(g_MaxMSASeqs);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		Seqs.SetSeqId(uSeqIndex, uSeqIndex);

	Tree GuideTree;
	GetGuideTree(Seqs, GuideTree);
	ProgressiveAlign(Seqs, GuideTree, Aln);
	}
