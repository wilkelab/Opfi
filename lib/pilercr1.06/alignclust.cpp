#include "pilercr.h"

#define TRACE	0

void AlignCluster(std::vector<ArrayAln *> &AAVec, const IntVec &Cluster,
  Seq &ConsSymbols)
	{
#if	TRACE
	Log("AlignCluster:\n");
#endif

	SeqVect Seqs;
	const size_t ClusterSize  = Cluster.size();
	for (size_t i = 0; i < ClusterSize ; ++i)
		{
		size_t ArrayIndex = Cluster[i];
		Seq *s = new Seq;
#if	TRACE
		ADVec[ArrayIndex]->ConsSeq.LogMe();
		Log("\n");
#endif
		const ArrayAln &AA = *(AAVec[ArrayIndex]); 
		s->Copy(AA.ConsSeq);
		if (AA.ClusterRevComp && ClusterSize > 1)
			s->RevComp();
		Seqs.push_back(s);
		}

	MSA Aln;
	MultipleAlign(Seqs, Aln);
#if	TRACE
	Aln.LogMe();
#endif

	GetConsSymbols(Aln, ConsSymbols);

	for (size_t i = 0; i < ClusterSize ; ++i)
		{
		size_t ArrayIndex = Cluster[i];
		Aln.GetSeq((unsigned) i, AAVec[ArrayIndex]->AlignedConsSeq);
#if	TRACE
		Log("Aligned=");
		AAVec[i]->AlignedConsSeq.LogMe();
#endif
		}
	}
