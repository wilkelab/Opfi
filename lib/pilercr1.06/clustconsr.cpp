#include "pilercr.h"

#define TRACE	0

// ClusterConsRC:
// Cluster consensus sequences considering the sequence
// or its reverse complement, but not both.

// Best entry is one with largest number
// of distances smaller than the threshold
// Indexes into DistFunc are 2i for i'th array,
// 2i+2 for reverse complement of i'th array.
static int FindBest(std::vector<ArrayAln *> &AAVec, const DistFunc &DF,
  BoolVec &Done, IntVec &Cluster)
	{
	Cluster.clear();

	const int ArrayCount = (int) AAVec.size();
	const int Count = DF.GetCount();
	assert(Count == 2*ArrayCount);

	IntVec NrUnderThresh;
	for (int i = 0; i < Count; ++i)
		NrUnderThresh.push_back(0);

	for (int i = 0; i < Count; ++i)
		{
		if (Done[i/2])
			continue;

		for (int j = 0; j < Count; ++j)
			{
			if (i == j || Done[j/2])
				continue;
#if	TRACE
			Log("Dist(");
			AAVec[i/2]->ConsSeq.LogMeSeqOnly();
			if (i%2 == 1)
				Log("-");
			else
				Log("+");
			Log(",");
			AAVec[j/2]->ConsSeq.LogMeSeqOnly();
			if (j%2 == 1)
				Log("-");
			else
				Log("+");
			Log(") = %.1f", DF.GetDist(i, j));
#endif
			if (DF.GetDist(i, j) <= g_ClusterMaxDist)
				{
#if	TRACE
				Log(" Under\n");
#endif
				++(NrUnderThresh[i]);
				}
			else
				{
#if	TRACE
				Log(" Over\n");
#endif
				;
				}
			}
		}
#if	TRACE
	{
	Log("FindBest:\n");
	for (int i = 0; i < Count; ++i)
		Log("  i=%d Done=%c NrUnder=%d\n", i, Done[i] ? 'T' : 'F', NrUnderThresh[i]);
	}
#endif

	int Best = -1;
	unsigned BestCount = 0;
	for (int i = 0; i < Count; ++i)
		{
		if (NrUnderThresh[i] > BestCount)
			{
			Best = i;
			BestCount = NrUnderThresh[i];
			}
		}

	if (Best >= 0)
		{
		for (int i = 0; i < Count; ++i)
			{
			if (Done[i/2])
				continue;
			if (DF.GetDist(Best, i) <= g_ClusterMaxDist)
				{
				ArrayAln &AA = *(AAVec[i/2]);
				AA.ClusterRevComp = (i%2 == 1);
				Done[i/2] = true;
				Cluster.push_back(i/2);
				}
			}
		}

	if (Cluster.size() == 0)
		return -1;

	return Best;
	}

// Enforce convention that first array is + strand.
static void FixStrand(std::vector<ArrayAln *> &AAVec, IntVec &Cluster)
	{
	int First = Cluster.front();
	const ArrayAln &FirstAA = *(AAVec[First]);
	if (!FirstAA.ClusterRevComp)
		return;

	const size_t N = Cluster.size();
	for (size_t i = 0; i < N; ++i)
		{
		int ArrayIndex = Cluster[i];
		assert(ArrayIndex < (int) AAVec.size());
		ArrayAln &AA = *(AAVec[ArrayIndex]);
		AA.ClusterRevComp = !AA.ClusterRevComp;
		}
	}

void ClusterConsRC(std::vector<ArrayAln *> &AAVec, IntVecVec &Clusters)
	{
	SeqVect Seqs;

	const size_t ArrayCount = AAVec.size();
	for (size_t i = 0; i < ArrayCount; ++i)
		{
		const ArrayAln &AA = *(AAVec[i]);
		Seqs.push_back((Seq *) &AA.ConsSeq);

		Seq *r = new Seq;
		r->Copy(AA.ConsSeq);
		r->RevComp();
		Seqs.push_back(r);
		}

	DistFunc DF;
	KmerDist(Seqs, DF);

#if TRACE
	DF.LogMe();
#endif

	BoolVec Done;
	for (size_t i = 0; i < ArrayCount; ++i)
		Done.push_back(false);

	for (;;)
		{
		IntVec Cluster;
		int Best = FindBest(AAVec, DF, Done, Cluster);
		if (Best == -1)
			break;
		FixStrand(AAVec, Cluster);
		Clusters.push_back(Cluster);
		}

	for (size_t i = 0; i < ArrayCount; ++i)
		{
		if (!Done[i])
			{
			IntVec Cluster;
			Cluster.push_back(i);
			FixStrand(AAVec, Cluster);
			Clusters.push_back(Cluster);
			}
		}

#if	TRACE
	{
	Log("%d clusters\n", (int) Clusters.size());
	for (size_t i = 0; i < (int) Clusters.size(); ++i)
		{
		const IntVec &Cluster = Clusters[i];
		Log("  ");
		size_t N = Cluster.size();
		for (size_t j = 0; j < N; ++j)
			Log(" %d", Cluster[j]);
		Log("\n");
		}
	}
#endif
	}
