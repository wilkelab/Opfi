#include "pilercr.h"

#define TRACE	0

// Best entry is one with largest number
// of distances smaller than the threshold
static int FindBest(const std::vector<ArrayData *> &ADVec, 
  const DistFunc &DF, BoolVec &Done, IntVec &Cluster)
	{
	Cluster.clear();

	const int Count = DF.GetCount();

	IntVec NrUnderThresh;
	for (int i = 0; i < Count; ++i)
		NrUnderThresh.push_back(0);

	for (int i = 0; i < Count; ++i)
		{
		if (Done[i])
			continue;

		for (int j = 0; j < Count; ++j)
			{
			if (i == j || Done[j])
				continue;
#if	TRACE
			Log("Dist(");
			ADVec[i]->ConsSeq.LogMeSeqOnly();
			Log(",");
			ADVec[j]->ConsSeq.LogMeSeqOnly();
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
			if (Done[i])
				continue;
			if (DF.GetDist(Best, i) <= g_ClusterMaxDist)
				{
				Done[i] = true;
				Cluster.push_back(i);
				}
			}
		}

	return Best;
	}

void ClusterCons(const std::vector<ArrayData *> &ADVec, IntVecVec &Clusters)
	{
	SeqVect Seqs;

	const size_t ArrayCount = ADVec.size();
	for (size_t i = 0; i < ArrayCount; ++i)
		{
		const ArrayData &AD = *(ADVec[i]);
		Seq *s = new Seq;
		s->Copy(AD.ConsSeq);
		Seqs.push_back(s);
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
		int Best = FindBest(ADVec, DF, Done, Cluster);
		if (Best == -1)
			break;
		Clusters.push_back(Cluster);
		}

	for (size_t i = 0; i < ArrayCount; ++i)
		{
		if (!Done[i])
			{
			IntVec Cluster;
			Cluster.push_back(i);
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
