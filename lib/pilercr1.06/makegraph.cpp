#include "pilercr.h"

#define TRACE	0

void MakeGraph(EdgeList &Edges)
	{
#if	TRACE
	Log("\n");
	Log("MakeGraph()\n");
	Log("  Hit  PileA  PileB\n");
	Log("=====  =====  =====\n");
#endif
	for (int HitIndex = 0; HitIndex < g_HitCount; ++HitIndex)
		{
		unsigned PileIndexA = g_HitIndexToPileIndexA[HitIndex];
		unsigned PileIndexB = g_HitIndexToPileIndexB[HitIndex];

		EdgeData Edge;
		Edge.Node1 = PileIndexA;
		Edge.Node2 = PileIndexB;
		Edge.Rev = false;
		Edges.push_back(Edge);
#if	TRACE
		Log("%5u  %5u  %5u  ", HitIndex, PileIndexA, PileIndexB);
		LogHit(HitIndex);
		Log("\n");
#endif
		}
#if	TRACE
	Log("\n");
#endif
	}
