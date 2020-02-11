#include "pilercr.h"
#include <algorithm>

#define	TRACE		0

static IntVec g_PIPos;

void LogPile(int PileIndex)
	{
	const PileData &Pile = GetPile(PileIndex);
	Log("Pile%d Lo=%d Hi=%d %s",
	  PileIndex, Pile.Lo, Pile.Hi, Pile.Deleted ? "DEL" : "");
	}

static int GetSpacer(int Index)
	{
	const int Count = (int) g_PIPos.size();
	assert(Index >= 0 && Index < Count);

	int PileIndex = g_PIPos[Index];
	const PileData &Pile = GetPile(PileIndex);

	for (int i = Index + 1; i < Count; ++i)
		{
		int NextPileIndex = g_PIPos[i];
		const PileData &NextPile = GetPile(NextPileIndex);
		if (NextPile.Deleted)
			continue;
		return NextPile.Lo - Pile.Hi;
		}
	return -1;
	}

static bool CmpPileIndexesLo(int PileIndex1, int PileIndex2)
	{
	const PileData &Pile1 = GetPile(PileIndex1);
	const PileData &Pile2 = GetPile(PileIndex2);
	return Pile1.Lo < Pile2.Lo;
	}

static bool CmpIndexesSpacer(int Index1, int Index2)
	{
	return GetSpacer(Index1) < GetSpacer(Index2);
	}

static void GetArray(int Index, IntVec &ArrayPileIndexes)
	{
#if	TRACE
	Log("GetArray(Index=%d, PileIndex=%d)\n", Index, g_PIPos[Index]);
#endif
	ArrayPileIndexes.clear();

	const int RootPileIndex = g_PIPos[Index];
	const PileData Pile = GetPile(RootPileIndex);
	if (Pile.Deleted)
		return;

	if (GetSpacer(Index) < g_DraftMinSpacerLength)
		return;

	int PileLength = GetPileLength(Pile);
	int PileSpacer = GetSpacer(Index);

	const int Count = (int) g_PIPos.size();

	ArrayPileIndexes.push_back(RootPileIndex);

	for (int PrevIndex = Index - 1; PrevIndex >= 0; --PrevIndex)
		{
		int PrevPileIndex = g_PIPos[PrevIndex];
		const PileData &PrevPile = GetPile(PrevPileIndex);
#if	TRACE
		Log("PrevPile = ");
		LogPile(PrevPileIndex);
		Log("\n");
#endif
		if (PrevPile.Deleted)
			continue;

		int PrevPileLength = GetPileLength(PrevPile);
		int PrevPileSpacer = GetSpacer(PrevIndex);

		if (PrevPileSpacer < g_DraftMinSpacerLength)
			break;

		double LengthRatio = GetRatio(PileLength, PrevPileLength);
		double SpacerRatio = GetRatio(PileSpacer + g_SpacerPadding,
		  PrevPileSpacer + g_SpacerPadding);
#if TRACE
		Log("  Spacer=%d PrevSpacer=%d LenR=%.1f SpacerR=%.1f ",
		  PileSpacer, PrevPileSpacer, LengthRatio, SpacerRatio);
#endif
		if (LengthRatio < g_DraftMinHitRatio || SpacerRatio < g_DraftMinSpacerRatio)
			{
#if	TRACE
			Log("Bad ratio\n");
#endif
			break;
			}

#if	TRACE
		Log("++ Add pile index %d\n", PrevPileIndex);
#endif
		ArrayPileIndexes.push_back(PrevPileIndex);
		}

	for (int NextIndex = Index + 1; NextIndex < Count; ++NextIndex)
		{
		int NextPileIndex = g_PIPos[NextIndex];

		if (NextIndex == Count - 1)
			{
			ArrayPileIndexes.push_back(NextPileIndex);
			break;
			}

		const PileData &NextPile = GetPile(NextPileIndex);
#if	TRACE
		Log("NextPile = ");
		LogPile(NextPileIndex);
		Log("\n");
#endif
		if (NextPile.Deleted)
			continue;

		int NextPileLength = GetPileLength(NextPile);
		int NextPileSpacer = GetSpacer(NextIndex);

		double LengthRatio = GetRatio(PileLength, NextPileLength);
		double SpacerRatio = GetRatio(PileSpacer + g_SpacerPadding,
		  NextPileSpacer + g_SpacerPadding);
#if TRACE
		Log("  Spacer=%d NextSpacer=%d LenR=%.1f SpacerR=%.1f ",
		  PileSpacer, NextPileSpacer, LengthRatio, SpacerRatio);
#endif

		if (LengthRatio < g_DraftMinHitRatio || SpacerRatio < g_DraftMinSpacerRatio)
			{
#if	TRACE
			Log("Bad ratio\n");
#endif
			break;
			}

#if	TRACE
		Log("++ Add pile index %d\n", NextPileIndex);
#endif
		ArrayPileIndexes.push_back(NextPileIndex);
		}

#if	TRACE
	Log("GetArray done, size=%d\n", (int) ArrayPileIndexes.size());
#endif
	if ((int) ArrayPileIndexes.size() < g_DraftMinArraySize)
		{
		DeletePile(RootPileIndex);
		ArrayPileIndexes.clear();
		}
	}

unsigned g_ArrayCount = 0;

void FindArrays(const CompData &Comp, std::vector<ArrayData *> &ADVec)
	{
	const int CompSize = (int) Comp.size();
	if (CompSize < g_DraftMinArraySize)
		return;

// g_PIPos is pile indexes sorted by genome pos
	g_PIPos.resize(CompSize);

	int Count = 0;
	for (CompData::const_iterator p = Comp.begin(); p != Comp.end(); ++p)
		{
		const CompMemberData &Node = *p;
		int PileIndex = Node.PileIndex;
		g_PIPos[Count] = PileIndex;
		++Count;
		}
	assert(Count == CompSize);

// Sort by genome position
	std::sort(g_PIPos.begin(), g_PIPos.end(), CmpPileIndexesLo);

#if	TRACE
	{
	Log("After sort by pos:\n");
	for (int i = 0; i < CompSize; ++i)
		{
		Log("Index=%d ", i);
		int PileIndex = g_PIPos[i];
		LogPile(PileIndex);
		Log(" Length=%d", GetPileLength(g_Piles[PileIndex]));
		Log(" Spacer=%d\n", GetSpacer(i));
		}
	}
#endif

// Sort by spacer lengths
// ISpacer is indexes into g_PIPos sorted by spacer length
	IntVec ISpacer;
	ISpacer.resize(Count);
	for (int i = 0; i < Count; ++i)
		ISpacer[i] = i;

	std::sort(ISpacer.begin(), ISpacer.end(), CmpIndexesSpacer);

	for (int i = 0; i < CompSize; ++i)
		{
		IntVec ArrayPileIndexes;
		GetArray(ISpacer[i], ArrayPileIndexes);
		if (ArrayPileIndexes.size() == 0)
			continue;

		ArrayData *AD = new ArrayData;
		AD->Id = ++g_ArrayCount;

		bool Ok = GetArrayData(ArrayPileIndexes, *AD);
		if (!Ok)
			{
			delete AD;
			continue;
			}

		ADVec.push_back(AD);

		for_CIntVec(ArrayPileIndexes, p)
			{
			int PileIndex = *p;
			DeletePile(PileIndex);
			}
		}
	}
