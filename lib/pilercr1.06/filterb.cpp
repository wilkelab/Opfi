#include "pilercr.h"
#include "forkmer.h"

// FilterB: Banded search of SeqQ
// (There is no SeqT).

#define	TRACE		0

// static char *SeqT;
static char *SeqQ;
static int Qlen;
static int MinMatch;
static int MaxError;
static int TubeOffset;
static int TubeWidth;
static int MinKmersPerHit;
static TubeState *Tubes;
static int MaxActiveTubes;
static int MaxKmerDist;

//==============================================================
// Start rolling index stuff
//==============================================================

static INDEX_ENTRY *Entries;
static INDEX_ENTRY **Heads;
static INDEX_ENTRY **Tails;
static INDEX_ENTRY *FreeEntries;
static int KmerIndexCount;
static int KmersInWindow;

static int FreeListSize()
	{
	int Size = 0;
	for (INDEX_ENTRY *I = FreeEntries; I; I = I->Next)
		++Size;
	return Size;
	}

static void AddToFreeList(INDEX_ENTRY *IE)
	{
#if	DEBUG
	IE->Kmer = -2;
	IE->Pos = -2;
#endif

	IE->Prev = 0;
	IE->Next = FreeEntries;
	if (0 != FreeEntries)
		FreeEntries->Prev = IE;
	FreeEntries = IE;
#if	TRACE
	Log("AddToFreeList(%x), size when done=%d\n", IE, FreeListSize());
#endif
	}

int GetKmersInWindow()
	{
	return g_Diameter - k + 1;
	}

static void AllocateIndex(int g_Diameter)
	{
	assert(KmersInWindow > 0);

	KmerIndexCount = pow4(k);

	int EntryCount = KmersInWindow + k;
	Entries = all(INDEX_ENTRY, EntryCount);
	zero(Entries, INDEX_ENTRY, EntryCount);

	Heads = all(INDEX_ENTRY *, KmerIndexCount);
	Tails = all(INDEX_ENTRY *, KmerIndexCount);

	zero(Heads, INDEX_ENTRY *, KmerIndexCount);
	zero(Tails, INDEX_ENTRY *, KmerIndexCount);

	for (int i = 0; i < EntryCount; ++i)
		AddToFreeList(&(Entries[i]));

#if	TRACE
	Log("KmersInWindow = g_Diameter=%d - k=%d + 1 = %d\n",
	  g_Diameter, k, KmersInWindow);
	Log("After AllocateIndex free list size = %d\n", FreeListSize());
#endif
	}

static void AddToIndex(int Kmer, int Pos)
	{
#if	TRACE
	Log("AddToIndex(Kmer=%x=%s, Pos=%d)\n", Kmer, CodeToString(Kmer, k), Pos);
#endif

	if (-1 == Kmer)
		return;

#if	TRACE
	Log("Free index has %d entries\n", FreeListSize());
#endif

	assert(Kmer >= 0 && Kmer < KmerIndexCount);
	assert(FreeEntries != 0);

	INDEX_ENTRY *Entry = FreeEntries;
	FreeEntries = FreeEntries->Next;
	if (FreeEntries != 0)
		FreeEntries->Prev = 0;

	INDEX_ENTRY *Tail = Tails[Kmer];

// Insert after tail of list
	Entry->Kmer = Kmer;
	Entry->Pos = Pos;
	Entry->Next = 0;
	Entry->Prev = Tail;

	if (0 == Tail)
		Heads[Kmer] = Entry;
	else
		Tail->Next = Entry;
	Tails[Kmer] = Entry;
	}

static inline int GetKmer(const char *Seq, int Pos)
	{
	assert(Pos + k <= SeqLengthQ);
	return StringToCode(Seq + Pos, k);
	}

#define DeclareListPtr(p)	INDEX_ENTRY *p

static inline INDEX_ENTRY *GetListPtr(int Kmer)
	{
	assert(Kmer >= 0 && Kmer < KmerIndexCount);
	return Heads[Kmer];
	}

static inline bool NotEndOfList(INDEX_ENTRY *p)
	{
	return p != 0;
	}

static inline INDEX_ENTRY *GetListNext(INDEX_ENTRY *p)
	{
	return p->Next;
	}

static inline int GetListPos(INDEX_ENTRY *p)
	{
	return p->Pos;
	}

static void ValidateIndex(const char *Seq, int WindowStart, int WindowEnd)
	{
	ProgressStart("Validating index");
	int FreeCount = 0;
	for (INDEX_ENTRY *p = FreeEntries; p != 0; p = p->Next)
		{
		if (++FreeCount > KmersInWindow)
			Quit("Validate index failed free count");

		if (p->Kmer != -2 || p->Pos != -2)
			Quit("Validate index failed free != -2");

		const INDEX_ENTRY *pNext = p->Next;
		if (0 != pNext && pNext->Prev != p)
			Quit("Validate index failed free pNext->Prev != p");

		const INDEX_ENTRY *pPrev = p->Prev;
		if (0 != pPrev && pPrev->Next != p)
			Quit("Validate index failed free pPrev->Next != p");
		}

	for (int Pos = WindowStart; Pos < WindowEnd; ++Pos)
		{
		const int Kmer = GetKmer(Seq, Pos);
		if (-1 == Kmer)
			continue;
		for (DeclareListPtr(p) = GetListPtr(Kmer); NotEndOfList(p); p = GetListNext(p))
			{
			const int HitPos = GetListPos(p);
			if (HitPos == Pos)
				goto Found;
			}
		Quit("Validate index failed, pos=%d not found for kmer %x=%s",
		  Pos, Kmer, CodeToString(Kmer, k));
	Found:;
		}

	int IndexedCount = 0;
	for (int Kmer = 0; Kmer < KmerIndexCount; ++Kmer)
		{
		INDEX_ENTRY *Head = Heads[Kmer];
		INDEX_ENTRY *Tail = Tails[Kmer];
		if (Head != 0 && Head->Prev != 0)
			Quit("Head->Prev != 0");
		if (Tail != 0 && Tail->Next != 0)
			Quit("Tail->Next != 0");
		if ((Head == 0) != (Tail == 0))
			Quit("Head / tail");
		int PrevHitPos = -1;
		int ListIndex = 0;
		for (DeclareListPtr(p) = GetListPtr(Kmer); NotEndOfList(p); p = GetListNext(p))
			{
			++IndexedCount;
			if (IndexedCount > KmersInWindow)
				Quit("Valiate index failed, count");

			const INDEX_ENTRY *pNext = p->Next;
			if (Kmer != p->Kmer)
				Quit("Validate index failed, kmer");

			if (0 != pNext && pNext->Prev != p)
				Quit("Validate index failed pNext->Prev != p");

			const INDEX_ENTRY *pPrev = p->Prev;
			if (0 != pPrev && pPrev->Next != p)
				Quit("Validate index failed pPrev->Next != p");

			const int HitPos = GetListPos(p);
			if (HitPos < WindowStart || HitPos > WindowEnd)
				Quit("ValidateIndex failed, hit not in window kmer=%d %s",
				  Kmer, CodeToString(Kmer, k));

			int IsTail = (p->Next == 0);
			if (HitPos < PrevHitPos)
				Quit("Validate index failed, sort order Kmer=%d HitPos=%d PrevHitPos=%d ListIndex=%d IsTail=%d",
				  Kmer, HitPos, PrevHitPos, ListIndex, IsTail);

			PrevHitPos = HitPos;
			++ListIndex;
			}
		}
	if (IndexedCount > KmersInWindow)
		Quit("Validate index failed, count [2]");

	ProgressDone();
	}

static void LogLocations(int Kmer)
	{
	Log("LogLocations(%d %s)", Kmer, CodeToString(Kmer, k));
	for (DeclareListPtr(p) = GetListPtr(Kmer); NotEndOfList(p); p = GetListNext(p))
		Log(" [%d]=%d", p->Pos, StringToCode(SeqQ + p->Pos, k));
	}

// Pos not required; used for sanity check
static void DeleteFirstInstanceFromIndex(int Kmer, int Pos)
	{
#if	TRACE
	Log("DeleteFirstInstanceFromIndex(Kmer=%x=%s, Pos=%d)\n",
	  Kmer, CodeToString(Kmer, k), Pos);
#endif

	if (-1 == Kmer)
		return;

	assert(Kmer >= 0 && Kmer < KmerIndexCount);

	INDEX_ENTRY *IE = Heads[Kmer];

	if (IE == 0)
		Quit("DFI Kmer=%d %s Pos=%d", Kmer, CodeToString(Kmer, k), Pos);
//	assert(IE != 0);
	assert(0 == IE->Prev);
	assert(Pos == IE->Pos);

// Delete from index
	INDEX_ENTRY *NewHead = IE->Next;
	if (NewHead == 0)
		{
		Heads[Kmer] = 0;
		Tails[Kmer] = 0;
		}
	else
		{
		assert(NewHead->Prev == IE);
		NewHead->Prev = 0;
		}
	Heads[Kmer] = NewHead;

	AddToFreeList(IE);
	}

//==============================================================
// End rolling index stuff
//==============================================================

#define CalcDiagIndex(t, q)			(Qlen - (t) + (q))
#define CalcTubeIndex(DiagIndex)	((DiagIndex)/TubeOffset)

static void AddHit(int TubeIndex, int qLo, int qHi)
	{
#if	TRACE
	Log("AddHit(Tube=%d, qLo=%d, qHi=%d)\n",
	  Qlen - TubeIndex*TubeOffset, qLo, qHi + k);
#endif

	SaveFilterHit(qLo, qHi + k, Qlen - TubeIndex*TubeOffset);
	}

// Called when end of a tube is reached
// A point in the tube -- the point with maximal q -- is (Qlen-1,q-1).
static void TubeEnd(int q)
	{
	if (q <= 0)
		return;
	int DiagIndex = CalcDiagIndex(Qlen - 1, q - 1);
	int TubeIndex = CalcTubeIndex(DiagIndex);

	TubeState *Tube = Tubes + TubeIndex%MaxActiveTubes;
#if	TRACE
	Log("TubeEnd(%d) DiagIndex=%d TubeIndex=%d Count=%d\n",
	  q, DiagIndex, TubeIndex, Tube->Count);
#endif
	if (Tube->Count >= MinKmersPerHit)
		AddHit(TubeIndex, Tube->qLo, Tube->qHi);

	Tube->Count = 0;
	}

// Called when q=Qlen - 1 to flush any hits in each tube.
static void TubeFlush(int TubeIndex)
	{
	TubeState *Tube = Tubes + TubeIndex%MaxActiveTubes;
#if	TRACE
	Log("TubeFlush(TubeIndex=%d) Count=%d\n",
	  TubeIndex, Tube->Count);
#endif
	if (Tube->Count < MinKmersPerHit)
		return;

	AddHit(TubeIndex, Tube->qLo, Tube->qHi);
	Tube->Count = 0;
	}

static void HitTube(int TubeIndex, int q)
	{
	TubeState *Tube = Tubes + TubeIndex%MaxActiveTubes;

#if	TRACE
	Log("HitTube(TubeIndex=%d, q=%d) Count=%d\n",
	  TubeIndex, q, Tube->Count);
#endif

	if (0 == Tube->Count)
		{
		Tube->Count = 1;
		Tube->qLo = q;
		Tube->qHi = q;
		return;
		}

	if (q - Tube->qHi > MaxKmerDist)
		{
		if (Tube->Count >= MinKmersPerHit)
			AddHit(TubeIndex, Tube->qLo, Tube->qHi);

		Tube->Count = 1;
		Tube->qLo = q;
		Tube->qHi = q;
		return;
		}

	++(Tube->Count);
	Tube->qHi = q;
	}

// Found a common k-mer
static inline void CommonKmer(int t, int q)
	{
	assert(t >= 0 && t < Qlen - k + 1);
	assert(q >= 0 && q < Qlen - k + 1);

	if (q <= t)
		return;

#if	TRACE
	Log("CommonKmer(%d,%d) SeqQ=%.*s SeqQ=%.*s\n",
	  q, t, k, SeqQ+q, k, SeqQ+t);
#endif

	int DiagIndex = CalcDiagIndex(t, q);
	int TubeIndex = CalcTubeIndex(DiagIndex);

#if	TRACE
	Log("HitTube(TubeIndex=%d, t=%d, q=%d)\n", TubeIndex, t, q);
#endif
	HitTube(TubeIndex, q);

// Hit in overlapping tube preceding this one?
	if (DiagIndex%TubeOffset < MaxError)
		{
		if (0 == TubeIndex)
			TubeIndex = MaxActiveTubes - 1;
		else
			--TubeIndex;
		assert(TubeIndex >= 0);
#if	TRACE
		Log("HitTube(TubeIndex=%d, t=%d, q=%d) [overlap]\n", TubeIndex, t, q);
#endif
		HitTube(TubeIndex, q);
		}
	}

void FilterBold(char *B_, int Qlen_, const FilterParams &FP)
	{
	SeqQ = B_;
	Qlen = Qlen_;

	MinMatch = FP.SeedLength;
	MaxError = FP.SeedDiffs;
	TubeOffset = FP.TubeOffset;

	const int Kmask = pow4(k) - 1;

// Ukonnen's Lemma
	MinKmersPerHit = MinMatch + 1 - k*(MaxError + 1);

// Maximum distance between SeqQ positions of two k-mers in a match
// (More stringent bounds may be possible, but not a big problem
// if two adjacent matches get merged).
	MaxKmerDist = MinMatch - k;

	TubeWidth = TubeOffset + MaxError;

	if (TubeOffset < MaxError)
		{
		Log("TubeOffset < MaxError\n");
		exit(1);
		}
	if (MinKmersPerHit <= 0)
		{
		Log("MinKmersPerHit <= 0\n");
		exit(1);
		}

	MaxActiveTubes = (Qlen + TubeWidth - 1)/TubeOffset + 1;
	Tubes = all(TubeState, MaxActiveTubes);
	zero(Tubes, TubeState, MaxActiveTubes);

// Ticker tracks cycling of circular list of active tubes.
	int Ticker = TubeWidth;

// Allocate memory for index
	ProgressStart("Allocating index");
	KmersInWindow = GetKmersInWindow();
	AllocateIndex(g_Diameter);
	ProgressDone();

/***
Scan query sequence.
Start is coordinate of first base in first k-mer in sliding window
End is coordinate of first base in last k-mer in sliding window

                   g_Diameter=10
                -----------------
               |                 |
               v                 v

             Start          End = Start + g_Diameter - k
               v____         v____
    |A|C|T|G|G|A|T|C|C|G|A|T|T|A|C|A|C|C|T|T|A|G|T|C|A|   SeqQ -->
     ^---^     ^___^             .
      k=3        ^___^           .
                   ^___^         . KmersInWindow =
                     ^___^       .   g_Diameter - k + 1
                       ^___^     .   = 8
                         ^___^   .
                           ^___^ .
                             ^___^
***/

	ProgressStart("Filtering");

// Position of last kmer in sequence
	const int LastKmerPos = Qlen - k;

// Find first valid kmer in window
	int FirstValidKmerPos = 0;
	int FirstValidKmer = 0;
	for (;;)
		{
		if (FirstValidKmerPos >= Qlen)
			Quit("No valid words in query");

		FirstValidKmer = GetKmer(SeqQ, FirstValidKmerPos);
		if (FirstValidKmer != -1)
			break;

		++FirstValidKmerPos;
		}
	AddToIndex(FirstValidKmer, FirstValidKmerPos);

// Last valid kmer in window
	int LastValidKmerPos = FirstValidKmerPos;

// Positions of first and last kmer in window
	int Start = -g_Diameter + k;
	int End = 0;

	int StartKmer = FirstValidKmer;
	int EndKmer = FirstValidKmer;

#if	TRACE
	Log("Start             %d\n", Start);
	Log("End               %d\n", End);
	Log("FirstValidKmer    %s\n", CodeToString(FirstValidKmer, k));
	Log("FirstValidKmerPos %d\n", FirstValidKmerPos);
	Log("LastKmerPos       %d\n", LastKmerPos);
#endif

	for (; Start <= LastKmerPos; ++Start, ++End)
		{
#if	TRACE
		if (Start < 0)
			Log("START  < 0\n");
		else
			Log("START  %4d  %.*s\n", Start, k, SeqQ+Start);
		if (End > LastKmerPos)
			Log("END    > End\n");
		else
			Log("END    %4d  %.*s\n", End, k, SeqQ+End);
#endif
		if (Start%1000000 == 0)
			ProgressStep(End, LastKmerPos);

		if (Start >= FirstValidKmerPos)
			{
			assert(StartKmer == GetKmer(SeqQ, Start));
			for (DeclareListPtr(p) = GetListPtr(StartKmer); NotEndOfList(p); p = GetListNext(p))
				{
				int HitPos = GetListPos(p);
				CommonKmer(Start, HitPos);
				}
			DeleteFirstInstanceFromIndex(StartKmer, Start);
			}

		if (0 == --Ticker)
			{
			TubeEnd(Start);
			Ticker = TubeOffset;
			}

		{
		if (Start >= FirstValidKmerPos)
			{
			char c = SeqQ[Start + k];
			int x = CharToLetter[c];
			if (x < 0)
				{
				StartKmer = 0;
				FirstValidKmerPos = Start + k + 1;
				}
			else
				StartKmer = ((StartKmer << 2) | x) & Kmask;
			}
		}

		{
		if (End <= LastKmerPos)
			{
			if (End >= LastValidKmerPos)
				{
#if	DEBUG
				if (EndKmer != GetKmer(SeqQ, End))
					{
					Log("EndKmer=%s", CodeToString(EndKmer, k));
					Log(" != GetKmer(End=%d) = %s\n",
					  End, CodeToString(GetKmer(SeqQ, End), k));
					Quit("EndKmer != GetKmer");
					}
#endif
				AddToIndex(EndKmer, End);
				}

			char c = SeqQ[End + k];
			int x = CharToLetter[c];
			if (x < 0)
				{
				EndKmer = 0;
				LastValidKmerPos = End + k + 1;
				}
			else
				EndKmer = ((EndKmer << 2) | x) & Kmask;
			}
		}
		}
	ProgressDone();

	TubeEnd(Qlen - 1);

	int DiagFrom = CalcDiagIndex(Qlen - 1, Qlen - 1) - TubeWidth;
	int DiagTo = CalcDiagIndex(0, Qlen - 1) + TubeWidth;

	int TubeFrom = CalcTubeIndex(DiagFrom);
	if (TubeFrom < 0)
		TubeFrom = 0;

	int TubeTo = CalcTubeIndex(DiagTo);

	for (int TubeIndex = TubeFrom; TubeIndex <= TubeTo; ++TubeIndex)
		TubeFlush(TubeIndex);

	freemem(Tubes);
	}
