#ifndef TYPES_H
#define TYPES_H

#include "intlist.h"

typedef int PILE_INDEX_TYPE;

struct FilterParams
	{
	int WordSize;
	int SeedLength;
	int SeedDiffs;
	int TubeOffset;
	};

struct DPParams
	{
	int g_MinHitLength;
	double MinId;
	};

struct TubeState
	{
	int qLo;
	int qHi;
	int Count;
	};

struct ContigData
	{
	int From;
	int Length;
	int Index;
	char *Label;
	};

struct FilterHit
	{
	int QFrom;
	int QTo;
	int DiagIndex;
	};

struct Trapezoid
	{
	Trapezoid *next;	// Organized in a list linked on this field
	int top, bot;		// B (query) coords of top and bottom of trapzoidal zone
	int lft, rgt;		// Left and right diagonals of trapzoidal zone
	};

struct DPHit
	{
	int aLo, bLo;	// Start coordinate of local alignment
	int aHi, bHi;	// End coordinate of local alignment
	int ldiag, hdiag;	// Alignment is between (anti)diagonals ldiag & hdiag
	int score;			// Score of alignment where match = 1, difference = -3
	float error;		// Lower bound on error rate of match
	bool Deleted;
	IntVec LinkedHits;
	};

struct PileData
	{
	int Lo;
	int Hi;
	bool Deleted;
//	IntVec LinkedPiles;
	IntVec HitIndexes;
//	int Spacer;
	};

struct ImageData
	{
	int Lo;
	int Hi;
	int HitIndex;
	bool IsA;
	};

struct EdgeData
	{
	int Node1;
	int Node2;
	bool Rev;
	};
typedef std::list<EdgeData> EdgeList;
typedef EdgeList::iterator PtrEdgeList;

struct CompMemberData
	{
	int PileIndex;
	bool Rev;
	};
typedef std::list<CompMemberData> CompData;
typedef CompData::iterator PtrFamData;

typedef std::list<CompData *> CompList;
typedef CompList::iterator PtrCompList;

enum FILEIO_MODE
	{
	FILEIO_MODE_Undefined = 0,
	FILEIO_MODE_ReadOnly,
	FILEIO_MODE_ReadWrite,
	FILEIO_MODE_WriteOnly
	};

#if	FILEIO_BINARY_MODE
#define FILIO_STRMODE_ReadOnly		"rb"
#define FILIO_STRMODE_WriteOnly		"wb"
#define FILIO_STRMODE_ReadWrite		"w+b"
#else
#define FILIO_STRMODE_ReadOnly		"r"
#define FILIO_STRMODE_WriteOnly		"w"
#define FILIO_STRMODE_ReadWrite		"w+"
#endif

struct INDEX_ENTRY
	{
	int Kmer;	// for debugging only
	int Pos;
	INDEX_ENTRY *Next;
	INDEX_ENTRY *Prev;
	};

enum TERMGAPS
	{
	TERMGAPS_Full,
	TERMGAPS_Half,
	TERMGAPS_Ext,
	};

#endif // TYPES_H
