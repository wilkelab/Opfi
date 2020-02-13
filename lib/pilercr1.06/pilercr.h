#ifdef _MSC_VER
#pragma warning(disable: 4996)	// deprecated functions
#define _CRT_SECURE_NO_DEPRECATE	1
#endif

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <errno.h>
#include <stdarg.h>
#include <algorithm>
#include <string>
#include "myassert.h"

#define PILERCR_LONG_VERSION	"pilercr v1.06"

#if	_MSC_VER
#pragma warning(disable:4800)	// don't warn about bool->int conversion
#endif

#ifdef _DEBUG
#define DEBUG	1
#endif

#if !defined(DEBUG) && !defined(NDEBUG)
#define	NDEBUG	1
#endif

#ifdef	WIN32
#define FILEIO_BINARY_MODE	1
#else
#define FILEIO_BINARY_MODE	0
#define stricmp strcasecmp
#define strnicmp strncasecmp
#endif

#include "types.h"
#include "progress.h"
#include "params.h"
#include "intlist.h"
#include "utils.h"
#include "multaln.h"

struct RowData
	{
	int RepeatLo;
	int RepeatLength;
	double PctId;
	};

typedef std::vector<std::string> StrVec;

struct ArrayAln
	{
	int Id;
	int Pos;	// of first repeat
	Seq ConsSeq;
	Seq AlignedConsSeq;
	bool ClusterRevComp;
	StrVec LeftFlanks;
	StrVec Repeats;
	StrVec Spacers;
	};

struct ArrayData
	{
	int Id;
	int Lo;
	int ArrayLength;
	int RepeatLength;
	int SpacerLength;
	IntVec PileIndexes;
	Seq ConsSeq;
	Seq AlignedConsSeq;
	bool ClusterRevComp;
	ArrayAln *AA;
	};

const int CONTIG_MAP_BIN_SIZE = 32;
const int MAX_CONTIGS = 1024*1024;

extern FILE *fDiscardedHits;
extern FILE *fPiles;
extern FILE *g_fOut;
extern std::vector<PileData> g_Piles;
extern int g_PileCount;
extern std::vector<ImageData> g_Images;
extern int g_ImageCount;
extern IntVec g_HitIndexToPileIndexA;
extern IntVec g_HitIndexToPileIndexB;
extern std::vector<DPHit> g_Hits;
extern int g_HitCount;

extern char *g_SeqQ;
extern int k;
extern int CharToLetter[256];
extern bool Banded;
extern int g_Diameter;
extern int SeqLengthQ;
extern ContigData *ContigsQ;
extern int *ContigMapQ;
extern int ContigCountQ;

extern unsigned g_DiscardedHitCount;
extern unsigned g_AcceptedHitCount;
extern unsigned g_NotAcceptedHitCount;

int pow4(int n);
double pow4d(int n);
double log4(double x);

void SetLog();
void Usage();
void Credits();
int GetElapsedSecs();

FILE *OpenStdioFile(const char *FileName, FILEIO_MODE Mode = FILEIO_MODE_ReadOnly);
int GetFileSize(FILE *f);

void ProcessArgVect(int argc, char *argv[]);
const char *ValueOpt(const char *Name);
void IntOpt(const char *Name, int *ptrValue);
void FloatOpt(const char *Name, double *ptrValue);
const char *RequiredValueOpt(const char *Name);
bool FlagOpt(const char *Name);
void CommandLineError(const char *Format, ...);

char *ReadMFA(const char FileName[], int *ptrLength, ContigData **ptrContigs,
  int *ptrContigCount, int **ptrContigMap);
void MakeContigMap(const ContigData *Contigs, int ContigCount, int **ptrMap);
void Complement(char *seq, int len);
void LogContigs(const ContigData *Contigs, int ContigCount);

double GetRAMSize();

int MinWordsPerFilterHit(int HitLength, int WordLength, int MaxErrors);
void GetParams(int SeqLengthQ, int g_Diameter, int Length, double MinId,
  FilterParams *ptrFP, DPParams *ptrDP);
double TotalMemRequired(int SeqLengthQ, const FilterParams &FP);
double AvgIndexListLength(int SeqLengthT, const FilterParams &FP);

void SaveFilterHit(int QFrom, int QTo, int DiagIndex);
void WriteDPHits(FILE *f, bool Comp);
void SetContigs(ContigData *ContigsT, ContigData *ContigsQ,
  int ContigCountT, int ContigCountQ, int *ContigMapT, int *ContigMapQ);
int GetFilterHitCount();
int GetFilterHitCountComp();

void Filter(int Tlen_, char *B_, int Qlen_, bool Self, bool Comp,
  const int *Finger, const int *Pos, const FilterParams &FP);
void FilterB(char *B_, int Qlen_, const FilterParams &FP);
void SetFilterOutFile(FILE *f);
void SetFilterOutFileComp(FILE *f);
FilterHit *ReadFilterHits(FILE *f, int Count);
void CloseFilterOutFile();
void CloseFilterOutFileComp();

Trapezoid *MergeFilterHits(const char *SeqT, int SeqLengthT, const char *SeqQ,
  int SeqLengthQ,  bool Self, const FilterHit *FilterHits,
  int FilterHitCount, const FilterParams &FP, int *ptrTrapCount);

void AlignTraps(char *A, int Alen, char *B, int Blen, Trapezoid *Traps,
  int TrapCount, int comp, const DPParams &DP);
int SumTrapLengths(const Trapezoid *Traps);
int SumDPLengths(const DPHit *Hits, int HitCount);

int StringToCode(const char s[], int len);
char *CodeToString(int code, int len);

void MakeIndex(char *S, int Slen, int **ptrFinger, int **ptrPos);
void CheckIndex(char *S, int Slen, const int Finger[], const int Pos[]);
void FreeIndex(int Finger[], int Pos[]);

// Memory wrappers.
// Macro hacks, but makes code more readable
// by hiding cast and sizeof.
#define	all(t, n)		(t *) allocmem((n)*sizeof(t))
#define zero(p,	t, n)	memset(p, 0, (n)*sizeof(t))
void *allocmem(int bytes);
void freemem(void *p);
unsigned GetPeakMemUseBytes();
unsigned GetMemUseBytes();
int RevCompKmer(int Kmer);
int GetKmersInWindow();

void *ckalloc(int size, const char *where);
void *ckrealloc(void *p, int size, const char *where);

void PILERCR();
bool AcceptHit(const DPHit &Hit);
void WriteDPHit(FILE *f, const DPHit &Hit, bool Comp, const char *Annot = "");
void FindArrays(const CompData &Comp, std::vector<ArrayData *> &ADVec);
double GetRatio(int x, int y);

int GetHitLength(const DPHit &Hit);
int GetHitLength(const PileData &Pile);
int GetSpacerLength(const DPHit &Hit);
int GetSpacerLength(const PileData &Pile1, const PileData &Pile2);
int GetOverlap(unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2);

inline int iabs(int x)
	{
	return x >= 0 ? x : -x;
	}
int GlobalToLocal(int Pos, char **ptrLabel);
int GlobalPosToContigIndex(int Pos);
void FindPiles();
const ContigData &GlobalPosToContig(int Pos);
bool GetArrayData(const IntVec &ArrayPileIndexes, ArrayData &AD);
void ClusterCons(const std::vector<ArrayData *> &ADVec, IntVecVec &Clusters);
void ClusterConsRC(std::vector<ArrayAln *> &AAVec, IntVecVec &Clusters);
void OutArrays(std::vector<ArrayAln *> &AAVec);
void AlignCluster(std::vector<ArrayAln *> &AAVec, const IntVec &Cluster,
  Seq &ConsSymbols);
void OutputArrayDetail(ArrayData &AD);
void WriteArrays(FILE *f, const std::vector<ArrayData *> &ADVec);
void PileToSeq(const std::vector<PileData> &Piles, int PileIndex, Seq &s,
  int *ptrPosLo, int *ptrPosHi);
void OutLocalAln(const char *A, unsigned LA, const char *B, unsigned LB,
  unsigned StartA, unsigned StartB, const char *Path_);
void LinkPiles();
int GetPileLength(const PileData &Pile);
int GetSpacerLength(const PileData &Pile1, const PileData &Pile2);
void DeletePile(int PileIndex);
void LogPile(int PileIndex);
void LogHit(const DPHit &Hit);
void LogHit(int HitIndex);
void LinkHits();
const ImageData &GetImage(int ImageIndex);
const DPHit &GetHit(int HitIndex);
const PileData &GetPile(int PileIndex);
PileData &GetModifiablePile(int PileIndex);
DPHit &GetModifiableHit(int HitIndex);
int GetHitLength(int ImageIndex);
int GetSpacerLengthFromImage(int ImageIndex);
int GetOverlap(const ImageData &Image1, const ImageData &Image2);
void HitsToImages();
void ExtendArray(IntVec &ArrayHitIndexes);
void GetBestArray(const IntVecVec &Arrays, IntVec &ArrayPileIndexes);
bool Cmp_Lo(const DPHit &Hit1, const DPHit &Hit2);
bool Cmp_Hi(const DPHit &Hit1, const DPHit &Hit2);
int GetBestLink(int FromHitIndex, const IntVec &Links);
void GetSubseq(int Lo, int Hi, Seq &s);
void LogHits();
void WritePiles();
void LogImage(const ImageData &Image);
bool IsDeletedHit(int HitIndex);
void DeleteHit(int HitIndex);
void SaveTraps(const Trapezoid *Traps);
void SaveFilterHits(const FilterHit *FHits, int FHitCount);
void WriteGFFRecord(FILE *f, int QueryFrom, int QueryTo, int TargetFrom,
  int TargetTo, bool Comp, char *Feature, int Score, char *Annot);
FILE *CreateStdioFile(const char *FileName);
void LogTraps(const Trapezoid *Traps);
int GetConnComps(EdgeList &Edges, CompList &Fams, int MinComponentSize);
void LogMx();
void ExpandPiles();
void Info();
void MakeGraph(EdgeList &Edges);
void RevComp(char *S, unsigned L);
double GetPctId(const char *A, const char *B, const char *Path);
int PathALength(const char *Path);
void AlignArray(const ArrayData &AD, ArrayAln &AA);
void FixBounds(ArrayAln &AA);
void LogAA(const ArrayAln &AA);
void ValidateAA(const ArrayAln &AA);
const int *GetCountsAA(const ArrayAln &AA, int ColIndex);
void StripBlanksAndGaps(std::string &s);
void GetConsSeqAA(const ArrayAln &AA, Seq &ConsSeq);
void StripBlanksAndGaps(Seq &s);
void OutDetailAA(const ArrayAln &AA);
int GetAvgSpacerLength(const ArrayAln &AA);
int GetAvgRepeatLength(const ArrayAln &AA);
int GetArrayLength(const ArrayAln &AA);
void AlignArrays(const std::vector<ArrayData *> &ADVec,
  std::vector<ArrayAln *> &AAVec);
int DistanceAA(const ArrayAln &AA1, const ArrayAln &AA2);
bool AcceptedAA(const ArrayAln &AA);
bool MergeAAs(ArrayAln &AA1, const ArrayAln &AA2);
int UnalignedLength(const std::string &s);
void OutPalindromes(const std::vector<ArrayAln *> &AAVec);
void SaveSeqs(FILE *f, const std::vector<ArrayAln *> AAVec);
void Cmp();
double GetFractId(const Seq &A, const Seq &B, unsigned StartA,
  unsigned StartB, std::string &Path);
void Options();
