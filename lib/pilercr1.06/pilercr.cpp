#include "pilercr.h"

// for getpid:
#if	WIN32
#include <process.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif

int k;
char *g_SeqQ = 0;
ContigData *ContigsQ = 0;
int SeqLengthQ = -1;
int *ContigMapQ = 0;
int ContigCountQ = -1;
bool FilterComp = false;
int g_MinHitLength = 16;
static double g_MinId = 0.94;
unsigned g_AcceptedHitCount = 0;
unsigned g_NotAcceptedHitCount = 0;
std::vector<PileData> g_Piles;
int g_PileCount;
std::vector<DPHit> g_Hits;
int g_HitCount;
IntVec g_HitIndexToPileIndexA;
IntVec g_HitIndexToPileIndexB;
std::vector<ImageData> g_Images;
int g_ImageCount;
int g_Diameter;
extern bool g_ShowProgress;

FILE *fDiscardedHits = 0;
FILE *fPiles = 0;
FILE *g_fOut = 0;
FILE *g_fSeqs = 0;

unsigned g_DiscardedHitCount = 0;

static const char *GetFilterOutFileName(const char *FilterOutPrefix, const char *Ext)
	{
	if (0 == FilterOutPrefix)
		{
		static char s[32];
		sprintf(s, "./_pf%d", getpid());
		FilterOutPrefix = s;
		}
	size_t n = strlen(FilterOutPrefix) + strlen(Ext) + 1;
	char *FileName = all(char, (int) n);
	strcpy(FileName, FilterOutPrefix);
	strcat(FileName, Ext);
	assert(strlen(FileName) == n - 1);
	return FileName;
	}

void PILERCR()
	{
	const char *InFileName = RequiredValueOpt("in");
	const char *OutFileName = RequiredValueOpt("out");

	const char *HitsFileName = ValueOpt("hits");
	const char *DiscardedHitsFileName = ValueOpt("discardedhits");
	const char *PilesFileName = ValueOpt("piles");
	const char *SeqsFileName = ValueOpt("seq");
	const char *FilterOutPrefix = ValueOpt("filterout");

	if (SeqsFileName != 0)
		g_fSeqs = CreateStdioFile(SeqsFileName);

	g_ShowHits = FlagOpt("showhits");
	g_LogHits = FlagOpt("loghits");
	g_LogImages = FlagOpt("logimages");
	g_LogAlns = FlagOpt("logalns");
	g_NoInfo = FlagOpt("noinfo");

	IntOpt("minhitlength", &g_MinHitLength);
	IntOpt("minrepeat", &g_MinRepeatLength);
	IntOpt("maxrepeat", &g_MaxRepeatLength);
	IntOpt("minspacer", &g_MinSpacerLength);
	IntOpt("maxspacer", &g_MaxSpacerLength);
	IntOpt("minarray", &g_MinArraySize);
	FloatOpt("minid", &g_MinId);
	FloatOpt("minrepeatratio", &g_MinRepeatRatio);
	FloatOpt("minspacerratio", &g_MinSpacerRatio);
	FloatOpt("mincons", &g_MinCons);

	bool Quiet = FlagOpt("quiet");
	if (Quiet)
		g_ShowProgress = false;

	g_DraftMinArraySize = g_MinArraySize - 2;
	if (g_DraftMinArraySize < 2)
		g_DraftMinArraySize = 2;

	Progress("Reading sequences from %s", InFileName);
	g_SeqQ = ReadMFA(InFileName, &SeqLengthQ, &ContigsQ,
	  &ContigCountQ, &ContigMapQ);
	ProgressDone();

	Progress("%d sequences, total length %d bases (%.0f Mb)",
		ContigCountQ,
		SeqLengthQ,
		SeqLengthQ/1e6);

	int g_MinHitLength = g_DraftMinRepeatLength;
	g_Diameter = 2*g_DraftMaxRepeatLength + g_DraftMaxSpacerLength;

	if (g_Diameter > SeqLengthQ)
		g_Diameter = SeqLengthQ;

	g_MaxDPHitLen = (3*g_DraftMaxRepeatLength)/2;

	FilterParams FP;
	DPParams DP;
	GetParams(SeqLengthQ, g_Diameter, g_MinHitLength, g_MinId, &FP, &DP);

	k = FP.WordSize;

	const double SeedPctId = (1.0 - (double) FP.SeedDiffs / (double) FP.SeedLength)*100.0;
	const double MemRequired = TotalMemRequired(SeqLengthQ, FP);
	const double AvgIndexList = AvgIndexListLength(SeqLengthQ, FP);

	Log("Filter parameters:\n");
	Log("   Word size       %d\n", FP.WordSize);
	Log("   Seed length     %d\n", FP.SeedLength);
	Log("   Seed diffs      %d\n", FP.SeedDiffs);
	Log("   Seed min id     %.1f%%\n", SeedPctId);
	Log("   Tube offset     %d\n", FP.TubeOffset);
	if (AvgIndexList > 2)
		Log("   Avg index list  %.1f\n", AvgIndexList);
	else
		Log("   Avg index list  %.2g\n", AvgIndexList);
	Log("DP parameters:\n");
	Log("   Min length      %d\n", DP.g_MinHitLength);
	Log("   Min id          %.0f%%\n", DP.MinId*100.0);
	Log("Estd. memory       %.0f Mb\n", MemRequired/1e6);
	Log("RAM                %.0f Mb\n", GetRAMSize()/1e6);

	int *Finger = 0;
	int *Pos = 0;

	const char *strFilterOutFileName = GetFilterOutFileName(FilterOutPrefix, ".f.tmp");
	FILE *fFilterOut = OpenStdioFile(strFilterOutFileName, FILEIO_MODE_ReadWrite);
	SetFilterOutFile(fFilterOut);

	FilterB(g_SeqQ, SeqLengthQ, FP);

	const int FilterHitCount = GetFilterHitCount();
	Progress("%d filter hits", FilterHitCount);

	FILE *fFilterOutComp = 0;
	int FilterHitCountComp = 0;
	const char *strFilterOutFileNameComp = 0;

	Progress("Read filter hits");
	FilterHit *FilterHits = ReadFilterHits(fFilterOut, FilterHitCount);
	CloseFilterOutFile();
	if (0 == FilterOutPrefix)
		remove(strFilterOutFileName);
	Progress("%d filter hits", FilterHitCount);
//	SaveFilterHits(FilterHits, FilterHitCount);

	Progress("Merge filter hits");
	int TrapCount;
	Trapezoid *Traps = MergeFilterHits(g_SeqQ, SeqLengthQ, g_SeqQ, SeqLengthQ,
	  true, FilterHits, FilterHitCount, FP, &TrapCount);
	LogTraps(Traps);
	SaveTraps(Traps);

	freemem(FilterHits);
	FilterHits = 0;

	const int SumLengths = SumTrapLengths(Traps);
	Progress("%d trapezoids, total length %d", TrapCount, SumLengths);

	if (DiscardedHitsFileName != 0)
		fDiscardedHits = OpenStdioFile(DiscardedHitsFileName,
		  FILEIO_MODE_WriteOnly);
	
	AlignTraps(g_SeqQ, SeqLengthQ, g_SeqQ, SeqLengthQ, Traps, TrapCount, false, DP);

	if (HitsFileName != 0)
		{
		FILE *fHits = OpenStdioFile(HitsFileName, FILEIO_MODE_WriteOnly);
		WriteDPHits(fHits, false);
		fclose(fHits);
		}

	Progress("Sorting hits");
	std::sort(g_Hits.begin(), g_Hits.end(), Cmp_Lo);

	if (g_LogHits)
		LogHits();

	if (fDiscardedHits != 0)
		{
		fclose(fDiscardedHits);
		fDiscardedHits = 0;
		}

	Progress("%d DP hits, %d accepted",
	  g_NotAcceptedHitCount + g_AcceptedHitCount, g_AcceptedHitCount);
	Progress("%u discarded hits", g_DiscardedHitCount);
	HitsToImages();

	Progress("Finding piles");
	FindPiles();
	ExpandPiles();
	WritePiles();
	Progress("%d piles found", (int) g_Piles.size());

	Progress("Building graph");
	EdgeList Edges;
	MakeGraph(Edges);

	Progress("Find connected components");
	CompList Comps;
	GetConnComps(Edges, Comps, g_DraftMinArraySize);
	Progress("%d connected components", (int) Comps.size());

	Progress("Find arrays");
	std::vector<ArrayData *> ADVec;
	for (CompList::const_iterator p = Comps.begin(); p != Comps.end(); ++p)
		{
		const CompData *Comp = *p;
		FindArrays(*Comp, ADVec);
		}

	g_fOut = CreateStdioFile(OutFileName);

	if (!g_NoInfo)
		Info();

	std::vector<ArrayAln *> AAVec;
	AlignArrays(ADVec, AAVec);

	const int ArrayCount = (int) AAVec.size();
	Progress("%d arrays found", (int) AAVec.size());

	Out(PILERCR_LONG_VERSION "\n");
	Out("By Robert C. Edgar\n\n");
	Out("%s: %d putative CRISPR arrays found.\n\n", InFileName, ArrayCount);

	if (ArrayCount == 0)
		return;

	OutArrays(AAVec);
	fclose(g_fOut);
	g_fOut = 0;

	if (g_fSeqs != 0)
		{
		SaveSeqs(g_fSeqs, AAVec);
		g_fSeqs = 0;
		}

	freemem(g_SeqQ);
	g_SeqQ = 0;

	free(Traps);
	Traps = 0;

	Progress("Done");
	}
