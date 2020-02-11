#include "pilercr.h"
#include <math.h>
#include <algorithm>

unsigned char CompChar[256];
int CharToLetter[256];

static void InitCharToLetter()
	{
	for (int i = 0; i < 256; ++i)
		CharToLetter[i] = 4;
	CharToLetter['A'] = 0;
	CharToLetter['a'] = 0;
	CharToLetter['C'] = 1;
	CharToLetter['c'] = 1;
	CharToLetter['G'] = 2;
	CharToLetter['g'] = 2;
	CharToLetter['T'] = 3;
	CharToLetter['t'] = 3;
	}

static void InitCompChar()
	{
	for (unsigned i = 0; i < 256; ++i)
		CompChar[i] = (unsigned char) i;

	CompChar['a'] = 't';
	CompChar['c'] = 'g';
	CompChar['g'] = 'c';
	CompChar['t'] = 'a';
	CompChar['n'] = 'n';
	CompChar['A'] = 'T';
	CompChar['C'] = 'G';
	CompChar['G'] = 'C';
	CompChar['T'] = 'A';
	}

static bool InitUtils()
	{
	InitCompChar();
	InitCharToLetter();
	return true;
	}

static bool UtilsInitialized = InitUtils();

char *strsave(const char *s)
	{
	if (0 == s)
		return 0;
	char *ptrCopy = strdup(s);
	if (0 == ptrCopy)
		Quit("Out of memory");
	return ptrCopy;
	}

static char *StdioStrMode(FILEIO_MODE Mode)
	{
	switch (Mode)
		{
	case FILEIO_MODE_ReadOnly:
		return FILIO_STRMODE_ReadOnly;
	case FILEIO_MODE_WriteOnly:
		return FILIO_STRMODE_WriteOnly;
	case FILEIO_MODE_ReadWrite:
		return FILIO_STRMODE_ReadWrite;
		}
	Quit("StdioStrMode: Invalid mode");
	return "r";
	}

FILE *OpenStdioFile(const char *FileName, FILEIO_MODE Mode)
	{
	char *strMode = StdioStrMode(Mode);
	FILE *f = fopen(FileName, strMode);
	if (0 == f)
		Quit("Cannot open %s, %s [%d]", FileName, strerror(errno), errno);
	return f;
	}

FILE *CreateStdioFile(const char *FileName)
	{
	FILE *f = fopen(FileName, "wb+");
	if (0 == f)
		Quit("Cannot open %s, %s [%d]", FileName, strerror(errno), errno);
	return f;
	}

int GetFileSize(FILE *f)
	{
	long CurrPos = ftell(f);
	if (CurrPos < 0)
		Quit("FileSize: ftell<0 (CurrPos), errno=%d", errno);

	int Ok = fseek(f, 0, SEEK_END);
	if (Ok != 0)
		Quit("FileSize fseek(END) != 0 errno=%d", errno);

	long Size = ftell(f);
	if (Size < 0)
		Quit("FileSize: ftell<0 (size), errno=%d", errno);

	Ok = fseek(f, CurrPos, SEEK_SET);
	if (Ok != 0)
		Quit("FileSize fseek(restore curr pos) != 0 errno=%d", errno);

	long NewPos = ftell(f);
	if (CurrPos < 0)
		Quit("FileSize: ftell=%ld != CurrPos=%ld", CurrPos, NewPos);

	return (int) Size;
	}

// 4^n
int pow4(int n)
	{
	assert(n >= 0 && n < 16);
	return (1 << (2*n));
	}

// 4^d, but much less likely to under or overflow
double pow4d(int n)
	{
	return pow(4.0, n);
	}

double log4(double x)
	{
	static double LOG4 = log(4.0);
	return log(x)/LOG4;
	}

// Ukonnen's Lemma
int MinWordsPerFilterHit(int HitLength, int WordLength, int MaxErrors)
	{
	return HitLength + 1 - WordLength*(MaxErrors + 1);
	}

int StringToCode(const char s[], int len)
	{
	int code = 0;
	for (int i = 0; i < len; ++i)
		{
		char c = s[i];
		switch (c)
			{
		case 'a': case 'A':
			code <<= 2;
			code |= 0x00;
			break;
		case 'c': case 'C':
			code <<= 2;
			code |= 0x01;
			break;
		case 'g': case 'G':
			code <<= 2;
			code |= 0x02;
			break;
		case 't': case 'T':
			code <<= 2;
			code |= 0x03;
			break;
		default:
			return -1;
			}
		}
	return code;
	}

char *CodeToString(int code, int len)
	{
	static char Str[100];
	static char Symbol[] = { 'A', 'C', 'G', 'T' };
	int i;

	assert(len < sizeof(Str));
	Str[len] = 0;
	for (i = len-1; i >= 0; --i)
		{
		Str[i] = Symbol[code & 0x3];
		code >>= 2;
		}
	return Str;
	}

void *ckalloc(int size, const char *where)
{ void *p;
  p = malloc(size);
  if (p == NULL) Quit("Out of memory (ckalloc)");
  return (p);
}

void *ckrealloc(void *p, int size, const char *where)
{ p = realloc(p,size);
  if (p == NULL) Quit("Out of memory (ckrealloc)");
  return (p);
}

int RevCompKmer(int Kmer)
	{
//	Log("Kmer=%s", CodeToString(Kmer, k));
	int RCKmer = 0;
	for (int i = 0; i < k; ++i)
		{
		int x = Kmer & 3;
		int cx = 3 - x;
		int shift = 2*(k - i - 1);
		RCKmer |= (cx << shift);
		Kmer >>= 2;
		}
//	Log(" Comp=%s\n", CodeToString(RCKmer, k));
	return RCKmer;
	}

int GetHitLength(const DPHit &Hit)
	{
	assert(Hit.aHi >= Hit.aLo);
// Plus only
	assert(Hit.bHi >= Hit.bLo);
	return (Hit.aHi - Hit.aLo + 1 + Hit.bHi - Hit.bLo + 1)/2;
	}

int GetSpacerLength(const DPHit &Hit)
	{
	assert(Hit.aHi >= Hit.aLo);
// Plus only
	assert(Hit.bHi >= Hit.bLo);
	return Hit.bLo - Hit.aHi;
	}

int GetOverlap(unsigned lo1, unsigned hi1, unsigned lo2, unsigned hi2)
	{
    int maxlo = std::max(lo1, lo2);
    int minhi = std::min(hi1, hi2);

    if (maxlo > minhi)
        return 0;

    return minhi - maxlo + 1;
	}

double GetRatio(int x, int y)
	{
	if (x == 0 && y == 0)
		return 0;
	if (x < y)
		return (double) x / (double) y;
	else
		return (double) y / (double) x;
	}

int GetPileLength(const PileData &Pile)
	{
	return Pile.Hi - Pile.Lo + 1;
	}

int GetSpacerLength(const PileData &Pile1, const PileData &Pile2)
	{
	assert(Pile2.Lo >= Pile1.Hi);
	return Pile2.Lo - Pile1.Hi;
	}

void DeletePile(int PileIndex)
	{
	assert(PileIndex >= 0 && PileIndex < g_PileCount);
	g_Piles[PileIndex].Deleted = true;
	}

const ImageData &GetImage(int ImageIndex)
	{
	assert(ImageIndex >= 0 && ImageIndex < g_ImageCount);
	return g_Images[ImageIndex];
	}

const DPHit &GetHit(int HitIndex)
	{
	assert(HitIndex >= 0 && HitIndex < g_HitCount);
	return g_Hits[HitIndex];
	}

const PileData &GetPile(int PileIndex)
	{
	assert(PileIndex >= 0 && PileIndex < g_PileCount);
	return g_Piles[PileIndex];
	}

PileData &GetModifiablePile(int PileIndex)
	{
	assert(PileIndex >= 0 && PileIndex < g_PileCount);
	return g_Piles[PileIndex];
	}

DPHit &GetModifiableHit(int HitIndex)
	{
	assert(HitIndex >= 0 && HitIndex < g_HitCount);
	return g_Hits[HitIndex];
	}

int GetHitLength(int ImageIndex)
	{
	const ImageData &Image = GetImage(ImageIndex);
	const DPHit &Hit = GetHit(Image.HitIndex);
	return GetHitLength(Hit);
	}

int GetSpacerLengthFromImage(int ImageIndex)
	{
	const ImageData &Image = GetImage(ImageIndex);
	const DPHit &Hit = GetHit(Image.HitIndex);
	return GetSpacerLength(Hit);
	}

int GetOverlap(const ImageData &Image1, const ImageData &Image2)
	{
	return GetOverlap(Image1.Lo, Image1.Hi, Image2.Lo, Image2.Hi);
	}

void LogHit(const DPHit &Hit)
	{
	Log("pos1=%d pos2=%d hitlen=%d space=%d",
	  Hit.aLo, Hit.bLo, GetHitLength(Hit), GetSpacerLength(Hit));
	}

void LogHit(int HitIndex)
	{
	const DPHit &Hit = GetHit(HitIndex);
	Log("Hit%d ", HitIndex);
	LogHit(Hit);
	if (IsDeletedHit(HitIndex))
		Log(" DEL");
	}

bool IsDeletedHit(int HitIndex)
	{
	assert(HitIndex >= 0 && HitIndex < g_HitCount);
	const DPHit &Hit = g_Hits[HitIndex];
	if (Hit.Deleted)
		return true;

	if (g_PileCount > 0)
		{
		int PileIndex = g_HitIndexToPileIndexA[HitIndex];
		const PileData &Pile = GetPile(PileIndex);
		if (Pile.Deleted)
			return true;
		}

	return false;
	}

void DeleteHit(int HitIndex)
	{
	assert(HitIndex >= 0 && HitIndex < g_HitCount);
	DPHit &Hit = g_Hits[HitIndex];
	Hit.Deleted = true;
	int PileIndex = g_HitIndexToPileIndexA[HitIndex];
	PileData &Pile = GetModifiablePile(PileIndex);
	Pile.Deleted = true;
	}
