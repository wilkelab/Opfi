#include "multaln.h"

const unsigned MAX_CHAR = 256;

unsigned g_CharToLetter[MAX_CHAR];
unsigned g_CharToLetterEx[MAX_CHAR];

char g_LetterToChar[MAX_ALPHA];
char g_LetterExToChar[MAX_ALPHA_EX];

char g_UnalignChar[MAX_CHAR];
char g_AlignChar[MAX_CHAR];

bool g_IsWildcardChar[MAX_CHAR];
bool g_IsResidueChar[MAX_CHAR];

#define Res(c, Letter)												\
	{																\
	const unsigned char Upper = (unsigned char) toupper(c);			\
	const unsigned char Lower = (unsigned char) tolower(c);			\
	g_CharToLetter[Upper] = Letter;									\
	g_CharToLetter[Lower] = Letter;									\
	g_CharToLetterEx[Upper] = Letter;								\
	g_CharToLetterEx[Lower] = Letter;								\
	g_LetterToChar[Letter] = Upper;									\
	g_LetterExToChar[Letter] = Upper;								\
	g_IsResidueChar[Upper] = true;									\
	g_IsResidueChar[Lower] = true;									\
	g_AlignChar[Upper] = Upper;										\
	g_AlignChar[Lower] = Upper;										\
	g_UnalignChar[Upper] = Lower;									\
	g_UnalignChar[Lower] = Lower;									\
	}

#define Wild(c, Letter)												\
	{																\
	const unsigned char Upper = (unsigned char) toupper(c);			\
	const unsigned char Lower = (unsigned char) tolower(c);			\
	g_CharToLetterEx[Upper] = Letter;								\
	g_CharToLetterEx[Lower] = Letter;								\
	g_LetterExToChar[Letter] = Upper;								\
	g_IsResidueChar[Upper] = true;									\
	g_IsResidueChar[Lower] = true;									\
	g_AlignChar[Upper] = Upper;										\
	g_AlignChar[Lower] = Upper;										\
	g_UnalignChar[Upper] = Lower;									\
	g_UnalignChar[Lower] = Lower;									\
	g_IsWildcardChar[Lower] = true;									\
	g_IsWildcardChar[Upper] = true;									\
	}

static void SetGapChar(char c)
	{
	unsigned char u = (unsigned char) c;

	g_CharToLetterEx[u] = NX_GAP;
	g_LetterExToChar[NX_GAP] = (char) u;
	g_AlignChar[u] = u;
	g_UnalignChar[u] = u;
	}

static void InitArrays()
	{
	memset(g_CharToLetter, 0xff, sizeof(g_CharToLetter));
	memset(g_CharToLetterEx, 0xff, sizeof(g_CharToLetterEx));

	memset(g_LetterToChar, '?', sizeof(g_LetterToChar));
	memset(g_LetterExToChar, '?', sizeof(g_LetterExToChar));

	memset(g_AlignChar, '?', sizeof(g_UnalignChar));
	memset(g_UnalignChar, '?', sizeof(g_UnalignChar));

	memset(g_IsWildcardChar, 0, sizeof(g_IsWildcardChar));
	}

static bool InitAlpha()
	{
	InitArrays();

	SetGapChar('.');
	SetGapChar('-');

	Res('A', NX_A)
	Res('C', NX_C)
	Res('G', NX_G)
	Res('U', NX_T)
	Res('T', NX_T)
	Wild('N', NX_N)

	return true;
	}

static bool Initialized = InitAlpha();
