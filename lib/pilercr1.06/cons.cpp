#include "pilercr.h"
#include "multaln.h"

#define TRACE	0

static const char Letter[5] = { 'A', 'C', 'G', 'T', '-'};

static void GetCounts(const MSA &Aln, int ColIndex, int Counts[5])
	{
	memset(Counts, 0, 5*sizeof(unsigned));
	const int SeqCount = Aln.GetSeqCount();
	for (int SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		char c = Aln.GetChar(SeqIndex, ColIndex);
		c = toupper(c);
		switch (c)
			{
		case 'a':
		case 'A':
			++(Counts[0]);
			break;
		case 'c':
		case 'C':
			++(Counts[1]);
			break;
		case 'g':
		case 'G':
			++(Counts[2]);
			break;
		case 't':
		case 'T':
			++(Counts[3]);
			break;
		case '-':
			++(Counts[4]);
			break;
			}
		}
	}

static double GetCons(const MSA &Aln, int Col, char *ptrc)
	{
	int Counts[5];
	GetCounts(Aln, Col, Counts);
	int MaxCount = -1;
	char c = '?';
	for (int i = 0; i < 5; ++i)
		{
		if (Counts[i] > MaxCount)
			{
			MaxCount = Counts[i];
			c = Letter[i];
			}
		}

	*ptrc = c;
	return (double) MaxCount / (double) Aln.GetSeqCount();;
	}

#if	TRACE
static void LogCol(const MSA &Aln, int Col)
	{
	int SeqCount = Aln.GetSeqCount();
	for (int i = 0; i < SeqCount; ++i)
		Log("%c", Aln.GetChar(i, Col));
	}
#endif

// Column is conserved if the two preceding and two following columns
// are conserved.
static void GetSmoothedCons(const MSA &Aln, double MinCons, BoolVec &ColCons)
	{
	int ColCount = (int) Aln.GetColCount();
	ColCons.resize(ColCount, false);

	for (int Col = 0; Col < ColCount; ++Col)
		{
		char c;
		ColCons[Col] = (GetCons(Aln, Col, &c) >= MinCons);
		}

	for (int Col = 0; Col < ColCount; ++Col)
		{
		if (ColCons[Col])
			continue;
		if (!ColCons[Col] && 
		  (Col > 0 && ColCons[Col-1]) &&
		  (Col > 1 && ColCons[Col-2]) &&
		  (Col < ColCount - 1 && ColCons[Col+1]) &&
		  (Col < ColCount - 2 && ColCons[Col+2]))
			ColCons[Col] = true;
		}
	}

static void FindStartEnd(const MSA &Aln, double MinCons, int *ptrStart, int *ptrEnd)
	{
	BoolVec ColCons;
	GetSmoothedCons(Aln, MinCons, ColCons);

	int BestStart = 0;
	int BestEnd = 0;
	int Length = 0;
	int Start = 0;

	int ColCount = (int) Aln.GetColCount();
	for (int Col = 0; Col <= ColCount; ++Col)
		{
#if	TRACE
		if (Col != ColCount)
			LogCol(Aln, Col);
#endif

		bool ConsOk = (Col != ColCount && ColCons[Col]);
		bool ContinueBlock = (Col != ColCount && ConsOk);
#if	TRACE
		if (ContinueBlock)
			Log("  cont\n");
		else
			Log("  ** END BLOCK **\n");
#endif

		if (ContinueBlock)
			{
			if (Start == -1)
				{
				Start = Col;
				Length = 1;
				}
			else
				++Length;
			}
		else
			{
			if (Length > BestEnd - BestStart + 1)
				{
				BestStart = Start;
				BestEnd = Col - 1;
				if (BestEnd - BestStart + 1 != Length)
					Quit("BestEnd");
				}
			Length = -1;
			Start = -1;
			}
		}
	*ptrStart = BestStart;
	*ptrEnd = BestEnd;
	}

static char ConsNextChar(const MSA &Aln, int EndCol)
	{
	const unsigned SeqCount = Aln.GetSeqCount();
	const unsigned ColCount = Aln.GetColCount();
	char Cons = '-';
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		for (unsigned Col = EndCol + 1; Col < ColCount; ++Col)
			{
			char c = Aln.GetChar(SeqIndex, Col);
			if (c == '-')
				continue;
			else
				{
				if (Col == EndCol + 1)
					Cons = c;
				else if (Cons != c)
					return '-';
				break;
				}
			}
		}
	return Cons;
	}

void GetConsSeq(const MSA &Aln, double MinCons, int *ptrStartCol,
  int *ptrEndCol, Seq &ConsSeq)
	{
	int SeqCount = Aln.GetSeqCount();
	if (SeqCount == 0)
		{
		ConsSeq.clear();
		*ptrStartCol = -1;
		*ptrEndCol = -1;
		return;
		}

	int ColCount = Aln.GetColCount();
	int StartCol = 0;
	int EndCol = ColCount - 1;

	FindStartEnd(Aln, MinCons, &StartCol, &EndCol);

	for (int Col = StartCol; Col <= EndCol; ++Col)
		{
		char c;
		GetCons(Aln, Col, &c);
		if (c != '-')
			ConsSeq.push_back(c);
		}

// Eggregious hack to compensate for missing the end
// of the concensus sequence sometimes. Check the
// following letter to see if it is conserved.
	char Nextc = ConsNextChar(Aln, EndCol);
	if (Nextc != '-')
		ConsSeq.push_back(Nextc);

	*ptrStartCol = StartCol;
	*ptrEndCol = EndCol;
	}

char GetConsSymbol(const MSA &Aln, unsigned Col)
	{
	int Counts[5];
	GetCounts(Aln, (int) Col, Counts);
	const int SeqCount = (int) Aln.GetSeqCount();
	int BigCount = 0;
	int BigLetter = 0;
	for (int i = 0; i < 5; ++i)
		{
		if (Counts[i] > BigCount)
			{
			BigCount = Counts[i];
			BigLetter = i;
			}
		}
	if (BigLetter == 4)
		return ' ';
	if (BigCount == SeqCount)
		return '*';
	if (((double) BigCount / (double) SeqCount) >= g_ColonThresh)
		return ':';
	if (Counts[BigLetter] == SeqCount - 1 && SeqCount > 5)
		return ':';
	return ' ';
	}

void GetConsSymbols(const MSA &Aln, Seq &ConsSymbols)
	{
	const unsigned ColCount = Aln.GetColCount();
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = GetConsSymbol(Aln, Col);
		ConsSymbols.push_back(c);
		}
	}
