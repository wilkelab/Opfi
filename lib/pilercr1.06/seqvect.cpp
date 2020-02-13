#include "multaln.h"
#include <errno.h>

const size_t MAX_FASTA_LINE = 16000;

const int BUFFER_BYTES = 16*1024;
const int CR = '\r';
const int NL = '\n';

#define ADD(c)															\
		{																\
		if (Pos >= BufferLength)										\
			{															\
			const int NewBufferLength = BufferLength + BUFFER_BYTES;	\
			char *NewBuffer	= new char[NewBufferLength];				\
			memcpy(NewBuffer, Buffer, BufferLength);					\
			delete[] Buffer;											\
			Buffer = NewBuffer;											\
			BufferLength = NewBufferLength;								\
			}															\
		Buffer[Pos++] = c;												\
		}

// Get next sequence from file.
static char *GetFastaSeq(FILE *f, unsigned *ptrSeqLength, char **ptrLabel,
  bool DeleteGaps)
	{
	unsigned BufferLength = 0;
	unsigned Pos = 0;
	char *Buffer = 0;

	int c = fgetc(f);
	if (EOF == c)
		return 0;
	if ('>' != c)
		Quit("Invalid file format, expected '>' to start FASTA label");

	for (;;)
		{
		int c = fgetc(f);
		if (EOF == c)
			Quit("End-of-file or input error in FASTA label");

	// Ignore CR (discard, do not include in label)
		if (CR == c)
			continue;

	// NL terminates label
		if (NL == c)
			break;

	// All other characters added to label
		ADD(c)
		}

// Nul-terminate label
	ADD(0)
	*ptrLabel = Buffer;

	BufferLength = 0;
	Pos = 0;
	Buffer = 0;
	int PreviousChar = NL;
	for (;;)
		{
		int c = fgetc(f);
		if (EOF == c)
			{
			if (feof(f))
				break;
			else if (ferror(f))
				Quit("Error reading FASTA file, ferror=TRUE feof=FALSE errno=%d %s",
				  errno, strerror(errno));
			else
				Quit("Error reading FASTA file, fgetc=EOF feof=FALSE ferror=FALSE errno=%d %s",
				  errno, strerror(errno));
			}

		if ('>' == c)
			{
			if (NL == PreviousChar)
				{
				ungetc(c, f);
				break;
				}
			else
				Quit("Unexpected '>' in FASTA sequence data");
			}
		else if (isspace(c))
			;
		else if (IsGapChar(c))
			{
			if (!DeleteGaps)
				ADD(c)
			}
		else if (isalpha(c))
			{
			c = toupper(c);
			ADD(c)
			}
		else if (isprint(c))
			{
			Warning("Invalid character '%c' in FASTA sequence data, ignored", c);
			continue;
			}
		else
			{
			Warning("Invalid byte hex %02x in FASTA sequence data, ignored", (unsigned char) c);
			continue;
			}
		PreviousChar = c;
		}

	if (0 == Pos)
		return GetFastaSeq(f, ptrSeqLength, ptrLabel, DeleteGaps);

	*ptrSeqLength = Pos;
	return Buffer;
	}

SeqVect::~SeqVect()
	{
	Clear();
	}

void SeqVect::Clear()
	{
	}

void SeqVect::PadToMSA(MSA &msa)
	{
	unsigned uSeqCount = Length();
	if (0 == uSeqCount)
		{
		msa.Clear();
		return;
		}

	unsigned uLongestSeqLength = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		unsigned uColCount = ptrSeq->Length();
		if (uColCount > uLongestSeqLength)
			uLongestSeqLength = uColCount;
		}
	msa.SetSize(uSeqCount, uLongestSeqLength);
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		msa.SetSeqName(uSeqIndex, ptrSeq->GetName());
		unsigned uColCount = ptrSeq->Length();
		unsigned uColIndex;
		for (uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			char c = ptrSeq->at(uColIndex);
			msa.SetChar(uSeqIndex, uColIndex, c);
			}
		while (uColIndex < uLongestSeqLength)
			msa.SetChar(uSeqIndex, uColIndex++, '.');
		}
	}

void SeqVect::Copy(const SeqVect &rhs)
	{
	clear();
	unsigned uSeqCount = rhs.Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = rhs.at(uSeqIndex);
		Seq *ptrSeqCopy = new Seq;
		ptrSeqCopy->Copy(*ptrSeq);
		push_back(ptrSeqCopy);
		}
	}

void SeqVect::StripGaps()
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->StripGaps();
		}
	}

void SeqVect::StripGapsAndWhitespace()
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->StripGapsAndWhitespace();
		}
	}

void SeqVect::ToUpper()
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->ToUpper();
		}
	}

bool SeqVect::FindName(const char *ptrName, unsigned *ptruIndex) const
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const Seq *ptrSeq = at(uSeqIndex);
		if (0 == stricmp(ptrSeq->GetName(), ptrName))
			{
			*ptruIndex = uSeqIndex;
			return true;
			}
		}
	return false;
	}

void SeqVect::AppendSeq(const Seq &s)
	{
	Seq *ptrSeqCopy = new Seq;
	ptrSeqCopy->Copy(s);
	push_back(ptrSeqCopy);
	}

void SeqVect::LogMe() const
	{
	unsigned uSeqCount = Length();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const Seq *ptrSeq = at(uSeqIndex);
		ptrSeq->LogMe();
		}
	}

const char *SeqVect::GetSeqName(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	const Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->GetName();
	}

unsigned SeqVect::GetSeqId(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	const Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->GetId();
	}

unsigned SeqVect::GetSeqIdFromName(const char *Name) const
	{
	const unsigned uSeqCount = GetSeqCount();
	for (unsigned i = 0; i < uSeqCount; ++i)
		{
		if (!strcmp(Name, GetSeqName(i)))
			return GetSeqId(i);
		}
	Quit("SeqVect::GetSeqIdFromName(%s): not found", Name);
	return 0;
	}

Seq &SeqVect::GetSeqById(unsigned uId)
	{
	const unsigned uSeqCount = GetSeqCount();
	for (unsigned i = 0; i < uSeqCount; ++i)
		{
		if (GetSeqId(i) == uId)
			return GetSeq(i);
		}
	Quit("SeqVect::GetSeqIdByUd(%d): not found", uId);
	return (Seq &) *((Seq *) 0);
	}

unsigned SeqVect::GetSeqLength(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	const Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->Length();
	}

Seq &SeqVect::GetSeq(unsigned uSeqIndex)
	{
	assert(uSeqIndex < size());
	return *at(uSeqIndex);
	}

const Seq &SeqVect::GetSeq(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < size());
	return *at(uSeqIndex);
	}

void SeqVect::SetSeqId(unsigned uSeqIndex, unsigned uId)
	{
	assert(uSeqIndex < size());
	Seq *ptrSeq = at(uSeqIndex);
	return ptrSeq->SetId(uId);
	}

void SeqVect::FromFASTAFile(FILE *f)
	{
	Clear();

	for (;;)
		{
		char *Label;
		unsigned uLength;
		char *SeqData = GetFastaSeq(f, &uLength, &Label, true);
		if (0 == SeqData)
			return;
		Seq *ptrSeq = new Seq;

		for (unsigned i = 0; i < uLength; ++i)
			{
			char c = SeqData[i];
			ptrSeq->push_back(c);
			}

		ptrSeq->SetName(Label);
		push_back(ptrSeq);

		delete[] SeqData;
		delete[] Label;
		}
	}
