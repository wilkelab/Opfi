#include "multaln.h"

const size_t MAX_FASTA_LINE = 16000;

static char GetWildcardChar()
	{
	return 'N';
	}

void InvalidLetterWarning(char c, char w)
	{
	if (isprint(c))
		Warning("Invalid char '%c' in sequence, replaced by '%c'", c, w);
	else
		Warning("Invalid char 0x02x in sequence, replaced by '%c'", c, w);
	}

void Seq::SetName(const char *ptrName)
	{
	delete[] m_ptrName;
	size_t n = strlen(ptrName) + 1;
	m_ptrName = new char[n];
	strcpy(m_ptrName, ptrName);
	}

void Seq::ExtractUngapped(MSA &msa) const
	{
	msa.Clear();
	unsigned uColCount = Length();
	msa.SetSize(1, 1);
	unsigned uUngappedPos = 0;
	for (unsigned n = 0; n < uColCount; ++n)
		{
		char c = at(n);
		if (!IsGapChar(c))
			msa.SetChar(0, uUngappedPos++, c);
		}
	msa.SetSeqName(0, m_ptrName);
	}

void Seq::Copy(const Seq &rhs)
	{
	clear();
	const unsigned uLength = rhs.Length();
	for (unsigned uColIndex = 0; uColIndex < uLength; ++uColIndex)
		push_back(rhs.at(uColIndex));
	const char *ptrName = rhs.GetName();
	if (ptrName == 0)
		Quit("Seq::Copy: Name=NULL");
	size_t n = strlen(ptrName) + 1;
	m_ptrName = new char[n];
	strcpy(m_ptrName, ptrName);
	SetId(rhs.GetId());
	}

void Seq::CopyReversed(const Seq &rhs)
	{
	clear();
	const unsigned uLength = rhs.Length();
	const unsigned uBase = rhs.Length() - 1;
	for (unsigned uColIndex = 0; uColIndex < uLength; ++uColIndex)
		push_back(rhs.at(uBase - uColIndex));
	const char *ptrName = rhs.GetName();
	size_t n = strlen(ptrName) + 1;
	m_ptrName = new char[n];
	strcpy(m_ptrName, ptrName);
	}

void Seq::StripGaps()
	{
	for (CharVect::iterator p = begin(); p != end(); )
		{
		char c = *p;
		if (IsGapChar(c))
			erase(p);
		else
			++p;
		}
	}

void Seq::StripGapsAndWhitespace()
	{
	for (CharVect::iterator p = begin(); p != end(); )
		{
		char c = *p;
		if (isspace(c) || IsGapChar(c))
			erase(p);
		else
			++p;
		}
	}

void Seq::ToUpper()
	{
	for (CharVect::iterator p = begin(); p != end(); ++p)
		{
		char c = *p;
		if (islower(c))
			*p = toupper(c);
		}
	}

unsigned Seq::GetLetter(unsigned uIndex) const
	{
	assert(uIndex < Length());
	char c = operator[](uIndex);
	return CharToLetter(c);
	}

bool Seq::EqIgnoreCase(const Seq &s) const
	{
	const unsigned n = Length();
	if (n != s.Length())
		return false;
	for (unsigned i = 0; i < n; ++i)
		{
		const char c1 = at(i);
		const char c2 = s.at(i);
		if (IsGapChar(c1))
			{
			if (!IsGapChar(c2))
				return false;
			}
		else
			{
			if (toupper(c1) != toupper(c2))
				return false;
			}
		}
	return true;
	}

bool Seq::Eq(const Seq &s) const
	{
	const unsigned n = Length();
	if (n != s.Length())
		return false;
	for (unsigned i = 0; i < n; ++i)
		{
		const char c1 = at(i);
		const char c2 = s.at(i);
		if (c1 != c2)
			return false;
		}
	return true;
	}

bool Seq::EqIgnoreCaseAndGaps(const Seq &s) const
	{
	const unsigned uThisLength = Length();
	const unsigned uOtherLength = s.Length();
	
	unsigned uThisPos = 0;
	unsigned uOtherPos = 0;

	int cThis;
	int cOther;
	for (;;)
		{
		if (uThisPos == uThisLength && uOtherPos == uOtherLength)
			break;

	// Set cThis to next non-gap character in this string
	// or -1 if end-of-string.
		for (;;)
			{
			if (uThisPos == uThisLength)
				{
				cThis = -1;
				break;
				}
			else
				{
				cThis = at(uThisPos);
				++uThisPos;
				if (!IsGapChar(cThis))
					{
					cThis = toupper(cThis);
					break;
					}
				}
			}

	// Set cOther to next non-gap character in s
	// or -1 if end-of-string.
		for (;;)
			{
			if (uOtherPos == uOtherLength)
				{
				cOther = -1;
				break;
				}
			else
				{
				cOther = s.at(uOtherPos);
				++uOtherPos;
				if (!IsGapChar(cOther))
					{
					cOther = toupper(cOther);
					break;
					}
				}
			}

	// Compare characters are corresponding ungapped position
		if (cThis != cOther)
			return false;
		}
	return true;
	}

unsigned Seq::GetUngappedLength() const
	{
	unsigned uUngappedLength = 0;
	for (CharVect::const_iterator p = begin(); p != end(); ++p)
		{
		char c = *p;
		if (!IsGapChar(c))
			++uUngappedLength;
		}
	return uUngappedLength;
	}

void Seq::LogMe() const
	{
	Log(">%s\n", m_ptrName);
	const unsigned n = Length();
	for (unsigned i = 0; i < n; ++i)
		Log("%c", at(i));
	Log("\n");
	}

void Seq::LogMeSeqOnly() const
	{
	const unsigned n = Length();
	for (unsigned i = 0; i < n; ++i)
		Log("%c", at(i));
	}

void Seq::FromString(const char *pstrSeq, const char *pstrName)
	{
	clear();
	const unsigned uLength = (unsigned) strlen(pstrSeq);
	for (unsigned uColIndex = 0; uColIndex < uLength; ++uColIndex)
		push_back(pstrSeq[uColIndex]);
	size_t n = strlen(pstrName) + 1;
	m_ptrName = new char[n];
	strcpy(m_ptrName, pstrName);
	}

bool Seq::HasGap() const
	{
	for (CharVect::const_iterator p = begin(); p != end(); ++p)
		{
		char c = *p;
		if (IsGapChar(c))
			return true;
		}
	return false;
	}

void Seq::FixAlpha()
	{
	for (CharVect::iterator p = begin(); p != end(); ++p)
		{
		char c = *p;
		if (!IsResidueChar(c))
			{
			char w = GetWildcardChar();
			InvalidLetterWarning(c, w);
			*p = w;
			}
		}
	}

void Seq::ToString(char String[], unsigned Bytes) const
	{
	const unsigned N = (unsigned) size();
	if (Bytes <= N)
		Quit("Seq::ToString, buffer too small");
	for (unsigned i = 0; i < N; ++i)
		String[i] = GetChar(i);
	String[N] = 0;
	}

void Seq::RevComp()
	{
	extern unsigned char CompMap[256];
	iterator s = begin();
	iterator t = end() - 1;
	while (s < t)
		{	
		unsigned char c = (unsigned char) *s;
		*s++ = CompMap[*t];
		*t-- = CompMap[c];
		}
	if (s == t)
		*s = CompMap[*s];
	}
