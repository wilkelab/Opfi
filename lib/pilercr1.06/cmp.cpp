#include "pilercr.h"
#include "sw.h"

#define TRACE	0

static double MinId = 0.899;

static void LogAln(const Seq &RefSeq, const Seq &TestSeq_, unsigned StartRef,
  unsigned StartTest, std::string &Path, bool Rev)
	{
	Seq TestSeq;
	TestSeq.Copy(TestSeq_);
	if (Rev)
		TestSeq.RevComp();

	unsigned AlnLength = (unsigned) Path.size();
	unsigned RefLength = RefSeq.Length();
	unsigned TestLength = TestSeq.Length();

// Ref seq
	if (StartTest > StartRef)
		for (unsigned i = StartRef; i < StartTest; ++i)
			Log(".");

	for (unsigned i = 0; i < StartRef; ++i)
		Log("%c", tolower(RefSeq[i]));

	unsigned RefPos = StartRef;
	for (unsigned i = 0; i < AlnLength; ++i)
		{
		switch (Path[i])
			{
		case 'M':
		case 'D':
			Log("%c", toupper(RefSeq[RefPos++]));
			break;
		case 'I':
			Log("-");
			break;
			}
		}

	while (RefPos < RefLength)
		Log("%c", tolower(RefSeq[RefPos++]));

	Log("\n");

// Test seq
	if (StartRef > StartTest)
		for (unsigned i = StartTest; i < StartRef; ++i)
			Log(".");

	for (unsigned i = 0; i < StartTest; ++i)
		Log("%c", tolower(TestSeq[i]));

	unsigned TestPos = StartTest;
	for (unsigned i = 0; i < AlnLength; ++i)
		{
		switch (Path[i])
			{
		case 'M':
		case 'I':
			Log("%c", toupper(TestSeq[TestPos++]));
			break;
		case 'D':
			Log("-");
			break;
			}
		}

	while (TestPos < TestLength)
		Log("%c", tolower(TestSeq[TestPos++]));

	Log("\n");
	}

double GetFractId(const Seq &A, const Seq &B, unsigned StartA,
  unsigned StartB, std::string &Path)
	{
	unsigned APos = StartA;
	unsigned BPos = StartB;
	unsigned Same = 0;
	unsigned AlnLength = (unsigned) Path.size();
	for (unsigned i = 0; i < AlnLength; ++i)
		{
		switch (Path[i])
			{
		case 'M':
			if (toupper(A[APos]) == toupper(B[BPos]))
				++Same;
			++APos;
			++BPos;
			break;
		case 'D':
			++APos;
			break;
		case 'I':
			++BPos;
			break;
		default:
			assert(false);
			}
		}
	unsigned LA = A.Length();
	unsigned LB = B.Length();
	if (LA == 0 || LB == 0)
		return 0.0;
	unsigned MaxLen = std::max(LA, LB);
	return (double) Same / (double) MaxLen;
	}

void Cmp()
	{
	const char *TestFileName = RequiredValueOpt("cmp");
	const char *RefFileName = RequiredValueOpt("ref");

	FILE *fTest = OpenStdioFile(TestFileName);
	FILE *fRef = OpenStdioFile(RefFileName);

	SeqVect TestSeqs;
	SeqVect RefSeqs;

	TestSeqs.FromFASTAFile(fTest);
	RefSeqs.FromFASTAFile(fRef);

	fclose(fTest);
	fclose(fRef);

	fTest = 0;
	fRef = 0;

	int TestSeqCount = (int) TestSeqs.Length();
	int RefSeqCount = (int) RefSeqs.Length();

	for (int RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
		{
		const Seq &RefSeq = *(RefSeqs[RefSeqIndex]);

		int BestTestSeqIndex = -1;
		bool Exact = false;
		bool Rev = false;
		double BestId = 0.0;
		unsigned BestStartTest;
		unsigned BestStartRef;
		std::string BestPath;
		for (int TestSeqIndex = 0; TestSeqIndex < TestSeqCount; ++TestSeqIndex)
			{
			const Seq &TestSeq = *(TestSeqs[TestSeqIndex]);
			unsigned StartTest;
			unsigned StartRef;
			std::string Path;
			score_t Score = SW(RefSeq, TestSeq, &StartRef, &StartTest, Path);
#if	TRACE
			Log("SW(%s,%s)=%.1f\n", TestSeq.GetName(), RefSeq.GetName(), Score);
			Log("Test=");
			for (unsigned i = 0; i < TestSeq.Length(); ++i)
				Log("%c", TestSeq[i]);
			Log("\n");
			Log("Ref =");
			for (unsigned i = 0; i < RefSeq.Length(); ++i)
				Log("%c", RefSeq[i]);
			Log("\n");
			if (Score > 0)
				{
				Log("Aln length = %d\n", (int) Path.size());
				LogAln(RefSeq, TestSeq, StartRef, StartTest, Path, false);
				}
#endif
			if (Score > 0)
				{
				double Id = GetFractId(TestSeq, RefSeq, StartTest, StartRef, Path);
#if	TRACE
				Log("Id = %.3f\n", Id);
#endif
				if (Id > 0.999)
					{
					Exact = true;
					Log("Exact: ref=%s, test=%s\n", RefSeq.GetName(), TestSeq.GetName());
					break;
					}
				if (Id >= MinId)
					{
					BestTestSeqIndex = TestSeqIndex;
					BestStartTest = StartTest;
					BestStartRef = StartRef;
					BestPath = Path;
					BestId = Id;
					Rev = false;
					}
				}

			Seq TestSeqRC;
			TestSeqRC.Copy(TestSeq);
			TestSeqRC.RevComp();

			Score = SW(RefSeq, TestSeqRC, &StartRef, &StartTest, Path);
#if	TRACE
			Log("SW(rev %s,%s)=%.1f\n", TestSeq.GetName(), RefSeq.GetName(), Score);
			Log("Test=");
			for (unsigned i = 0; i < TestSeqRC.Length(); ++i)
				Log("%c", TestSeqRC[i]);
			Log("\n");
			Log("Ref =");
			for (unsigned i = 0; i < RefSeq.Length(); ++i)
				Log("%c", RefSeq[i]);
			Log("\n");
			if (Score > 0)
				{
				Log("Aln length = %d\n", (int) Path.size());
				LogAln(RefSeq, TestSeq, StartRef, StartTest, Path, true);
				}
#endif
			if (Score > 0)
				{
				double Id = GetFractId(TestSeqRC, RefSeq, StartTest, StartRef, Path);
#if	TRACE
				Log("Id = %.3f\n", Id);
#endif
				if (Id > 0.999)
					{
					Exact = true;
					Log("Exact: rev ref=%s, test=%s\n", RefSeq.GetName(), TestSeq.GetName());
					break;
					}
				if (Id >= MinId)
					{
					BestTestSeqIndex = TestSeqIndex;
					BestId = Id;
					BestStartTest = StartTest;
					BestStartRef = StartRef;
					BestPath = Path;
					Rev = true;
					}
				}
			}
		if (Exact)
			continue;
		if (BestTestSeqIndex == -1)
			Log("None: best id %.1f%% ref=%s\n", BestId*100.0, RefSeq.GetName());
		else
			{
			const Seq &TestSeq = *(TestSeqs[BestTestSeqIndex]);
			Log("Close: %.1f%%: rev=%c ref=%s, test=%s\n", BestId*100.0, Rev ? 'T' : 'F',
			  RefSeq.GetName(), TestSeq.GetName());
			LogAln(RefSeq, TestSeq, BestStartRef, BestStartTest, BestPath, Rev);
			}
		}
	}
