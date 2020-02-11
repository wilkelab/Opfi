#include "pilercr.h"
#include "sw.h"

static double CmpSeqsFwd(const Seq &s1, const Seq &s2)
	{
	unsigned Start1;
	unsigned Start2;
	std::string Path;
	score_t Score = SW(s1, s2, &Start1, &Start2, Path);
	if (Score <= 0)
		return 0;
	return GetFractId(s1, s1, Start1, Start2, Path);
	}

static double CmpSeqsBwd(const Seq &s1, const Seq &s2)
	{
	unsigned Start1;
	unsigned Start2;
	std::string Path;

	Seq s1rev;
	s1rev.Copy(s1);
	s1rev.RevComp();
	score_t Score = SW(s1rev, s2, &Start1, &Start2, Path);
	if (Score <= 0)
		return 0;
	return GetFractId(s1rev, s1, Start1, Start2, Path);
	}

static double CmpSeqs(const Seq &s1, const Seq &s2)
	{
	double IdFwd = CmpSeqsFwd(s1, s2);
	if (IdFwd > 0.899)
		return IdFwd;
	return CmpSeqsBwd(s1, s2);
	}

void SaveSeqs(FILE *f, const std::vector<ArrayAln *> AAVec)
	{
	bool TrimSeqs = FlagOpt("trimseqs");
	const size_t ArrayCount = AAVec.size();
	for (size_t ArrayIndex = 0; ArrayIndex < ArrayCount; ++ArrayIndex)
		{
		const ArrayAln &AA = *(AAVec[ArrayIndex]);
		const Seq &s = AA.ConsSeq;
		if (TrimSeqs)
			{
			bool Discard = false;
			for (size_t ArrayIndex2 = 0; ArrayIndex2 < ArrayIndex; ++ArrayIndex2)
				{
				const ArrayAln &AA2 = *(AAVec[ArrayIndex2]);
				const Seq &s2 = AA2.ConsSeq;
				if (CmpSeqs(s, s2) > 0.899)
					{
					Discard = true;
					break;
					}
				}
			if (Discard)
				continue;
			}

		int ArrayLength = GetArrayLength(AA);
		char *Label;
		int LocalPos = GlobalToLocal(AA.Pos, &Label);
		fprintf(f, ">%s[Array%d;Pos=%d]", Label, AA.Id, LocalPos + 1);
		unsigned L = s.Length();
		for (unsigned i = 0; i < L; ++i)
			{
			if (i%60 == 0)
				fprintf(f, "\n");
			fprintf(f, "%c", s[i]);
			}
		fprintf(f, "\n");
		}
	}
