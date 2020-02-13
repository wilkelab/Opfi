#include "pilercr.h"

static bool CmpAA(const ArrayAln *AA1, const ArrayAln *AA2)
	{
	return AA1->Pos < AA2->Pos;
	}

void AlignArrays(const std::vector<ArrayData *> &ADVec,
  std::vector<ArrayAln *> &AAVec)
	{
	AAVec.clear();

// Create array alignments
	const size_t ArrayDataCount = ADVec.size();
	for (size_t Index = 0; Index < ArrayDataCount; ++Index)
		{
		ArrayData &AD = *(ADVec[Index]);
		ArrayAln &AA = *new ArrayAln;
		AlignArray(AD, AA);

		if (AA.Repeats.size() > 0 && AcceptedAA(AA))
			AAVec.push_back(&AA);
		}

// Sort by position
	std::sort(AAVec.begin(), AAVec.end(), CmpAA);

// Check for arrays where we missed some repeats, these got split
// and can be merged.
	size_t ArrayIndex = 0;
	for (;;)
		{
		size_t ArrayCount = AAVec.size();
		if (ArrayIndex + 1 >= ArrayCount)
			break;
		bool Merged = MergeAAs(*AAVec[ArrayIndex], *AAVec[ArrayIndex+1]);
		if (Merged)
			{
		// Delete entry from array and continue.
		// This allows merger of three or more arrays.
			assert(ArrayCount > 1);
			for (size_t i = ArrayIndex + 1; i < ArrayCount - 1; ++i)
				AAVec[i] = AAVec[i+1];
			AAVec.resize(ArrayCount - 1);
			}
		++ArrayIndex;
		}

// Assign ids
	size_t ArrayCount = AAVec.size();
	int Id = 1;
	for (size_t Index = 0; Index < ArrayCount; ++Index)
		{
		ArrayAln &AA = *(AAVec[Index]);
		AA.Id = Id++;
		}
	}
