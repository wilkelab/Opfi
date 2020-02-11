#include "multaln.h"

SCORE AlignProfiles(
  const ProfPos *PA, unsigned uLengthA,
  const ProfPos *PB, unsigned uLengthB,
  PWPath &Path, ProfPos **ptrPout, unsigned *ptruLengthOut)
	{
	assert(uLengthA < 100000);
	assert(uLengthB < 100000);

	SCORE Score = GlobalAlign(PA, uLengthA, PB, uLengthB, Path);

	AlignTwoProfsGivenPath(Path, PA, uLengthB, PB, uLengthB, ptrPout, ptruLengthOut);

	return Score;
	}
