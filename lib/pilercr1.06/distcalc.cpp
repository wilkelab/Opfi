#include "multaln.h"

void DistCalcDF::Init(const DistFunc &DF)
	{
	m_ptrDF = &DF;
	}

void DistCalcDF::CalcDistRange(unsigned i, dist_t Dist[]) const
	{
	for (unsigned j = 0; j < i; ++j)
		Dist[j] = m_ptrDF->GetDist(i, j);
	}

unsigned DistCalcDF::GetCount() const
	{
	return m_ptrDF->GetCount();
	}

unsigned DistCalcDF::GetId(unsigned i) const
	{
	return m_ptrDF->GetId(i);
	}

const char *DistCalcDF::GetName(unsigned i) const
	{
	return m_ptrDF->GetName(i);
	}
