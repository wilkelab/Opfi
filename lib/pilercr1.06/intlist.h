#ifndef IntList_H
#define IntList_H

#include <list>
#include <vector>
#include <set>
#include <map>

typedef std::list<unsigned> IntList;
typedef IntList::iterator IntListIter;
typedef IntList::const_iterator CIntListIter;

typedef std::vector<unsigned> IntVec;
typedef IntVec::iterator IntVecIter;
typedef IntVec::const_iterator CIntVecIter;

typedef std::set<unsigned> IntSet;
typedef IntSet::iterator IntSetIter;
typedef IntSet::const_iterator CIntSetIter;

typedef std::map<unsigned, unsigned> IntMap;
typedef IntMap::iterator IntMapIter;
typedef IntMap::const_iterator CIntMapIter;

typedef std::vector<IntVec> IntVecVec;
typedef std::vector<IntList> IntListList;

typedef std::vector<bool> BoolVec;

// Warning: sorts A and B in place:
void IntersectIntLists2(IntList &A, IntList &B, IntList &Common);

// More expensive but safer (makes tmp copies of A and B):
void IntersectIntLists(const IntList &A, const IntList &B, IntList &Common);

void LogList(const IntList &List);

#define	for_IntList(L, p)		for (IntListIter p = L.begin(); p != L.end(); ++p)
#define	for_CIntList(L, p)		for (CIntListIter p = L.begin(); p != L.end(); ++p)

#define	for_IntVec(L, p)		for (IntVecIter p = L.begin(); p != L.end(); ++p)
#define	for_CIntVec(L, p)		for (CIntVecIter p = L.begin(); p != L.end(); ++p)

#define	for_IntSet(L, p)		for (IntSetIter p = L.begin(); p != L.end(); ++p)
#define	for_CIntSet(L, p)		for (CIntSetIter p = L.begin(); p != L.end(); ++p)

#define	for_IntMap(L, p)		for (IntMapIter p = L.begin(); p != L.end(); ++p)
#define	for_CIntMap(L, p)		for (CIntMapIter p = L.begin(); p != L.end(); ++p)

#endif	// IntList_H
