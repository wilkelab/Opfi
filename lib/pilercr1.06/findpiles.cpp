#include "pilercr.h"

#define	TRACE	0

void LogImage(const ImageData &Image)
	{
	Log("Lo=%d Hi=%d Len=%d Hit=%d %c",
	  Image.Lo, Image.Hi, Image.Hi - Image.Lo + 1, Image.HitIndex,
	  Image.IsA ? 'A' : 'B');
	}

// Requires sorted hits as input
void FindPiles()
	{
	g_HitIndexToPileIndexA.clear();
	g_HitIndexToPileIndexB.clear();
	g_Piles.clear();

	if (g_ImageCount == 0)
		return;

// Just for speed. g_ImageCount/4 is arbitrary guess.
	g_Piles.reserve(g_ImageCount/4);

// Maps aLo-aHi segment of hit to pile.
	g_HitCount = g_ImageCount/2;
	g_HitIndexToPileIndexA.resize(g_HitCount);
	g_HitIndexToPileIndexB.resize(g_HitCount);

	int PileIndex = 0;

	const ImageData &FirstImage = g_Images[0];

	PileData Pile;
	Pile.Lo = FirstImage.Lo;
	Pile.Hi = FirstImage.Hi;
	Pile.Deleted = false;

	if (FirstImage.IsA)
		{
		g_HitIndexToPileIndexA[FirstImage.HitIndex] = 0;
		Pile.HitIndexes.push_back(FirstImage.HitIndex);
		}
	else
		g_HitIndexToPileIndexB[FirstImage.HitIndex] = 0;

#if	TRACE
	Log("Image 0 =");
	LogImage(FirstImage);
	Log("\n");
#endif
	
	for (int ImageIndex = 1; ImageIndex < g_ImageCount; ++ImageIndex)
		{
		const ImageData &Image = g_Images[ImageIndex];
#if TRACE
		Log("Pile Lo=%d Hi=%d Image%d ", Pile.Lo, Pile.Hi, ImageIndex);
		LogImage(Image);
		Log("\n");
#endif
		if (Image.Lo > Pile.Hi)
			{
#if	DEBUG
		// Hack to validate coords
			GlobalPosToContig(Pile.Lo);
			GlobalPosToContig(Pile.Hi);
#endif
			g_Piles.push_back(Pile);
#if	TRACE
			Log("Image.Lo > Pile.Hi\n");
			Log("  New pile ");
			++g_PileCount;
			LogPile(PileIndex);
			Log("\n");
#endif
			Pile.Lo = Image.Lo;
			Pile.Hi = Image.Hi;
			Pile.HitIndexes.clear();
			++PileIndex;
			}
		else
			Pile.Hi = Image.Hi;

		if (Image.IsA)
			{
			g_HitIndexToPileIndexA[Image.HitIndex] = PileIndex;
			Pile.HitIndexes.push_back(Image.HitIndex);
			}
		else
			g_HitIndexToPileIndexB[Image.HitIndex] = PileIndex;
		}

#if	DEBUG
// Hack to validate coords
	GlobalPosToContig(Pile.Lo);
	GlobalPosToContig(Pile.Hi);
#endif
	g_Piles.push_back(Pile);
	g_PileCount = (int) g_Piles.size();

	int SizeA = (int) g_HitIndexToPileIndexA.size();
	int SizeB = (int) g_HitIndexToPileIndexB.size();
	if (SizeA != g_HitCount || SizeB != g_HitCount)
		Quit("SizeA=%d SizeB=%d != HitCount=%d", SizeA, SizeB, g_HitCount);
	}
