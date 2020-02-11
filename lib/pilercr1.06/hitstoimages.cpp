#include "pilercr.h"

// Sort by Lo
static bool Cmp_ImageLo(const ImageData &Image1, const ImageData &Image2)
	{
	if (Image1.Lo < Image2.Lo)
		return true;
	return false;
	}

void HitsToImages()
	{
	g_ImageCount = 2*g_HitCount;
	g_Images.resize(g_ImageCount);

	Progress("Converting hits to images");
	for (int HitIndex = 0; HitIndex < g_HitCount; ++HitIndex)
		{
		const DPHit &Hit = g_Hits[HitIndex];

		ImageData &ID1 = g_Images[2*HitIndex];
		ImageData &ID2 = g_Images[2*HitIndex + 1];

		ID1.Lo = Hit.aLo;
		ID1.Hi = Hit.aHi;
		ID1.HitIndex = HitIndex;
		ID1.IsA = true;

		ID2.Lo = Hit.bLo;
		ID2.Hi = Hit.bHi;
		ID2.HitIndex = HitIndex;
		ID2.IsA = false;
		}

	Progress("Sorting images");
	std::sort(g_Images.begin(), g_Images.end(), Cmp_ImageLo);

	if (!g_LogImages)
		return;

	Log("\n");
	Log("\n");
	Log("Image         Lo          Hi    Hit  A\n");
	Log("=====  =========  ==========  =====  =\n");
	for (int ImageIndex = 0; ImageIndex < g_ImageCount; ++ImageIndex)
		{
		const ImageData &Image = g_Images[ImageIndex];
		Log("%5d  %10d  %10d  %5d  %c",
		  ImageIndex, Image.Lo, Image.Hi, Image.HitIndex, Image.IsA ? 'A' : 'B');
		Log("     ");
		LogHit(Image.HitIndex);
		Log("\n");
		}
	}
