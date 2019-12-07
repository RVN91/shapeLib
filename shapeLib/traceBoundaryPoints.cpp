/*
	Tracing Boundary in 2D Image Using Moore Neighborhood Approach. 
	Output is a set of continuous points along the boundary of object 
	in 2D image. The current implementation supports only one connected 
	object in the image.

	https://www.codeproject.com/Articles/1105045/Tracing-Boundary-
	in-D-Image-Using-Moore-Neighborho

	Edited and corrected:
	- Inverted y-axis (normal in image processing).
	- X and Y coordinates are switched to follow row and
	  column notation (normal for arrays/matrices).

	NOTES!
	- Particles with weird corners fail if the corners
	  are detected during the first scan.

	Author: Rasmus Vest Nielsen
*/

#include "stdafx.h"
#include "traceBoundaryPoints.h"
#include <vector>

TraceBoundaryPoints::TraceBoundaryPoints()
{
	// Do nothing
}
TraceBoundaryPoints::~TraceBoundaryPoints()
{
	//Do nothing
}
/*
* Description - Get the continous boundary points
* Parameters
* InputImage	- Input image
* Width_i		- Width of the image
* Height_i		- Height of Image
* BoundaryPoints - Vector of boundary points (output)
*/
void TraceBoundaryPoints::GetContinousBoundaryPoints(std::vector < unsigned int > InputImage, 
	int Width_i, 
	int Height_i,
	std::vector<Point2D> &BoundaryPoints)
{
	int nImageSize = Width_i * Height_i;
	if (InputImage.size() != 0)
	{
		int Offset[8][2] = {
			{ -1, -1 },       //  +----------+----------+----------+
			{ 0, -1 },		  //  |          |          |          |
			{ 1, -1 },        //  |(x-1,y-1) | (x,y-1)  |(x+1,y-1) |
			{ 1, 0 },         //  +----------+----------+----------+
			{ 1, 1 },         //  |(x-1,y)   |  (x,y)   |(x+1,y)   |
			{ 0, 1 },         //  |          |          |          |
			{ -1, 1 },        //  +----------+----------+----------+
			{ -1, 0 }         //  |          | (x,y+1)  |(x+1,y+1) |
		};                    //  |(x-1,y+1) |          |          |
							  //  +----------+----------+----------+
		const int NEIGHBOR_COUNT = 8;
		Point2D BoundaryPixelCord;
		Point2D BoundaryStartingPixelCord;
		Point2D BacktrackedPixelCord;
		int BackTrackedPixelOffset[1][2] = { { 0,0 } };
		bool bIsBoundaryFound = false;
		bool bIsStartingBoundaryPixelFound = false;
		for (int Idx = 0; Idx < nImageSize; ++Idx) // getting the starting pixel of boundary
		{
			if (0 != InputImage[Idx])
			{
				BoundaryPixelCord.Y= Idx % Width_i;
				BoundaryPixelCord.X = Idx / Width_i;
				BoundaryStartingPixelCord = BoundaryPixelCord;
				BacktrackedPixelCord.Y = (Idx - 1) % Width_i;
				BacktrackedPixelCord.X = (Idx - 1) / Width_i;
				BackTrackedPixelOffset[0][0] = BacktrackedPixelCord.Y - BoundaryPixelCord.Y;
				BackTrackedPixelOffset[0][1] = BacktrackedPixelCord.X - BoundaryPixelCord.X;
				BoundaryPoints.push_back(BoundaryPixelCord);
				bIsStartingBoundaryPixelFound = true;
				break;
			}
		}
		Point2D CurrentBoundaryCheckingPixelCord;
		Point2D PrevBoundaryCheckingPixxelCord;
		if (!bIsStartingBoundaryPixelFound)
		{
			BoundaryPoints.pop_back();
		}
		while (true && bIsStartingBoundaryPixelFound)
		{
			int CurrentBackTrackedPixelOffsetInd = -1;
			for (int Ind = 0; Ind < NEIGHBOR_COUNT; ++Ind)
			{
				if (BackTrackedPixelOffset[0][0] == Offset[Ind][0] &&
					BackTrackedPixelOffset[0][1] == Offset[Ind][1])
				{
					CurrentBackTrackedPixelOffsetInd = Ind;// Finding the bracktracked pixel's offset index
					break;
				}
			}
			int Loop = 0;
			while (Loop < (NEIGHBOR_COUNT - 1) && CurrentBackTrackedPixelOffsetInd != -1)
			{
				int OffsetIndex = (CurrentBackTrackedPixelOffsetInd + 1) % NEIGHBOR_COUNT;
				CurrentBoundaryCheckingPixelCord.Y = BoundaryPixelCord.Y + Offset[OffsetIndex][0];
				CurrentBoundaryCheckingPixelCord.X = BoundaryPixelCord.X + Offset[OffsetIndex][1];
				int ImageIndex = CurrentBoundaryCheckingPixelCord.X * Width_i + CurrentBoundaryCheckingPixelCord.Y;
				if (0 != InputImage[ImageIndex])// finding the next boundary pixel
				{
					BoundaryPixelCord = CurrentBoundaryCheckingPixelCord;
					BacktrackedPixelCord = PrevBoundaryCheckingPixxelCord;
					BackTrackedPixelOffset[0][0] = BacktrackedPixelCord.Y - BoundaryPixelCord.Y;
					BackTrackedPixelOffset[0][1] = BacktrackedPixelCord.X - BoundaryPixelCord.X;
					BoundaryPoints.push_back(BoundaryPixelCord);
					break;
				}
				PrevBoundaryCheckingPixxelCord = CurrentBoundaryCheckingPixelCord;
				CurrentBackTrackedPixelOffsetInd += 1;
				Loop++;
			}
			if (BoundaryPixelCord.Y == BoundaryStartingPixelCord.Y &&
				BoundaryPixelCord.X == BoundaryStartingPixelCord.X) // if the current pixel = starting pixel
			{
				BoundaryPoints.pop_back();
				bIsBoundaryFound = true;
				break;
			}
		}
		if (!bIsBoundaryFound) // If there is no connected boundary clear the list
		{
			std::cout << "No connected border!" << std::endl;
			BoundaryPoints.clear();
		}
	}
}