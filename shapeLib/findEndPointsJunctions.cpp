/*
Finds end points and junctions in skeletonized (binary) 
images.

https://answers.opencv.org/question/104324/skeleton-image/

Minor changes to fit the goals of this project.
- No opencv container.
- Reduced to only use of standard libs.

Author: Rasmus Vest Nielsen
*/

#include "stdafx.h"

#include <vector>
#include <iostream>
#include "convexHull.h"

std::vector <Point> findEndPoints(std::vector< std::vector< unsigned int > > skeleton)
{
	// Let's collect "interesting points", 
	// if there's only one neighbour, it's an endpoint, 
	// if it has 3, it's a bifurcation(or T-junction), 
	// 4 neighbours denote an X-like crossing.
	std::vector<Point> one, three, four;

	unsigned int nRows = skeleton.size();
	unsigned int nCols = skeleton[0].size();
	//std::cout << nRows << ", " << nCols << std::endl;
	// since we're checking 1 pixel neighbourhood, 
	// we need to spare 1 pixel border on each side: 
	for (unsigned int r = 1; r < nRows - 1; r++) {
		for (unsigned int c = 1; c < nCols - 1; c++) {
			unsigned int cen = skeleton[r][c];
			if (cen == 0) continue; // Background.
									// Now we just walk in a circle around the center pixel, 
									// and collect neighbours(starting top-left):    
			int neighbours = skeleton[r - 1][c - 1] + skeleton[r - 1][c] + skeleton[r - 1][c + 1]
						   + skeleton[r    ][c - 1]					     + skeleton[r    ][c + 1]
						   + skeleton[r + 1][c - 1] + skeleton[r + 1][c] + skeleton[r + 1][c + 1];

			if (neighbours == 1)
			{
				one.push_back(Point{ (double)r,
					(double)c });
			}

			if (neighbours == 3)
			{
				three.push_back(Point{ (double)r,
					(double)c });
			}

			if (neighbours == 4)
			{
				four.push_back(Point{ (double)r,
					(double)c });
			}
		}
	}
	
	return one;
}
