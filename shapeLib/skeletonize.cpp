/*
	Code for thinning a binary image using Zhang-Suen algorithm.
	https://answers.opencv.org/question/3207/what-is-a-good-
	thinning-algorithm-for-getting-the-skeleton-of-characters-for-ocr/
	
	Minor rework:
	- Changed variable types.
	- Error correction.
	- Now only standard libs.
	
	Author: Rasmus Vest Nielsen
*/

#include "stdafx.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

void ThinSubiteration1(std::vector< std::vector< unsigned int > > & pSrc,
	std::vector< std::vector< unsigned int > > & pDst, int nRows, int nCols) {
	
	pDst = pSrc;
	
	for (int i = 0; i < nRows; i++) 
	{
		for (int j = 0; j < nCols; j++)
		{
			if (pSrc[i][j] == 1)
			{
				// Get 8 neighbors.
				// Calculate C(p).
				int neighbor0 = (int)pSrc[i - 1][j - 1 ];
				int neighbor1 = (int)pSrc[i - 1][ j    ];
				int neighbor2 = (int)pSrc[i - 1][ j + 1];
				int neighbor3 = (int)pSrc[i    ][ j + 1];
				int neighbor4 = (int)pSrc[i + 1][ j + 1];
				int neighbor5 = (int)pSrc[i + 1][ j    ];
				int neighbor6 = (int)pSrc[i + 1][ j - 1];
				int neighbor7 = (int)pSrc[i    ][ j - 1];
				int C = int(~neighbor1 & (neighbor2 | neighbor3)) +
					int(~neighbor3 & (neighbor4 | neighbor5)) +
					int(~neighbor5 & (neighbor6 | neighbor7)) +
					int(~neighbor7 & (neighbor0 | neighbor1));
				if (C == 1) 
				{
					// Calculate N.
					int N1 = int(neighbor0 | neighbor1) +
						int(neighbor2 | neighbor3) +
						int(neighbor4 | neighbor5) +
						int(neighbor6 | neighbor7);
					int N2 = int(neighbor1 | neighbor2) +
						int(neighbor3 | neighbor4) +
						int(neighbor5 | neighbor6) +
						int(neighbor7 | neighbor0);
					int N = min(N1, N2);
					if ((N == 2) || (N == 3)) 
					{
						// Calculate criteria 3
						int c3 = (neighbor1 | neighbor2 | ~neighbor4) & neighbor3;
						if (c3 == 0) {
							pDst[i][j] = 0;
						}
					}
				}
			}
		}		
	}
	
}

void ThinSubiteration2(std::vector< std::vector< unsigned int > > & pSrc,
	std::vector< std::vector< unsigned int > > & pDst, int nRows, int nCols) {
	
	pDst = pSrc;
	for (int i = 0; i < nRows; i++) 
	{
		for (int j = 0; j < nCols; j++) 
		{
			if (pSrc[i][j] == 1) 
			{
				// Get 8 neighbors.
				// Calculate C(p).
				int neighbor0 = (int)pSrc[i - 1][j - 1];
				int neighbor1 = (int)pSrc[i - 1][j    ];
				int neighbor2 = (int)pSrc[i - 1][j + 1];
				int neighbor3 = (int)pSrc[i    ][j + 1];
				int neighbor4 = (int)pSrc[i + 1][j + 1];
				int neighbor5 = (int)pSrc[i + 1][j    ];
				int neighbor6 = (int)pSrc[i + 1][j - 1];
				int neighbor7 = (int)pSrc[i    ][j - 1];
				int C = int(~neighbor1 & (neighbor2 | neighbor3)) +
					int(~neighbor3 & (neighbor4 | neighbor5)) +
					int(~neighbor5 & (neighbor6 | neighbor7)) +
					int(~neighbor7 & (neighbor0 | neighbor1));
				if (C == 1) 
				{
					// Calculate N.
					int N1 = int(neighbor0 | neighbor1) +
						int(neighbor2 | neighbor3) +
						int(neighbor4 | neighbor5) +
						int(neighbor6 | neighbor7);
					int N2 = int(neighbor1 | neighbor2) +
						int(neighbor3 | neighbor4) +
						int(neighbor5 | neighbor6) +
						int(neighbor7 | neighbor0);
					int N = min(N1, N2);
					if ((N == 2) || (N == 3)) 
					{
						int E = (neighbor5 | neighbor6 | ~neighbor0) & neighbor7;
						if (E == 0) {
							pDst[i][j] = 0;
						}
					}
				}
			}
		}
	}
}

void skeletonize(std::vector< std::vector< unsigned int > > &inArray,
	std::vector< std::vector< unsigned int > > &outArray, int nRows, int nCols)
{
	bool bDone = false;
	
	// Pad source by delta to avoid boundary pixels that
	// cannot have 8 neighbours.
	unsigned int delta = 2; // A value of 2 is arbitrary.
	std::vector< std::vector< unsigned int > >
		p_enlarged_src(nRows + delta, std::vector< unsigned int >(nCols + delta, 0));
	
	for (int i = 0; i < (nRows + delta); i++) {
		p_enlarged_src[i][0] = 0;
		p_enlarged_src[i][nCols + 1] = 0;
	}
	for (int j = 0; j < (nCols + delta); j++) {
		p_enlarged_src[0][j] = 0;
		p_enlarged_src[nRows + 1][j] = 0;
	}
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols; j++) {
			if (inArray[i][j] >= 1) {
				p_enlarged_src[i + 1][j + 1] = 1;
			}
			else
				p_enlarged_src[i + 1][j + 1] = 0;
		}
	}
	
	/*
	// Display binary array.
	std::cout << "Padded zero array: " << endl;
	for (int i = 0; i < nRows + delta; i++)
	{
		for (int j = 0; j < nCols + delta; j++)
		{
			std::cout << static_cast<unsigned>
				(p_enlarged_src[i][j]) << ", ";
		}
		std::cout << std::endl;
	}
	*/

	// Start to thin.
	std::vector< std::vector< unsigned int > >
		p_thinMat1(nRows + delta, std::vector< unsigned int >(nCols + delta, 0));
	std::vector< std::vector< unsigned int > >
		p_thinMat2(nRows + delta, std::vector< unsigned int >(nCols + delta, 0));
	std::vector< std::vector< unsigned int > >
		p_cmp(nRows + delta, std::vector< unsigned int >(nCols + delta, 0));
	
	while (bDone != true) {
		// Sub-iteration 1.
		ThinSubiteration1(p_enlarged_src, p_thinMat1, nRows, nCols);

		// Sub-iteration 2.
		ThinSubiteration2(p_thinMat1, p_thinMat2, nRows, nCols);
		
		// Check if any pixels were set during the last
		// 2 iterations.
		if (p_enlarged_src == p_thinMat2) 
		{
			bDone = true;
		}
		// Copy.
		p_enlarged_src = p_thinMat2;
	}

	// Copy results.
	outArray = p_enlarged_src;

	/* C-style arrays!
	unsigned char ** p_enlarged_src = new unsigned char *[nRows + delta];
	for (int i = 0; i < nRows + delta; i++)
	{
	p_enlarged_src[i] = new unsigned char[nCols + delta];
	}

	// Start to thin.
	unsigned char ** p_thinmat1 = new unsigned char *[nRows + delta];
	unsigned char ** p_thinmat2 = new unsigned char *[nRows + delta];
	unsigned char ** p_cmp = new unsigned char *[nRows + delta];
	for (int i = 0; i < nRows + delta; i++)
	{
	p_thinmat1[i] = new unsigned char[nCols + delta];
	p_thinmat2[i] = new unsigned char[nCols + delta];
	p_cmp[i] = new unsigned char[nCols + delta];
	}
	*/
}

