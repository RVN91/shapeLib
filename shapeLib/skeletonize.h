/*
* Skeletonization.h
*
*  Code for thinning a binary image using Zhang-Suen algorithm.
*  https://answers.opencv.org/question/3207/what-is-a-good-
*  thinning-algorithm-for-getting-the-skeleton-of-characters-for-ocr/
*/

// Set guard.
#ifndef SKELETONIZE_H
#define SKELETONIZE_H

#include <iostream>

//void ThinSubiteration2(Mat & pSrc, Mat & pDst);
//void ThinSubiteration1(Mat & pSrc, Mat & pDst);
void skeletonize(std::vector< std::vector< unsigned int > > & inArray,
	std::vector< std::vector< unsigned int > > & outArray, int rows, int cols);


#endif // End of guard