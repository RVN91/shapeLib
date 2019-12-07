// shapeLib.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// My libraries.
#include "loadParticle.h"
#include "traceBoundaryPoints.h"

// The extern "C" construct prevents the compiler to add decorations on 
// the functions' names in the DLL/SO and it is necessary while using C++.
#define DESCRIPTORS extern "C" __declspec(dllexport)

/*
Creates a resultfile with headers.
*/

DESCRIPTORS const char * createResultFile(bool debug)
{
	const char * inputFileNameChar;
	const char * csvFileNameChar;

	if (debug == true)
	{
		inputFileNameChar = "T1_1.pmp";
		csvFileNameChar = "results.csv";
	}
	else
	{
		std::cout << "Enter file name: " << std::endl;
		std::string inputFileName;
		std::cin >> inputFileName;
		inputFileNameChar = inputFileName.c_str();
		std::cout << "Name of input file: " << inputFileName << std::endl;
		std::cout << "Set name for output file name (add .csv to the end of the file name!): " << std::endl;
		std::string csvFileName;
		std::cin >> csvFileName;
		csvFileNameChar = csvFileName.c_str();
	}

	// Remove previous files.
	remove(csvFileNameChar);

	// Write headers to .csv file.
	std::ofstream file(csvFileNameChar, std::ios::app);
	if (file.is_open())
	{
		file << "Particle_no" << ",";
		file << "Area" << ",";
		file << "Convex area" << ",";
		file << "Perimeter" << ",";
		file << "Convex_Perimeter" << ",";
		file << "Max_Feret" << ",";
		file << "Min_Feret" << ",";
		file << "Aspect_ratio" << ",";
		file << "Feret_aspect_ratio" << ",";
		file << "Circularity" << ",";
		file << "Convexity" << ",";
		file << "Solidity" << ",";
		file << "Equi_diameter" << ",";
		file << "Compactness" << ",";
		file << "Fiber_length" << ",";
		file << "Curl" << ",";
		file << "Average_fiber_diameter" << ",";
		file << "Fiber_elongation" << std::endl;
		file.close();
	}
	return inputFileNameChar;
}

DESCRIPTORS int getParticleCount(const char *inputFileNameChar, long *bytePosition)
{
	int nParticles = 0;

	countParticles(inputFileNameChar, &nParticles, bytePosition);
	std::cout << "Number of particles: " << nParticles << std::endl;
	
	return nParticles;
}

DESCRIPTORS void getPixelsCount(const char *inputFileNameChar, 
	long *bytePosition,
	int *nPixels)
{
	getPixelCount(inputFileNameChar, bytePosition, nPixels);
	std::cout << *nPixels << std::endl;
}

DESCRIPTORS void getParticle(const char *inputFileNameChar, 
	long *bytePosition,
	int *nPixels,
	int *xPixels,
	int *yPixels)
{

	// Get particle pixels and particle dimensions.
	// Minor and major dimension of particles are already calculated
	// in sIMPLe sofware, hence they can be used to calculate aspect ratio.
	float *minorDimensions = 0;
	float *majorDimensions = 0;
	
	std::vector<unsigned int> xPixelsVec; 
	std::vector<unsigned int> yPixelsVec;
	
	getParticlePixels(inputFileNameChar, xPixelsVec, yPixelsVec, minorDimensions,
		majorDimensions, bytePosition, nPixels);

	// Insert pixel coordinates into pointers.
	for (int i = 0; i < *nPixels; i++)
	{
		xPixels[i] = xPixelsVec[i];
		yPixels[i] = yPixelsVec[i];
	}
}

DESCRIPTORS void traceBoundary(int *particlePointers, 
	int *nPixels,
	int nRows,
	int nCols,
	int **boundaryPointsArrayX,
	int **boundaryPointsArrayY,
	int *size)
{
	// Convert array of pointers to pointers into
	// vector.
	std::vector <unsigned int> particleVector;
	
	for (int i = 0; i < *nPixels; i++)
	{
		particleVector.push_back(particlePointers[i]);
		//std::cout << particleVector[i] << std::endl;
	}
	
	std::vector <Point2D> boundaryPoints;
	TraceBoundaryPoints boundObj;
	boundObj.GetContinousBoundaryPoints(particleVector, 
		nCols, nRows, boundaryPoints);
	
	// Write to array.
	*size = boundaryPoints.size();
	*boundaryPointsArrayX = new int[*size];
	*boundaryPointsArrayY = new int[*size];
	for (int i = 0; i < *size; i++)
	{
		(*boundaryPointsArrayX)[i] = boundaryPoints[i].X;
		(*boundaryPointsArrayY)[i] = boundaryPoints[i].Y;
	}
}

DESCRIPTORS void convexHull(int *xPixels, int *yPixels, int nPixels,
	int **xConvex, int **yConvex, int *size)
{
	// Calculate convex hull.
	// Transform particle pixels vectors into points.
	// The reverse order of x and y is due to images are defined
	// with x as rows and y as columns.
	std::vector < Point > points;
	
	// Set pixel values in Point vector.
	// Point object requires casting to double.
	for (int k = 0; k < nPixels; k++)
	{
		points.push_back(Point{ (double)xPixels[k],
			(double)yPixels[k] });
		//std::cout << "x, y: " << xPixels[k] << ", " << yPixels[k] << std::endl;
	}

	// Conversion from unsigned int to double of points 
	// (NEEDS TO BE CHANGED).
	std::vector <Point> convexHull = makeConvexHullPresorted(points);

	// Insert pixel coordinates into pointers.
	*size = convexHull.size();
	*xConvex = new int[*size];
	*yConvex = new int[*size];
	for (int i = 0; i < *size; i++)
	{
		(*xConvex)[i] = (int)points[i].x;
		(*yConvex)[i] = (int)points[i].y;
	}
}


// Test function!
DESCRIPTORS void testArray(int **array, int *size)
{
	*size = 10;
	*array = new int[*size];
	for (int i = 0; i < *size; i++)
	{
		(*array)[i] = i;
	}
}

DESCRIPTORS void free_mem_float(float** a)
{
	delete[] a;
}
