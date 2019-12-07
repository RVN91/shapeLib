/*
Declarations of functions for loadParticle.cpp.

Author: Rasmus Vest Nielsen
*/

#include <iostream>
#include <vector>
#include "convexHull.h"

// Set guard
#ifndef LOADPARTICLE_H
#define LOADPARTICLE_H

// Function declaration.
char readU8(std::istream& file, long *bytePosition);

unsigned int readU32(std::istream& file, long *bytePosition);

signed int read32(std::istream& file, long *bytePosition);

float read32float(std::istream& file, long *bytePosition);

void endianTest();

void countParticles(const char *fileName, int *nParticles,
	long *bytePosition);

void getPixelCount(const char * fileName,
	long *bytePosition,
	int *nPixels);

void getParticlePixels(const char * fileName,
	std::vector<unsigned int> &xPixels,
	std::vector<unsigned int> &yPixels,
	float *minorDim,
	float *majorDim,
	long *bytePosition,
	int *nPixels);

std::vector <Point> convertParticlesPixelsToPoints(std::pair < std::vector < unsigned int >,
	std::vector < unsigned int > > particlePixels);

void offsetPixelVectors(std::pair < std::vector < unsigned int >,
	std::vector < unsigned int > > &particlePixels,
	std::pair < std::vector < unsigned int >, std::vector < unsigned int > > &offsetPixels);

void calculatePixelBoundaries(std::pair < std::vector < unsigned int >,
	std::vector < unsigned int > > &particlePixels,
	std::pair <int, int> &boundaryXY);

void pixelsTo2DBinaryVector(
	std::pair < std::vector < unsigned int >,
	std::vector < unsigned int > > & pixelsXY,
	std::vector< std::vector< unsigned int > > & particle2DVector);

std::vector< std::vector< unsigned int > > vectorTo2DBinaryVector(
	std::vector < Point > & pixelsXY,
	std::vector< std::vector< unsigned int > > & particle2DVector);

void matrixTo1DArray(
	std::vector< std::vector< unsigned int > > & matrix,
	std::vector< unsigned int > & array);

#endif // End of guard
