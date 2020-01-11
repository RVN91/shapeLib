/*
	Load and assign particle information to variables from sIMPLE binary
	file.

	Author: Rasmus Vest Nielsen
*/

#include "stdafx.h"

#include <vector>
#include <stdlib.h> // To exit if file is not found.
#include <fstream>
#include <algorithm>
#include "convexHull.h"
#include "traceBoundaryPoints.h"

/*
	Safely reads and casts 8-bit (1 byte) to char from file stream assuming
	little endianness. Reading 1 byte memmory blocks are chosen to ensure
	independence from system types.
*/
char readU8(std::istream& file, long *bytePosition)
{
	file.seekg(*bytePosition);
	char * memmoryBlock;
	int nBytes = 1;
	memmoryBlock = new char[nBytes];
	file.read(memmoryBlock, nBytes);

	*bytePosition += 1;

	return * memmoryBlock;
}

/*
	Safely reads and casts unsigned 32-bit (4 bytes) to integer from file stream
	assuming little endianness. Reading 1 byte memmory blocks are chosen to ensure
	independence from system types.
*/
unsigned int readU32(std::istream& file, long *bytePosition)
{
	file.seekg(*bytePosition);
	char * memmoryBlock;
	int nBytes = 4;
	memmoryBlock = new char [nBytes];
	file.read(memmoryBlock, nBytes);
	unsigned int * val = (unsigned int *) memmoryBlock;

	*bytePosition += 4;

	return *val;
}

/*
	Safely reads and casts signed 32-bit (4 bytes) to integer from file stream
	assuming little endianness. Reading 1 byte memmory blocks are chosen to ensure
	independence from system types.
*/
signed int read32(std::istream& file, long *bytePosition)
{
	file.seekg(*bytePosition);
	char * memmoryBlock;
	int nBytes = 4;
	memmoryBlock = new char [nBytes];
	file.read(memmoryBlock, nBytes);
	signed int * val = (signed int *) memmoryBlock;

	*bytePosition += 4;

	return *val;
}

/*
	Safely reads and casts sined 32-bit (4 bytes) to float from file stream
	assuming little endianness. Reading 1 byte memmory blocks are chosen to ensure
	independence from system types.
*/
float read32float(std::istream& file, long *bytePosition)
{
	file.seekg(*bytePosition);
	char * memmoryBlock;
	int nBytes = 4;
	memmoryBlock = new char [nBytes];
	file.read(memmoryBlock, nBytes);

	//bitset<32> x(*memmoryBlock);
	//cout << x << endl;

	float * val = (float*) memmoryBlock;

	*bytePosition += 4;

	return *val;
}

/*
	Checks if local machine is small or big endian. Most modern CPU's use small
	endian so it might be a waste of time.

	Currently not in use...
*/
void endianTest()
{
	unsigned char test[2] = {1,0};
	short x;
	x = *(short*) test;

	if(x == 1)
	{
		std::cout << "System is little endian!" << std::endl;
	}
	else
	{
		std::cout << "System is big endian!" << std::endl;
	}
}


/*
	Loads particle count
*/
void countParticles(const char *fileName, int *nParticles, 
	long *bytePosition)
{
	std::ifstream particleFile;
	particleFile.open(fileName, std::ios::binary);
	if (!particleFile.is_open())
	{
		std::cout << "ERROR: Cannot open file!" << std::endl;
		exit(EXIT_FAILURE);
	} else if (*bytePosition != 0)
	{
		std::cout << "ERROR: Byte position needs to be 0!" << std::endl;
		exit(EXIT_FAILURE);
	}

	*nParticles = readU32(particleFile, bytePosition);
}


/*
	Loads the pixel numbers of the current particle.
*/
void getPixelCount(const char * fileName, 
	long *bytePosition,
	int *nPixels)
{
	std::ifstream particleFile;
	particleFile.open(fileName, std::ios::binary);
	if (!particleFile.is_open())
	{
		std::cout << "Cannot open file!" << std::endl;
		exit(EXIT_FAILURE);
	}

	unsigned int rsC = readU32(particleFile, bytePosition);
	char booleanChar = readU8(particleFile, bytePosition);
	*nPixels = readU32(particleFile, bytePosition);
	particleFile.close();
}

/*
	Loads pixels from binary file into vectors of vectors.
*/
void getParticlePixels(const char * fileName,
	std::vector<unsigned int> &xPixels, 
	std::vector<unsigned int> &yPixels,
	float *minorDim,
	float *majorDim,
	long *bytePosition,
	int *nPixels)
{
	std::ifstream particleFile;
	particleFile.open(fileName, std::ios::binary);
	if (!particleFile.is_open())
	{
		std::cout << "Cannot open file!" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	/*
	unsigned int rsC = readU32(particleFile, bytePosition);
	char booleanChar = readU8(particleFile, bytePosition);
	unsigned int nPixels = readU32(particleFile, bytePosition);
	*/
	
	// x and y pixels coordinates.
	for (unsigned int j = 0; j <= *nPixels; j++)
	{
		//xPixels[j] = read32(particleFile, bytePosition);
		xPixels.push_back(read32(particleFile, bytePosition));
	}

	for (unsigned int j = 0; j <= *nPixels; j++)
	{
		//yPixels[j] = read32(particleFile, bytePosition);
		yPixels.push_back(read32(particleFile, bytePosition));
	}

	unsigned int xMean = readU32(particleFile, bytePosition);
	unsigned int yMean = readU32(particleFile, bytePosition);
	float majorDimension = read32float(particleFile, bytePosition);
	float minorDimension = read32float(particleFile, bytePosition);
	float volume = read32float(particleFile, bytePosition);
	float mass = read32float(particleFile, bytePosition);
	float pno = readU32(particleFile, bytePosition);
	float mpno = readU32(particleFile, bytePosition);
	
	*minorDim = minorDimension;
	*majorDim = majorDimension;

	// Delphi strings are not fixed in size but a comes with a byte
	// length defined in the first byte of the string.
	char strLength = readU8(particleFile, bytePosition);
	//printf("Str. length: %d \n", strLength)
	*bytePosition += strLength;
	
	particleFile.close();
}

/*
	Convert particle pixel vectors to Point.
*/
std::vector <Point> convertParticlesPixelsToPoints(std::pair < std::vector < unsigned int >,
	std::vector < unsigned int > > particlePixels)
{
	// Extract particle vectors
	std::vector< unsigned int > xPixels = particlePixels.first;
	std::vector< unsigned int > yPixels = particlePixels.second;
	std::size_t nPixels = particlePixels.first.size();

	// The reverse order of x and y is due to images are defined
	// with x as rows and y as columns.
	std::vector < Point > points;

	// Set pixel values in Point vector.
	// Point object requires casting to double.
	for (int k = 0; k <= nPixels - 1; k++)
	{
		points.push_back(Point{ (double)xPixels[k],
			(double)yPixels[k] });
		//std::cout << "x, y: " << xPixels[k] << ", " << yPixels[k] << std::endl;
	}

	//std::cout << "Size: " << points.size() << std::endl;

	return points;
}

/*
Convert particle pixels vectors to Point2D.
*/
std::vector <Point2D> convertParticlesPixelsToPoints2D(std::pair < std::vector < unsigned int >,
	std::vector < unsigned int > > particlePixels)
{
	// Extract particle vectors
	std::vector< unsigned int > xPixels = particlePixels.first;
	std::vector< unsigned int > yPixels = particlePixels.second;
	std::size_t nPixels = particlePixels.first.size();

	// The reverse order of x and y is due to images are defined
	// with x as rows and y as columns.
	std::vector < Point2D > points;

	// Set pixel values in Point vector.
	for (int k = 0; k <= nPixels - 1; k++)
	{
		points.push_back(Point2D{ (int)xPixels[k],
			(int)yPixels[k] });
		//std::cout << "x, y: " << xPixels[k] << ", " << yPixels[k] << std::endl;
	}

	//std::cout << "Size: " << points.size() << std::endl;

	return points;
}

/*
	Offset pixel vectors by 10.
*/
void offsetPixelVectors(std::pair < std::vector < unsigned int >, 
	std::vector < unsigned int > > &particlePixels, 
	std::pair < std::vector < unsigned int >, std::vector < unsigned int > > &offsetPixels)
{
	// Extract particle vectors
	std::vector< unsigned int > xPixels = particlePixels.first;
	std::vector< unsigned int > yPixels = particlePixels.second;
	int nPixels = particlePixels.first.size();

	// Find minimum and maximum values.
	auto xMin = std::min_element(std::begin(xPixels),
		std::end(xPixels));
	auto xMax = std::max_element(std::begin(xPixels),
		std::end(xPixels));
	auto yMin = std::min_element(std::begin(yPixels),
		std::end(yPixels));
	auto yMax = std::max_element(std::begin(yPixels),
		std::end(yPixels));

	// Offset pixel values
	std::vector < unsigned int >  xPixelsOffset(nPixels);
	std::vector < unsigned int >  yPixelsOffset(nPixels);

	for (unsigned int j = 0; j <= nPixels - 1; j++)
	{
		xPixelsOffset[j] = xPixels[j] - *xMin + 5;
	}

	for (unsigned int j = 0; j <= nPixels - 1; j++)
	{
		yPixelsOffset[j] = yPixels[j] - *yMin + 5;
	}

	// Free memory.
	xPixels.clear();
	yPixels.clear();

	// Pair up offset vectors.
	offsetPixels = std::make_pair(xPixelsOffset, yPixelsOffset);
}

/*
	Calculate new pixel boundaries of offset pixel vectors.
*/
void calculatePixelBoundaries(std::pair < std::vector < unsigned int >,
	std::vector < unsigned int > > &particlePixels,
	std::pair <int, int> &boundaryXY)
{	
	// Extract particle vectors
	std::vector< unsigned int > xPixels = particlePixels.first;
	std::vector< unsigned int > yPixels = particlePixels.second;
	
	// Find minimum and maximum values.
	auto xMin = std::min_element(std::begin(xPixels),
		std::end(xPixels));
	auto xMax = std::max_element(std::begin(xPixels),
		std::end(xPixels));
	auto yMin = std::min_element(std::begin(yPixels),
		std::end(yPixels));
	auto yMax = std::max_element(std::begin(yPixels),
		std::end(yPixels));

	// Calculate new boundaries.
	int boundarySizeX = *xMax - *xMin + 10;
	int boundarySizeY = *yMax - *yMin + 10;

	boundaryXY = std::make_pair(boundarySizeX, boundarySizeY);
}

/*
	Accepts a pair of vectors of x, y coordinates and transforms
	them into a binary 2D vector of vectors (2D matrix).
*/
void pixelsTo2DBinaryVector(
	std::pair < std::vector < unsigned int >,
	std::vector < unsigned int > > & pixelsXY,
	std::vector< std::vector< unsigned int > > & particle2DVector)
{
	// Set particle pixels to 1.
	for (unsigned int i = 0; i < pixelsXY.first.size(); i++)
	{
		particle2DVector[pixelsXY.first[i]][pixelsXY.second[i]] = 1;
	}
}

/*
	Transforms 2D array into a binary 1D vector.
*/
void matrixTo1DArray(
	std::vector< std::vector< unsigned int > > & matrix,
	std::vector< unsigned int > & array)
{
	// Convert to 1D array.
	for (unsigned int i = 0; i < matrix.size(); i++)
	{
		for (unsigned int j = 0; j < matrix[0].size(); j++)
		{
			array.push_back(matrix[i][j]);
		}
	}
}

/*
	Accepts a vector of x, y coordinates and transforms
	them into a binary 2D vector of vectors (2D matrix).
*/
std::vector< std::vector< unsigned int > > vectorTo2DBinaryVector(
	std::vector < Point > & pixelsXY,
	std::vector< std::vector< unsigned int > > & particle2DVector,
	std::size_t nRows, std::size_t nCols)
{
	// Set particle pixels to 1.
	for (unsigned int i = 0; i < nRows; i++)
	{
		particle2DVector[pixelsXY[i].x][pixelsXY[i].y] = 1;
	}
	return particle2DVector;
}