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
#include "skeletonize.h"

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
	int *yPixels,
	float *minorDimension,
	float *majorDimension)
{

	// Get particle pixels and particle dimensions.
	// Minor and major dimension of particles are already calculated
	// in sIMPLe sofware, hence they can be used to calculate aspect ratio.
	std::vector<unsigned int> xPixelsVec; 
	std::vector<unsigned int> yPixelsVec;
	
	getParticlePixels(inputFileNameChar, xPixelsVec, yPixelsVec, minorDimension,
		majorDimension, bytePosition, nPixels);

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
	std::vector <Point> convexHull = makeConvexHull(points);

	// Insert pixel coordinates into pointers.
	*size = convexHull.size();
	*xConvex = new int[*size];
	*yConvex = new int[*size];
	for (int i = 0; i < *size; i++)
	{
		(*xConvex)[i] = (int)convexHull[i].x;
		(*yConvex)[i] = (int)convexHull[i].y;
	}
}

/*
Distance between two points.
*/
double dist(Point a, Point b)
{
	return sqrt((a.x - b.x) * (a.x - b.x)
		+ (a.y - b.y) * (a.y - b.y));
}

/*
Perimeter of a set of points describing the contour
of a particle.
*/
DESCRIPTORS void findPerimeter(int *x, int *y, 
	int *nPoints, double *perimeter)
{
	// Transform pointers into vector,
	std::vector<Point> pointsVec;
	
	// Set pixel values in Point vector.
	// Point object requires casting to double.
	for (int k = 0; k < *nPoints; k++)
	{
		pointsVec.push_back(Point{ (double)x[k],
			(double)y[k] });
		//std::cout << "x, y: " << xPixels[k] << ", " << yPixels[k] << std::endl;
	}
	
	*perimeter = 0.0;
	// Find the distance between adjacent points 
	for (int i = 0; i < pointsVec.size() - 1; i++) {
		*perimeter += dist(pointsVec[i],
			pointsVec[i + 1]);
	}

	// Add the distance between first and last point 
	*perimeter += dist(pointsVec[0],
		pointsVec[pointsVec.size() - 1]);
}

/*
Calculate area of polygon from points.

https://www.mathopenref.com/coordpolygonarea2.html
*/
double polygonArea(std::vector <Point> &points)
{
	// Number of points.
	std::size_t n = points.size();

	// Initialze area.
	double area = 0.0;

	// Calculate value of shoelace formula.
	std::size_t j = n - 1;
	for (int i = 0; i < n; i++)
	{
		area += (points[j].x + points[i].x) *
			(points[j].y - points[i].y);
		j = i;  // j is previous vertex to i.
	}

	// Return absolute value.
	return abs(area / 2.0);
}


DESCRIPTORS void calculateArea(
	int *x, 
	int *y, 
	int *nPoints, 
	double *area)
{
	// Transform pointers into vector,
	std::vector<Point> points;

	// Set pixel values in Point vector.
	// Point object requires casting to double.
	for (int k = 0; k < *nPoints; k++)
	{
		points.push_back(Point{ (double)x[k],
			(double)y[k] });
		//std::cout << "x, y: " << xPixels[k] << ", " << yPixels[k] << std::endl;
	}
	
	*area = polygonArea(points);
}

DESCRIPTORS void calculateSolidity(
	double *area, 
	double *convexArea,
	double *solidity)
{
	*solidity = *area / *convexArea;
}

/*
	Calculate the absolute area from 3 vertices in the 
	2D plane.

			   C
			   1
			  1 1
			 1   1
			1     1
		   1       1
		  1         1
		 1           1
		111111111111111
	   A               B

	https://www.mathopenref.com/coordtrianglearea.html
*/
double areaTriangle(Point a, Point b, Point c)
{
	double area;
	area = 0.5 * (a.x * (b.y - c.y) + b.x *
		(c.y - a.y) + c.x * (a.y - b.y));
	double absArea = std::abs(area);
	return absArea;
}

/*
Calculates the distance between two vertices
in the 2-dimensional plane.
A ----- B
*/
double distance(Point a, Point b)
{
	double distance;
	distance = std::sqrt(std::pow((b.x - a.x), 2) + std::pow((b.y - a.y), 2));
	return distance;
}

/*
Calculate the minimum and maximum Feret diameters by from the
antipodal points of the convex hull.

Best knowledge of original reference mentioning Feret diameters:
L.R. Feret, “La Grosseur des Grains” (Assoc. Intern. Essais Math.
2D, Zurich, 1931).

The minimum Feret diameter is equal to the to the size measured
when passing the object through sieve analysis.

The maximum Feret diameter is useful as an ESTIMATE of the length
of elongated objects.

This algorithm is based on the rotating calipers algorithm but
the problem with this algorithm is that the width estimated for
very elongated objects is not very accurate. The orientation that
produces the shortest projection could be up to 1 degree away from
the optimal orientation, meaning that the estimated width is
length*sin(pi/180) too large. This does not sound like much, but if
the aspect ratio is 100, meaning the length is 100 times the width,
we can overestimate the width by up to 175%!

Antipodal points are a pair of two points on the boundary of the
convex hull, such there exists two parallel lines, one through
each point, and every other point of the polygon lies within these
two lines.

Area is used as a proxy for distance between two points (square
roots and shit is a computationally expensive task compared to
area which just requires multiplication and addition).

The actual minumum and maximum feret are found by scanning the
antipodal points and finding the minimum and maximum distance
between the antipodal points.

Notes and improvements:
- Specific cases where it breaks (holes and more).
- The scan of antipodal points is by calculating the distance
between two points. Could be made faster by calculating the
height of a triangle instead?
*/
std::pair < double, double > getMinMaxFeretDiameters(std::vector < Point > & convexHull)
{
	// Resulting minimum and maximum Feret diameters
	std::pair < double, double > minMaxFeretDiameters;

	// In a convex polygon the antipodal points for vertex p_i
	// are guaranteed to be between the vertex farthest from the
	// edge p_i-1 -- p_i and the vertex farthest from p_i -- 
	// p_i+1 when transversing the boundary of the convex polygon.
	std::vector < std::pair < Point, Point > > antiPodalPairs;

	// Calculate area from triangle formed from 3 vertices.
	// Compare the area of first and second triangle. 
	long n = (long)convexHull.size(); // Number of vertices.
	long k = 1; // Top vertex of triangle.
	long m = n - 1; // Left vertex of triangle.
	long i = 0; // Right vertex of triangle.

	if (n == 0) // Case of no input.
	{
		minMaxFeretDiameters.first = 0;
		minMaxFeretDiameters.second = 0;

		return minMaxFeretDiameters;
	}
	else if (n == 2) // Case of only two points.
	{
		minMaxFeretDiameters.first = 1;
		minMaxFeretDiameters.second = distance(convexHull[0],
			convexHull[1]);

		return minMaxFeretDiameters;
	}
	else if (n >= 3) // Minimum is a triangle...
	{

		// Find furthest point, Pk, from edge Pm-P1-P2. 
		while (areaTriangle(convexHull[m], convexHull[i],
			convexHull[k + 1]) > areaTriangle(convexHull[m], convexHull[i],
				convexHull[k]))
		{
			k++; // Jump to next vertice.
		}

		// From Pk find the next antipodal points by transversing the 
		// boundary until the area formed by the triangle Pi-Pi+1-Pj+1
		// is greater than the area formed by the triangle Pi-Pi+1-Pj, 
		// which is the vertex, Pl, furthest from the edge Pi-Pi+1.
		i = 0;
		long j = k;

		// Advance to the next edge until Pi-Pk (i<=k) and the antipodal 
		// Pk, Pm is reached (j<=m).
		while (i <= k && j <= m)
		{
			// First antipodal point
			antiPodalPairs.push_back(std::make_pair(Point(convexHull[i]),
				Point(convexHull[j])));

			while (j < m && areaTriangle(convexHull[i], convexHull[i + 1],
				convexHull[j + 1]) > areaTriangle(convexHull[i],
					convexHull[i + 1], convexHull[j]))
			{
				antiPodalPairs.push_back(std::make_pair(Point(convexHull[i]),
					Point(convexHull[j + 1])));
				j++; // Jump to next vertex
			}
			i++;
		}

		// Scan for maximum and minimum distances.
		double maxDistance = 0.0;
		double minDistance = 9999; // float_max instead?

		for (int i = 0; i < antiPodalPairs.size(); i++)
		{
			double dist = distance(antiPodalPairs[i].first,
				antiPodalPairs[i].second);

			//std::cout << dist << std::endl;
			if (dist > maxDistance)
			{
				maxDistance = dist;
			}
			if (dist < minDistance)
			{
				minDistance = dist;
			}
		}

		minMaxFeretDiameters.first = minDistance;
		minMaxFeretDiameters.second = maxDistance;
	}

	return minMaxFeretDiameters;
}

DESCRIPTORS void feretDiameters(int *xPixels, int *yPixels, int nPixels,
	double *minFeret, double *maxFeret)
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

	// Calculate minimum and maximum Feret diameters.
	std::pair <double, double> minMaxFeretDiameters;
	minMaxFeretDiameters = getMinMaxFeretDiameters(points);
	*minFeret = minMaxFeretDiameters.first;
	*maxFeret = minMaxFeretDiameters.second;
}

/*
Skeletonize.
*/
DESCRIPTORS void getSkeleton(
	int *x, 
	int *y,
	int *nPixels,
	int *nRows, 
	int *nCols,
	int **resultX,
	int **resultY,
	int *skelSize)
{
	// Zero matrix.
	std::vector< std::vector< unsigned int > >
		skeleton(*nRows, std::vector< unsigned int >(*nCols, 0));
	
	// Convert to array.
	std::vector< std::vector< unsigned int > >
		particleMatrix(*nRows, std::vector< unsigned int >(*nCols, 0));
	int size = *nRows * *nCols;
	
	for (int i = 0; i < *nPixels; i++)
	{
		//std::cout << i << std::endl;
		//std::cout << x[i] << ", " << y[i] << std::endl;
		particleMatrix[ x[i] ][ y[i] ] = 1;
	}
	
	skeletonize(particleMatrix, skeleton, *nRows, *nCols);

	// Convert vector to array.
	int k = 0;
	for (int i = 0; i < skeleton.size(); i++)
	{
		for (int j = 0; i < skeleton[0].size(); j++)
		{
			std::cout << i << ", " << j << std::endl;
			if (skeleton[i][j] == 1)
			{
				k++;
			}
		}
	}
	*skelSize = k;
	
	*resultX = new int[*skelSize];
	*resultY = new int[*skelSize];

	k = 0;
	for (int i = 0; i < skeleton.size(); i++)
	{
		for (int j = 0; i < skeleton[0].size(); j++)
		{
			if (skeleton[i][j] == 1)
			{
				(*resultX)[k] = i;
				(*resultY)[k] = j;
				k++;
			}
		}
	}
}

DESCRIPTORS void free_mem_float(float** a)
{
	delete[] a;
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
