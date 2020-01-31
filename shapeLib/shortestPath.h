/*
Declarations of functions for shortestPath.cpp.

Author: Rasmus Vest Nielsen
*/

#include <iostream>
#include <vector>

// Set guard
#ifndef SHORTESTPATH_H
#define SHORTESTPATH_H

// Function declaration.
bool isValid(int row, int col, int nRows, int nCols);
int breadthFirstSearch(std::vector< std::vector< unsigned int > > & skeleton,
	Point src, Point dest);

#endif // End of guard