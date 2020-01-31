/*
	Find the largest path in a binary matrix given start and
	end points using first breadh-first search.
	
	https://www.geeksforgeeks.org/shortest-path-in-a-binary-maze/
	
	Minor rewrite:
	- Extended the code to handle special cases.
	- Made it easier to understand and maintain.

	Author: Rasmus Vest Nielsen
*/

#include "stdafx.h"

#include "convexHull.h"
#include <iostream>
#include <queue> 
#include <vector>
#include "shortestPath.h"

using namespace std;

// A Data Structure for queue used in BFS. 
struct queueNode
{
	Point pt;  // The cordinates of a cell 
	int dist;  // cell's distance of from the source.
};

// Check whether given cell (row, col) is a valid 
// cell or not. 
bool isValid(int row, int col, int nRows, int nCols)
{
	// Return true if row number and column number 
	// is in range.
	return (row >= 0) && (row < nRows) &&
		(col >= 0) && (col < nCols);
}

// Function to find the shortest path between 
// a given source cell to a destination cell. 
int breadthFirstSearch(std::vector< std::vector< unsigned int > > & inMat, Point src, Point dest)
{
	// Check source and destination cell 
	// of the matrix have value 1.
	if (!inMat[src.x][src.y] || !inMat[dest.x][dest.y])
		return -1;
	
	// These arrays are used to get row and column 
	// numbers of 8 neighbours of a given cell.
	int rowNum[] = { -1, -1, -1,  0, 0,  1, 1, 1 };
	int colNum[] = { -1,  0,  1, -1, 1, -1, 0, 1 };

	int nRows = inMat.size();
	int nCols = inMat[0].size();

	// Bolean binary array. 0 = false; 1 = true.
	std::vector< std::vector< unsigned int > >
		visited(nRows, std::vector< unsigned int >(nCols, 0));

	// Mark the source cell as visited.
	visited[src.x][src.y] = 1;

	// Create a queue for BFS.
	queue<queueNode> q;

	// Distance of source cell is 0. 
	queueNode s = { src, 0 };
	q.push(s);  // Enqueue source cell.

	// Do a BFS starting from source cell.
	while (!q.empty())
	{
		queueNode curr = q.front();
		Point pt = curr.pt;

		// If we have reached the destination cell, 
		// we are done. 
		if (pt.x == dest.x && pt.y == dest.y)
			return curr.dist;

		// Otherwise dequeue the front cell in the queue 
		// and enqueue its adjacent cells.
		q.pop();

		/*
		// Print current cell
		std::cout << "Current cell: " << pt.x << ", " << pt.y <<
			", value: " << inMat[pt.x][pt.y] << std::endl;
		*/

		for (int i = 0; i < 8; i++)
		{
			int row = pt.x + rowNum[i];
			int col = pt.y + colNum[i];

			/*
			// Print adjacent cell
			int value = inMat[row][col];
			cout << "Adjacent cell: " << row << ", " << col << 
			", value: " << value << endl;
			*/

			// If adjacent cell is valid, has path and 
			// not visited yet, enqueue it. 
			if (isValid(row, col, nRows, nCols) && inMat[row][col] 
				== 1 && !visited[row][col])
			{
				// Mark cell as visited and enqueue it.
				visited[row][col] = true;
				queueNode Adjcell = { { row, col },
					curr.dist + 1 };
				q.push(Adjcell);
			}
		}
	}

	// Return -1 if destination cannot be reached. 
	return -1;
}