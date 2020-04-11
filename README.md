# shapeLib
Library for analysing particle size and shape from sIMPe.

Calculates the following shape discriptors:
* Convexity, Area, Perimeter length, Maximum
	and minimum Feret Diameters, Solidity /
	Circularity, Aspect ratio, Circular equivalent diameter,
	Curl / Straightness, Fiber average diameter,
	Fiber elongation, Solidity

Example of Python 3.8 interface to shapeLib.dll.

shapeLib.dll is written in C-idiomatic style C++ code, making it "easily" 
integrated with Python using Ctypes.

All functions return Python- or numpy-like variables and arrays (handled by 
hapeLibPyUtils.py).

NOTES:

* However, since Pythonic created pointers do not run out of scope, memory 
  leaks are to be expected, and the PC should be restarted after this script 
  has been run (this should be easy to fix).

* Dynamic C-arrays not running out of scope in shapeLib.dll results in memory
  leakage on the library side. This is fixed by using C++ vectors instead 
  (should be completed, but have to be checked).

* Beware of dynamic memmory allocations in shapeLib.dll needs to be deallocated 
  using the function free_mem(arr)! However only for variables in the DLL, 
  which are NOT to be available on the Python side. This is an artifact of 
  previous bad decisions and will be fixed.

PROBLEMS:

* Boundary algorithm cannot handle hollow particles, and I do not want to 
  implement that functionality.

* Calculations of fiber curl, length, average diameter, or elongation fails for 
  some particles. Need to figure out why.
 
