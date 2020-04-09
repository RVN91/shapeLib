"""
Utility functions for shapeLib.py

Author: Rasmus Vest Nielsen
"""

import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from ctypes import *
import ctypes

# Load library. Full path must be specified.
dll_c = cdll.LoadLibrary(("C:/Users/rasmus/Desktop/microplastic/stuff/"
                         "shapeLib/Debug/shapeLib.dll"))

# "Global" variables...
# If debug is enabled, file names are set to default. If disabled, user is 
# asked for file names.
DEBUG = True
BYTE_POSITION = 0 # Keeps track of cursor in binary file.

# Log...
def log(var):
    """
    Prints log messages if DEBUG is enabled.
    """
    if DEBUG:
            print(var)

# Create result file. 
dll_c.createResultFile.argtypes = [ctypes.c_bool]
dll_c.createResultFile.restype = ctypes.c_char_p
inputFileName = dll_c.createResultFile(DEBUG)

def load_particle():
    """
    Loads a particle from a sIMPLe binary file
    and streaming it to a numpy array.
    """
    global BYTE_POSITION

    # Get the number of particles.
    dll_c.getParticleCount.restype = ctypes.c_int
    byte_position = ctypes.c_long(BYTE_POSITION) 

    # The number of particles is stored at byte_position = 0.
    # Any other input will fail!
    n_particles = dll_c.getParticleCount(inputFileName, 
            ctypes.byref(byte_position))

    # Update byte position.
    BYTE_POSITION = byte_position.value

    return n_particles
    
def read_particle_coords():
    """
    Reads the pixel count of the current particle and then creates a numpy 
    array for shapeLib.dll to place the pixel coordinates.

    Allocates a numpy array, initializes a pointer to numpy array and 
    associates the point to an object.
    
    This way we do not need to manage memory manually(?).
    """
    global BYTE_POSITION
    
    byte_position = ctypes.c_long(BYTE_POSITION) 
    
    n_pixels = ctypes.c_int(0)
    dll_c.getPixelsCount(inputFileName, 
            ctypes.byref(byte_position), 
            ctypes.byref(n_pixels))
    
    n_pixels = n_pixels.value

    # Update byte position.
    byte_position = ctypes.c_long(byte_position.value)

    # Create empty numpy array.
    x_pixels = np.zeros(n_pixels, dtype = ctypes.c_int)
    y_pixels = np.zeros(n_pixels, dtype = ctypes.c_int)

    # Create pointer for empty array object. This lets Python handle memory 
    # allocation and therefore we do not need to deallocate memory manually.
    ptr_x_pixels = ctypes.POINTER(ctypes.c_int)
    ptr_y_pixels = ctypes.POINTER(ctypes.c_int)

    # Create object struct with pointers to pixel arrays.
    obj_x_pixels = x_pixels.ctypes.data_as(ptr_x_pixels)
    obj_y_pixels = y_pixels.ctypes.data_as(ptr_y_pixels)
    
    # Create pointer for minor and major dimensions.
    minor_dim = ctypes.c_float(0.0)
    major_dim = ctypes.c_float(0.0)

    # Get x and y pixel coordinates into numpy arrays.
    n_pixels = ctypes.c_int(n_pixels)
    dll_c.getParticle(inputFileName, 
            ctypes.byref(byte_position), 
            ctypes.byref(n_pixels),
            obj_x_pixels,
            obj_y_pixels,
            ctypes.byref(minor_dim),
            ctypes.byref(major_dim))

    # Change dtype to numpy's version of long.
    x_pixels = x_pixels.astype(np.int_)
    y_pixels = y_pixels.astype(np.int_)

    # Update byte position.
    BYTE_POSITION = byte_position.value
    
    return x_pixels, y_pixels, minor_dim.value, major_dim.value

def dynamic_array():
    """
    Initializes a dynamic array as a pointers to pointers of vectors native to 
    the C language.

    Returns a pointer to a dynamic and a pointer to the unspecified size.
    """
    # Create pointer to pointer.
    dll_c.testArray.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),
            ctypes.POINTER(ctypes.c_int)]

    mem = ctypes.POINTER(ctypes.c_int)()
    size = ctypes.c_int(0)

    dll_c.testArray(ctypes.byref(mem), ctypes.byref(size))

    return mem, size

def offset_particle_coords(x_pixels, y_pixels):
    """
    Offsets particle pixels from sIMPLe to zero axis reference. Makes plotting 
    of figures easier.
    
    Pads particle pixels by +1 in x and y directions
    """
    # Calculate particle rows and columns and pad by +3.
    # The offset is needed for the 8-directional boundary 
    # tracer, where pixels on the boundary of the matrix
    # are not allowed.
    n_rows = max(y_pixels) - min(y_pixels) + 3  
    n_cols = max(x_pixels) - min(x_pixels) + 3

    # Offset particle pixels by -min(x) and -min(y) to 
    # avoid large arrays. This sets the particle to be
    # atleast one pixel away from the border.
    x_offset = x_pixels - min(x_pixels) + 1 
    y_offset = y_pixels - min(y_pixels) + 1

    return x_offset, y_offset, n_rows, n_cols

def generate_particle_matrix(x, y, n_rows, n_cols):
    """
    Generates a matrix from x and y coordinates for plotting purposes.
    """
    # Insert into matrix.
    matrix = np.zeros((n_rows, n_cols))
    
    # Reverse x and y? Error in C++ code...
    for i, j in zip(y, x):
        matrix[j, i] = 1

    return matrix

def generate_1d_array(matrix):
    """
    Generates 1D array from padded matrix.

    Used for the boundary tracer.
    """ 
    # Generate 1D array.
    array = matrix.flatten()

    return array

def boundary_points(x, y, n_cols, n_rows):
    """
    Initializes a dynamic array as a pointers to pointers of vectors native to 
    the C language.

    Returns a pointer to a dynamic and a pointer to the unspecified size.
    """
    # Generate particle vector. Needed for boundary tracer.
    matrix = generate_particle_matrix(x, y, n_rows, n_cols)
    array  = generate_1d_array(matrix)
    
    # Convert imput data to native C format.
    n_rows_c = ctypes.c_int(n_rows)
    n_cols_c = ctypes.c_int(n_cols)
    
    # Convert numpy particle pixels array to C style
    # 32 bit "int".
    c_int_ptr = ctypes.POINTER(ctypes.c_int)   # Pointer.
    data      = array.astype(np.int32)         # np array to 32 bit ints.
    data_ptr  = data.ctypes.data_as(c_int_ptr) # np array to Ctype 32 bit ints.

    n_pixels = ctypes.c_int(array.shape[0])
    
    boundary_points_x_ptr = ctypes.POINTER(ctypes.c_int)()
    boundary_points_y_ptr = ctypes.POINTER(ctypes.c_int)()
    
    # Dynamic size.
    size = ctypes.c_int(0)

    dll_c.traceBoundary(data_ptr, ctypes.byref(n_pixels), 
            n_rows_c, n_cols_c, 
            ctypes.byref(boundary_points_x_ptr), 
            ctypes.byref(boundary_points_y_ptr), 
            ctypes.byref(size))
    
    # Convert to numpy array.
    boundary_points_x = np.empty(size.value, dtype = np.intc)
    boundary_points_y = np.empty(size.value, dtype = np.intc)
    
    for i in range(0, size.value):
        boundary_points_x[i] = boundary_points_x_ptr[i]        
        boundary_points_y[i] = boundary_points_y_ptr[i]        

    return boundary_points_x, boundary_points_y, size

def convex_hull(x, y):        
    """
    Calculates the convex hull.

    NOTE: The input coordinates MUST be the raw particle coordinates and not 
          padded!
    """
    # Calculate convex hull.
    n_pixels = ctypes.c_int(x.shape[0])
    
    # Set argument types.
    dll_c.convexHull.argtypes = [ndpointer(ctypes.c_int, 
        flags="C_CONTIGUOUS"), ndpointer(ctypes.c_int, 
        flags="C_CONTIGUOUS"), ctypes.c_int, 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int)), 
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int)), 
        ctypes.POINTER(ctypes.c_int)]
    
    # Set pointers.
    convex_points_x_ptr = ctypes.POINTER(ctypes.c_int)()
    convex_points_y_ptr = ctypes.POINTER(ctypes.c_int)()
    size = ctypes.c_int(0)
    
    # Calculate convex hull.
    dll_c.convexHull(x, y, n_pixels, 
            ctypes.byref(convex_points_x_ptr), 
            ctypes.byref(convex_points_y_ptr), 
            ctypes.byref(size))

    # Convert to numpy array.
    convex_points_x = np.empty(size.value, dtype = np.intc)
    convex_points_y = np.empty(size.value, dtype = np.intc)
    
    for i in range(0, size.value):
        convex_points_x[i] = convex_points_x_ptr[i]        
        convex_points_y[i] = convex_points_y_ptr[i]

    return convex_points_x, convex_points_y

def calc_perimeter(x, y):
    """
    Calculates the perimeter enclosed by a given set of points.
    """
    # Calculate the perimeter.
    n_points  = ctypes.c_int(x.shape[0])
    perimeter = ctypes.c_double(0.0)
    
    # Set argument types.
    dll_c.findPerimeter.argtypes = [ 
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double)]
    
    # Calculate convex hull.
    dll_c.findPerimeter(x, y, ctypes.byref(n_points), ctypes.byref(perimeter))
    
    return perimeter.value

def calc_area(x, y):
    """
    Calculates the area as enclosed by a given set of points.
    """
    # converting pythonic types to c types.
    n_points  = ctypes.c_int(x.shape[0])
    area      = ctypes.c_double(0.0)
    
    # set argument types.
    dll_c.calculateArea.argtypes = [ 
        ndpointer(ctypes.c_int, flags="c_contiguous"),
        ndpointer(ctypes.c_int, flags="c_contiguous"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double)]
    
    # calculate convex hull.
    dll_c.calculateArea(x, y, 
            ctypes.byref(n_points), 
            ctypes.byref(area))
    
    return area.value

def calc_solidity(area, convex_area):
    """
    Calculates the solidity as: solidity = area / convex_area
    """
    # Converting Pythonic types to C types.
    area = ctypes.c_double(area)
    convex_area = ctypes.c_double(convex_area)
    solidity = ctypes.c_double(0.0)
    
    dll_c.calculateSolidity(
            ctypes.byref(area),
            ctypes.byref(convex_area),
            ctypes.byref(solidity))

    return solidity.value

def calc_feret_diameter(x, y):
    """
    Calculates the minimum and maximum Feret diameters.
    """
    # converting pythonic types to c types.
    n_points  = ctypes.c_int(x.shape[0])
    min_feret = ctypes.c_double(0.0)
    max_feret = ctypes.c_double(0.0)
    
    # Set argument types.
    dll_c.feretDiameters.argtypes = [ 
        ndpointer(ctypes.c_int, flags="c_contiguous"),
        ndpointer(ctypes.c_int, flags="c_contiguous"),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double)]
    
    # calculate convex hull.
    dll_c.feretDiameters(x, y, n_points, 
            ctypes.byref(min_feret),
            ctypes.byref(max_feret))
    
    return min_feret.value, max_feret.value

def skeletonize(x, y, n_cols, n_rows):
    """
    Initializes a dynamic array as a pointers to pointers of vectors native to 
    the C language.

    Returns a pointer to a dynamic and a pointer to the unspecified size.
    """ 
    # Set argument types.
    dll_c.getSkeleton.argtypes = [ 
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),
        ctypes.POINTER(ctypes.c_int)]
    
    # Convert imput data to native C format.
    n_pixels = ctypes.c_int(x.shape[0])
    
    n_rows_c = ctypes.c_int(n_rows)
    n_cols_c = ctypes.c_int(n_cols)
    
    skeleton_x_ptr = ctypes.POINTER(ctypes.c_int)()
    skeleton_y_ptr = ctypes.POINTER(ctypes.c_int)()
    
    # Dynamic size for skeleton points.
    size = ctypes.c_int(0)

    dll_c.getSkeleton(
            x, y, 
            ctypes.byref(n_pixels), 
            ctypes.byref(n_rows_c), 
            ctypes.byref(n_cols_c), 
            ctypes.byref(skeleton_x_ptr), 
            ctypes.byref(skeleton_y_ptr), 
            ctypes.byref(size))
    
    # Convert to numpy array.
    skeleton_points_x = np.empty(size.value, dtype = np.intc)
    skeleton_points_y = np.empty(size.value, dtype = np.intc)
    
    for i in range(0, size.value):
        skeleton_points_x[i] = skeleton_x_ptr[i]        
        skeleton_points_y[i] = skeleton_y_ptr[i]        

    return skeleton_points_x, skeleton_points_y, size

def get_end_points(x, y, n_cols, n_rows):
    """
    Finds the end-points of the skeleton matrix and returns the matrix indeces.
    """

    # Set argument types.
    dll_c.getEndPoints.argtypes = [ 
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),
        ctypes.POINTER(ctypes.c_int)]

    # Convert imput data to native C format.
    n_pixels = ctypes.c_int(x.shape[0])
    
    n_rows_c = ctypes.c_int(n_rows)
    n_cols_c = ctypes.c_int(n_cols)
    
    end_points_x_ptr = ctypes.POINTER(ctypes.c_int)()
    end_points_y_ptr = ctypes.POINTER(ctypes.c_int)()
    
    # Dynamic size for end points.
    size = ctypes.c_int(0)
    
    dll_c.getEndPoints(
            x, y, 
            ctypes.byref(n_pixels), 
            ctypes.byref(n_rows_c), 
            ctypes.byref(n_cols_c), 
            ctypes.byref(end_points_x_ptr), 
            ctypes.byref(end_points_y_ptr), 
            ctypes.byref(size))
    
    # Convert to numpy array.
    end_points_x = np.empty(size.value, dtype = np.intc)
    end_points_y = np.empty(size.value, dtype = np.intc)
    
    for i in range(0, size.value):
        end_points_x[i] = end_points_x_ptr[i]        
        end_points_y[i] = end_points_y_ptr[i]        

    return end_points_x, end_points_y, size

def get_fiber_length(x_end_points, y_end_points, 
        x_skeleton, y_skeleton, n_rows, n_cols):
    """
    Finds the maximum distance between the end-points of the particle skeleton.
    """

    # Set argument types.
    dll_c.getMaxDistanceEndPoints.argtypes = [ 
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)]

    # Convert imput data to native C format.
    max_distance = ctypes.c_int(0);
    n_end_points_c = ctypes.c_int(x_end_points.shape[0])
    n_rows_c = ctypes.c_int(n_cols) # Why do I need to    
    n_cols_c = ctypes.c_int(n_rows) # switch rows and columns?
    n_skeleton_points_c = ctypes.c_int(x_skeleton.shape[0])
   
    """
    print(n_rows)
    print(n_cols)
    print(x_skeleton)
    print(y_skeleton)
    """

    dll_c.getMaxDistanceEndPoints(
            x_end_points, y_end_points, 
            ctypes.byref(n_end_points_c), 
            ctypes.byref(max_distance),
            x_skeleton, y_skeleton,
            ctypes.byref(n_rows_c), 
            ctypes.byref(n_cols_c),
            ctypes.byref(n_skeleton_points_c))

    return max_distance.value

