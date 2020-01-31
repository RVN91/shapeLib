"""
Use example of Python 3.8 interface to shapeLib.dll.

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
  (should be done, but have to be checked).

* Beware of dynamic memmory allocations in shapeLib.dll needs to be deallocated 
  using the function free_mem(arr)! However only for variables in the DLL, 
  which are NOT to be available on the Python side.

Author: Rasmus Vest Nielsen
"""

import numpy as np
import pyUtils
import matplotlib.pyplot as plt

# main().
def main():
    """
    Main function. Run functions from here.
    """
    # Load particles
    n_particles = pyUtils.load_particle()
    
    # Read x, y particle coordinates. 
    x_coords, y_coords, minor_dim, major_dim = pyUtils.read_coords()
   
    # Zero reference the coordinates, add +1 padding.
    x_offset, y_offset, n_rows, n_cols = pyUtils.offset_coords(x_coords, y_coords) 
    
    # Particle size threshold. I have arbitrarely chosen 5 pixels as the 
    # minimum pixel size. This choice is mainly based on what I think makes
    # sense to include in shape analysis. It makes sense for particle size. Let 
    # me know what you think about that. Secondly, to calulate convex hull and 
    # Feret Diameters, a minimum particle size of 3 pixels is required. 
    particle_size_threshold = 5 # [pixels]
    
    if x_coords.shape[0] > particle_size_threshold:
        
        # Trace boundary contour.
        boundary_x, boundary_y, size = pyUtils.boundary_points(x_offset, 
                y_offset, n_rows, n_cols)
        
        # Calculate convex hull.
        convex_x, convex_y = pyUtils.convex_hull(x_offset, 
                y_offset)
        
        # Calculate the particle perimenter.
        perimeter = pyUtils.calc_perimeter(boundary_x, 
                boundary_y)

        convex_perimeter = pyUtils.calc_perimeter(convex_x, 
                convex_y);
       
        # Calculate the area.
        area        = pyUtils.calc_area(boundary_x, boundary_y)
        convex_area = pyUtils.calc_area(convex_x, convex_y)
        
        # Calculate solidity.
        solidity = area /convex_area 
        
        # Calculate convexity.
        convexity = convex_perimeter / perimeter
        
        # Calculate equivalent circular diameter.
        pi = 3.14159265358979323846
        equivalent_diameter = np.sqrt(4 * area * pi) / convex_perimeter**2

        # Calculate circularity and compactness.
        circularity = (4 * pi * area) / convex_perimeter**2
        
        compactness = (4 * pi * area) / perimeter**2
        
        # Calculate minumim and maximum Feret diameters.
        min_feret, max_feret = pyUtils.calc_feret_diameter(convex_x, 
                convex_y)
        
        # Aspect ratio.
        aspect_ratio = major_dim / minor_dim

        # Feret aspect ratio.
        feret_ratio = min_feret / max_feret
        
        # Fiber-related parameters.
        # Get skeleton.
        skeleton_x, skeleton_y, n_points_skeleton = pyUtils.skeletonize(
                x_offset, y_offset, 
                n_rows, n_cols)

        # Find end points. 
        end_points_x, end_points_y, n_end_points = pyUtils.get_end_points(
                skeleton_x, skeleton_y, 
                n_rows, n_cols)
       
        # Find the fiber length.
        fiber_length = pyUtils.get_fiber_length(end_points_x, end_points_y,
                skeleton_x, skeleton_y, n_rows, n_cols)

        # Calculate curl or straigthness (Max. Feret diameter 
        # fiber length).
        # 0 <= x <= 1; where 1 is a straight shape.
        # A measure of the curl up a fiberlike shape.
        curl = max_feret / fiber_length;

        # Calculate average diameter (fiberlike shapes).
        # area / Lenght = width (assuming the fiber is
        # straightened out to a rectangle).
        fiber_average_diameter = area / fiber_length;

        # Calculate fiber elongation (average diameter /
        # fiber length).
        # 0 <= x <= 1; where => 0 the particle is more elongated.
        # A potential measure of the fiberlike characteristics.
        fiber_elongation = fiber_average_diameter / fiber_length;

        print(perimeter, convex_perimeter, area, convex_area, 
                solidity, convexity, equivalent_diameter, circularity, 
                compactness, min_feret, max_feret, minor_dim, major_dim,
                aspect_ratio, feret_ratio, fiber_length, curl, 
                fiber_average_diameter, fiber_elongation)

        # Plot results...
        # Zero reference the coordinates, add +1 padding.
        convex_x, convex_y, n_convex_rows, n_convex_cols = pyUtils.offset_coords(convex_x, convex_y) 
    
        # Insert into matrix.
        particle_matrix = pyUtils.generate_particle_matrix(y_offset, x_offset,
                n_rows, n_cols)
        convex_matrix   = pyUtils.generate_particle_matrix(convex_y, convex_x, 
                n_rows, n_cols)
        skeleton_matrix = pyUtils.generate_particle_matrix(skeleton_y, 
                skeleton_x, n_rows, n_cols)
        end_points_matrix = pyUtils.generate_particle_matrix(end_points_y, 
                end_points_x, n_rows, n_cols)

        # Plot figures.
        fig, axes = plt.subplots(2, 2)
        
        # Particle shape.
        axes[0, 0].imshow(particle_matrix)
        axes[0, 0].set_title("Particle shape.")
        
        # Convex particle shape.
        axes[0, 1].scatter(x_offset, y_offset, c = "y", marker = "s")
        axes[0, 1].plot(convex_x, convex_y, c = "r", marker = "s")
        axes[0, 1].set_ylim(axes[0, 1].get_ylim()[::-1]) # Reverse yaxis.
        axes[0, 1].set_title("Convex hull.")
        
        # Particle skeleton.
        axes[1, 0].imshow(skeleton_matrix)
        axes[1, 0].set_title("Skeleton.")
        
        # End points.
        axes[1, 1].imshow(skeleton_matrix)
        axes[1, 1].imshow(end_points_matrix, alpha = 0.4)
        axes[1, 1].set_title("End Points.")
        plt.show()

if __name__ == "__main__": 
    main()
