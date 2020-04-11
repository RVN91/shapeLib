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
  (should be completed, but have to be checked).

* Beware of dynamic memmory allocations in shapeLib.dll needs to be deallocated 
  using the function free_mem(arr)! However only for variables in the DLL, 
  which are NOT to be available on the Python side. This is an artifact of 
  previous bad decisions and will be fixed.

Author: Rasmus Vest Nielsen
"""

import numpy as np
import pyUtils
import pandas as pd
import os
import matplotlib.pyplot as plt

# File names.
INPUT_FILE_NAME = "T1_1.pmp"
OUTPUT_FILE_NAME = "results.csv"

# Extra options.
PLOT_FIGS = False  # If true, plot particle shape in "/plots".
PLOT_PARTICLE_N = 56 # Plots particle shape of particle N. 
DEBUG = True      # If true, extra information is available in the terminal.

def main():
    """
    Computes particle shapes and sizes, and prints them to a CSV file.
    
    Run functions from here.
    """
    # Print enabled options.
    print_options()

    # Load particles
    n_particles = pyUtils.load_particle(INPUT_FILE_NAME)
    
    # Iterate through each particle and store particle size and shape in Data 
    # Frame.
    df_list = []
    for n_particle in range(0, n_particles): 
        log("Particle {0} of {1}...".format(n_particle, n_particles))
        
        df_tmp, particle_coords = get_particle_shape_size(n_particle)
        df_list.append(df_tmp)
        
        if PLOT_FIGS == True and df_tmp["area"].to_list() != [0]:
            plot_shape(n_particle, particle_coords)

        if PLOT_PARTICLE_N == n_particle:
            plot_shape(n_particle, particle_coords)

    # Concatenate Data Frames.
    df = pd.concat(df_list)
    df.index.name = "particle_no"
    df.to_csv("shape_sizes.csv")
    print(df)

def get_particle_shape_size(n_particle):
    """
    Reads particles and calculates particle shape and size. 
    
    Parameters:
        None.

    Returns:
        df : Pandas Data Frame 
            Particle size and shape descriptors.
    """
    # Read x, y particle coordinates. 
    x_coords, y_coords, minor_dim, major_dim = pyUtils.read_particle_coords(INPUT_FILE_NAME)
   
    # Zero reference the coordinates, add +1 padding.
    x_offset, y_offset, n_rows, n_cols = \
        pyUtils.offset_particle_coords(x_coords, y_coords) 
    
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
        equivalent_diameter = np.sqrt(4 * area / pi)

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

        try:
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

        except:
            log("Ups. Something went wrong...")
            log(max_feret)
            log(fiber_length)
            curl = 0
            fiber_average_diameter = 0
            fiber_elongation = 0


        # Store data in CSV file.
        d = {
            "perimeter": perimeter, 
            "convex_perimeter": convex_perimeter, 
            "area": area, 
            "convex_area": convex_area, 
            "solidity": solidity, 
            "convexity": convexity, 
            "equivalent_diameter": equivalent_diameter, 
            "circularity": circularity, 
            "compactness": compactness, 
            "min_feret": min_feret, 
            "max_feret": max_feret, 
            "minor_dim": minor_dim, 
            "major_dim": major_dim, 
            "aspect_ratio": aspect_ratio, 
            "feret_ratio": feret_ratio, 
            "fiber_length": fiber_length, 
            "curl": curl, 
            "fiber_average_diameter": fiber_average_diameter, 
            "fiber_elongation": fiber_elongation
        }
        
        df = pd.DataFrame(data = d, index = [n_particle]) 
        
        if n_particle == 179:
            print(df)

        #print(df)
        #print(perimeter, convex_perimeter, area, convex_area, 
        #        solidity, convexity, equivalent_diameter, circularity, 
        #        compactness, min_feret, max_feret, minor_dim, major_dim,
        #        aspect_ratio, feret_ratio, fiber_length, curl, 
        #        fiber_average_diameter, fiber_elongation)
        
        if PLOT_FIGS:
            particle_coords = {"convex_x": convex_x,
                               "convex_y": convex_y,
                               "x_offset": x_offset,
                               "y_offset": y_offset,
                               "skeleton_x": skeleton_x,
                               "skeleton_y": skeleton_y,
                               "end_points_x": end_points_x,
                               "end_points_y": end_points_y,
                               "n_rows": n_rows,
                               "n_cols": n_cols,
                               "boundary_x": boundary_x,
                               "boundary_y": boundary_y
                    }
        else:
            particle_coords = {"convex_x": convex_x,
                               "convex_y": convex_y,
                               "x_offset": x_offset,
                               "y_offset": y_offset,
                               "skeleton_x": skeleton_x,
                               "skeleton_y": skeleton_y,
                               "end_points_x": end_points_x,
                               "end_points_y": end_points_y,
                               "n_rows": n_rows,
                               "n_cols": n_cols,
                               "boundary_x": boundary_x,
                               "boundary_y": boundary_y
                    }
    else: # Not enough pixels.
        d = {
            "perimeter": 0, 
            "convex_perimeter": 0, 
            "area": 0, 
            "convex_area": 0, 
            "solidity": 0, 
            "convexity": 0, 
            "equivalent_diameter": 0, 
            "circularity": 0, 
            "compactness": 0, 
            "min_feret": 0, 
            "max_feret": 0, 
            "minor_dim": 0, 
            "major_dim": 0, 
            "aspect_ratio": 0, 
            "feret_ratio": 0, 
            "fiber_length": 0, 
            "curl": 0, 
            "fiber_average_diameter": 0, 
            "fiber_elongation": 0
        }
        
        df = pd.DataFrame(data = d, index = [n_particle]) 
        
        particle_coords = {"convex_x": 0,
                           "convex_y": 0,
                           "x_offset": 0,
                           "y_offset": 0,
                           "skeleton_x": 0,
                           "skeleton_y": 0,
                           "end_points_x": 0,
                           "end_points_y": 0,
                           "n_rows": 0,
                           "n_cols": 0,
                           "boundary_x": 0,
                           "boundary_y": 0
        }

    return df, particle_coords 


def plot_shape(n_particle, coords):
    """
    Plot particle shape and saves it in subdirectory "/plots".

    Parametes:
        n_particle : int
        df:          pandas Data Frame

    Returns:
        None.
    """
    # Zero reference the coordinates, add +1 padding.
    #boundary_x, boundary_y, n_boundary_rows, n_boundary_cols = \
    #        pyUtils.offset_particle_coords(coords["boundary_x"], coords["boundary_y"]) 
    # Insert into matrix.
    particle_matrix = pyUtils.generate_particle_matrix(coords["y_offset"], 
            coords["x_offset"], coords["n_rows"], coords["n_cols"])
    boundary_matrix = pyUtils.generate_particle_matrix(coords["boundary_y"], 
            coords["boundary_x"], coords["n_rows"], coords["n_cols"])
    skeleton_matrix = pyUtils.generate_particle_matrix(coords["skeleton_y"], 
            coords["skeleton_x"], coords["n_rows"], coords["n_cols"])
    end_points_matrix = pyUtils.generate_particle_matrix(coords["end_points_y"], 
            coords["end_points_x"], coords["n_rows"], coords["n_cols"])

    # Plot figures.
    fig, axes = plt.subplots(2, 2)
    
    # Particle shape.
    axes[0, 0].imshow(particle_matrix + boundary_matrix)
    axes[0, 0].set_title("Particle shape.")
    
    # Convex particle shape.
    axes[0, 1].scatter(coords["x_offset"], coords["y_offset"], 
            c = "y", marker = "s")
    axes[0, 1].plot(coords["convex_x"], coords["convex_y"], 
            c = "r", marker = "s")
    axes[0, 1].set_ylim(axes[0, 1].get_ylim()[::-1]) # Reverse yaxis.
    axes[0, 1].set_title("Convex hull.")
    
    # Particle skeleton.
    axes[1, 0].imshow(skeleton_matrix)
    axes[1, 0].set_title("Skeleton.")
    
    # End points.
    axes[1, 1].imshow(skeleton_matrix + end_points_matrix)
    axes[1, 1].set_title("End Points.")
    
    # Adjust figure.
    fig.subplots_adjust(
        left = 0.04,  
        right = 0.99,   
        bottom = 0.06,              
        top = 0.95,    
        wspace = 0.16,  
        hspace = 0.29
    )
    
    # Check if subdirectory for plots exists, if not, create it.
    if not os.path.exists('plots'):
        os.mkdir('plots')
    
    plt.savefig("plots/particle_{0}.pdf".format(n_particle))
    plt.clf()
    plt.close()
    #plt.show()

def log(var):
    """
    Prints log messages if DEBUG is enabled.

    Parameters:
        var: Int, string, whatever...

    Returns:
        None.
    """
    if DEBUG:
            print(var)

def print_options():
    """
    Prints enabled options by the user in the terminal.

    Parameters:
        None
    
    Returns:
        None
    """
    if PLOT_FIGS:
        print("Plotting all particles to '/plots'.\n" + 
              "Disable in {0} for faster run times.".format(__name__) +
              "Or plot specific particle with the 'PLOT_PARTICLE_N' varibale.")
    
    if DEBUG:
        print("Debug mode enabled. More information is printed to terminal.")

if __name__ == "__main__": 
    main()
