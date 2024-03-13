# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies using the mammalian crania 

# Code for voxelisation of mesh to allow for segementation 

# Author: James M. Mulqueeney

# Date Last Modified: 13/03/2024

# Load in libraries 
import os
import trimesh
import numpy as np
import pyvista as pv
import SimpleITK as sitk

# Define a function to calculate adaptive pitch based on mesh properties
def calculate_adaptive_pitch(mesh):
    # Example: Calculate pitch based on the bounding box dimensions
    bounding_box = mesh.bounds
    max_dimension = max(bounding_box[1] - bounding_box[0])
    adaptive_pitch = max_dimension / 500.0  # Adjust the divisor as needed
    return adaptive_pitch

# Step 1: Specify the directories for input meshes and output files
input_directory = 'path/to/your/meshes/directory'
output_directory = 'path/to/your/output/directory'

# Step 2: Get a list of .ply mesh files in the input directory
mesh_files = [file for file in os.listdir(input_directory) if file.endswith('.ply')]

# Step 3: Process each mesh file
for mesh_file in mesh_files:
    # Step 3a: Load the .ply mesh data
    mesh_path = os.path.join(input_directory, mesh_file)
    mesh = trimesh.load_mesh(mesh_path)   
    # Step 3b: Calculate the adaptive pitch for this mesh
    adaptive_pitch = calculate_adaptive_pitch(mesh)    
    # Step 3c: Transform the mesh to the correct position and orientation
    # Step 3d: Convert the mesh to a closed structure using adaptive pitch
    volume = mesh.voxelized(pitch=adaptive_pitch)    
    # Step 3e: Convert the voxel grid to a binary voxel map
    voxel_map = volume.matrix.astype(np.uint8)
    voxel_map = np.flip(voxel_map, axis=0)   
    # Calculate the spacing based on the voxel grid dimensions and bounds
    dimensions = np.array(voxel_map.shape)
    min_bounds, max_bounds = volume.bounds
    spacing = (max_bounds - min_bounds) / (dimensions - 1)    
    # Create the SimpleITK image with correct spacing
    binary_image = sitk.GetImageFromArray(voxel_map)
    binary_image.SetSpacing(spacing.tolist())
    binary_image.SetOrigin(min_bounds)  # Set the origin to match the mesh    
    # Step 3f: Apply binary morphological operations to fill the inner region
    filled_image = sitk.BinaryFillhole(binary_image)    
    # Step 3g: Export the binary voxel map with the filled inner region as a multipage .tiff stack
    filled_voxel_map_array = sitk.GetArrayFromImage(filled_image)
    filled_image = sitk.GetImageFromArray(filled_voxel_map_array)    
    # Set the spacing and origin for the final output
    filled_image.SetSpacing(spacing.tolist())
    filled_image.SetOrigin(min_bounds)    
    # Save the output with the same filename but different extension
    output_path = os.path.join(output_directory, os.path.splitext(mesh_file)[0] + '.tif')
    sitk.WriteImage(filled_image, output_path, useCompression=True)
