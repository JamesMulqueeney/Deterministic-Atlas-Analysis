# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies using the mammalian crania 

# Code for decimating and smoothing meshes 

# Author: James M. Mulqueeney

# Date Last Modified: 13/03/2024

# Load in Libraries 
import os
import numpy as np
from skimage import io
from skimage import measure
import trimesh

# Define input and output directories 
input_dir = "E:\CTData\James Mulqueeney\Mammalian Data\Placental Mammalian Data\Original Files\Aligned Mesh Files"
output_dir = r"E:\CTData\James Mulqueeney\Mammalian Data\Placental Mammalian Data\Decimated Meshes\100,000 Faces"

# Create output directory if it does not exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Loop through input directory
for filename in os.listdir(input_dir):
    if filename.endswith('.ply'):
        # Load the input mesh
        mesh = trimesh.load(os.path.join(input_dir, filename))
        # Smooth the mesh
        smoothed_mesh = trimesh.smoothing.filter_laplacian(mesh, lamb=0.5, iterations=10) 
        # Decimate the mesh
        decimated_mesh = mesh.simplify_quadratic_decimation(50000)
        # Define output file path
        output_filepath = os.path.join(output_dir, os.path.splitext(filename)[0] + '.ply')
        # Save the processed mesh
        decimated_mesh.export(output_filepath)

