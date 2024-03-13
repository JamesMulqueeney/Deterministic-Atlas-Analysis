# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies using the mammalian crania 

# Code for converting binary mesh files into ASCII format 

# Author: James M. Mulqueeney

# Date Last Modified: 13/03/2024

# Load in libraries 
import os
from plyfile import PlyData

# Converting ASCII Code 
def convert_to_ascii(input_file, output_file):
    # Read the non-ASCII PLY file
    plydata = PlyData.read(input_file)    
    # Write the data in ASCII format
    with open(output_file, 'w') as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex {}\n".format(len(plydata.elements[0])))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("element face {}\n".format(len(plydata.elements[1])))
        f.write("property list uchar int vertex_indices\n")
        f.write("end_header\n")       
        for vertex in plydata.elements[0]:
            f.write("{:.6f} {:.6f} {:.6f}\n".format(vertex['x'], vertex['y'], vertex['z']))       
        for face in plydata.elements[1]:
            f.write("{} {} {} {}\n".format(len(face['vertex_indices']), *face['vertex_indices']))

# Use ASCII code to make into a batch process 
def batch_convert_to_ascii(input_directory, output_directory):
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)
    # List all .ply files in the input directory and sort them alphabetically
    ply_files = sorted([f for f in os.listdir(input_directory) if f.endswith(".ply")])
    for filename in ply_files:
        input_path = os.path.join(input_directory, filename)
        output_path = os.path.join(output_directory, filename)
        convert_to_ascii(input_path, output_path)

# Example usage
input_directory = 'input_directory'  # Replace with the path to your input directory
output_directory = 'output_directory'  # Replace with the path to your output directory

batch_convert_to_ascii(input_directory, output_directory)
