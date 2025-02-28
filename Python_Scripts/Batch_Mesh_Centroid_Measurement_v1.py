import vtk
import numpy as np
import os
import csv

def centroid(vtk_points):
    cp = [0] * 3
    np_points = vtk_points.GetNumberOfPoints()
    for i in range(np_points):
        p = vtk_points.GetPoint(i)
        cp[0] += p[0]; cp[1] += p[1]; cp[2] += p[2];
    cp[0] /= np_points; cp[1] /= np_points; cp[2] /= np_points;
    return cp

def centroid_size(vtk_points):
    cp = centroid(vtk_points)
    n = vtk_points.GetNumberOfPoints()
    S = 0
    for i in range(n):
        p = vtk_points.GetPoint(i)
        S += vtk.vtkMath.Distance2BetweenPoints(p, cp)
    return np.sqrt(S)

def get_mesh_files(directory, extension=".vtk"):
    return [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(extension)]

def batch_process(directory):
    mesh_files = get_mesh_files(directory)
    centroid_sizes = []

    for f in mesh_files:
        r = vtk.vtkPolyDataReader()
        r.SetFileName(f)
        r.Update()
        centroid_size_value = centroid_size(r.GetOutput().GetPoints())
        centroid_sizes.append((os.path.basename(f), centroid_size_value))

    return centroid_sizes

def save_to_csv(centroid_sizes, output_file):
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Mesh File", "Centroid Size"])
        for mesh_file, size in centroid_sizes:
            writer.writerow([mesh_file, size])

# Example usage
directory = 'E:\James Mulqueeney\Paper 2- Mammal Shape\Final Results and Code\Large Results\Final Poisson Analysis\Arctictis_binturong_atlas\Kernel 20.0'
output_csv = r'E:\James Mulqueeney\Paper 2- Mammal Shape\Final Results and Code\Deterministic-Atlas-Analysis-main (2)\Deterministic-Atlas-Analysis-main\Data\New Data\Centroid Data\Mesh_centroid_sizes.csv'

centroid_sizes = batch_process(directory)
save_to_csv(centroid_sizes, output_csv)

print(f"Centroid sizes saved to {output_csv}")
