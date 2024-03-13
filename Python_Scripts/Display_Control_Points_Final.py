# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies using the mammalian crania 

# Code for displaying control points on the intial atlas template 

# Author: James M. Mulqueeney

# Date Last Modified: 13/03/2024

import vtk
import numpy as np
import os

# Change to selected directory
os.chdir("name/of/directory")

# Load the mesh from the .vtk file
vtk_file_path = "Inputs/initial_template.vtk"
reader = vtk.vtkPolyDataReader()
reader.SetFileName(vtk_file_path)
reader.Update()
mesh = reader.GetOutput()

# Parse the control points from the text file
control_points_file_path = "DeterministicAtlas__EstimatedParameters__ControlPoints.txt"
control_points = np.loadtxt(control_points_file_path)

# Create a KDTree for efficient nearest neighbor search
kd_tree = vtk.vtkKdTreePointLocator()
kd_tree.SetDataSet(mesh)
kd_tree.BuildLocator()

# Map control points onto the mesh
mapped_points = []
for point in control_points:
    # Transform control point to mesh coordinate system
    pointId = kd_tree.FindClosestPoint(point[:3])  # Use only the first 3 components (x, y, z)
    mesh_point = mesh.GetPoint(pointId)
    mapped_points.append(mesh_point)

# Create a vtkPoints object for control points
control_point_points = vtk.vtkPoints()
for point in mapped_points:
    control_point_points.InsertNextPoint(point)

control_point_polydata = vtk.vtkPolyData()
control_point_polydata.SetPoints(control_point_points)

# Save mapped points to a text file
output_file_path = "mapped_points.txt"
np.savetxt(output_file_path, mapped_points)

# Create spheres for control points
sphere_source = vtk.vtkSphereSource()
sphere_source.SetRadius(1.0)  # Adjust the sphere radius as needed

glyph = vtk.vtkGlyph3D()
glyph.SetInputData(control_point_polydata)
glyph.SetSourceConnection(sphere_source.GetOutputPort())

# Create a mapper for the control points
mapper_points = vtk.vtkPolyDataMapper()
mapper_points.SetInputConnection(glyph.GetOutputPort())

# Create an actor for the control points
actor_points = vtk.vtkActor()
actor_points.SetMapper(mapper_points)
actor_points.GetProperty().SetColor(1.0, 0.0, 0.0)  # Set color to red

# Create a mapper for the mesh
mapper_mesh = vtk.vtkPolyDataMapper()
mapper_mesh.SetInputData(mesh)

# Create an actor for the mesh
actor_mesh = vtk.vtkActor()
actor_mesh.SetMapper(mapper_mesh)
actor_mesh.GetProperty().SetColor(0.87, 0.79, 0.69)  # Set color to #DECAB0 (RGB values)

# Set specular reflection properties for a shiny appearance
actor_mesh.GetProperty().SetSpecular(0.5)  # Adjust the specular intensity as needed
actor_mesh.GetProperty().SetSpecularPower(30)  # Adjust the specular power as needed

renderer = vtk.vtkRenderer()
renderer.AddActor(actor_mesh)
renderer.AddActor(actor_points)
renderer.SetBackground(1.0, 1.0, 1.0)  # Set background to white

# Create a render window and set its size
render_window = vtk.vtkRenderWindow()
render_window.SetSize(800, 600)  # Set the desired size (width, height)

# Set the renderer to the render window
render_window.AddRenderer(renderer)

# Create a render window interactor
render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)

# Render the scene and start the interaction
render_window.Render()
render_window_interactor.Start()
