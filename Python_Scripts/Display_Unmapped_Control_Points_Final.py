import vtk
import numpy as np
import os

# Change to selected directory
os.chdir("E:\James Mulqueeney\Paper 2- Mammal Shape\Final Results and Code\Large Results\Final Poisson Analysis\Arctictis_binturong_atlas\Kernel 40.0\output")

# Load the mesh from the .vtk file
vtk_file_path = "initial_template.vtk"
reader = vtk.vtkPolyDataReader()
reader.SetFileName(vtk_file_path)
reader.Update()
mesh = reader.GetOutput()

# Parse the control points from the text file
control_points_file_path = "DeterministicAtlas__EstimatedParameters__ControlPoints.txt"
control_points = np.loadtxt(control_points_file_path)

# Create a vtkPoints object for control points (original)
original_points_vtk = vtk.vtkPoints()
for point in control_points:
    original_points_vtk.InsertNextPoint(point[:3])  # Insert original control points (x, y, z)

# Create a vtkPolyData object for original points
original_point_polydata = vtk.vtkPolyData()
original_point_polydata.SetPoints(original_points_vtk)

# Create spheres for control points
sphere_source = vtk.vtkSphereSource()
sphere_source.SetRadius(1.0)  # Adjust the sphere radius as needed

# Create glyphs for original points
glyph_original = vtk.vtkGlyph3D()
glyph_original.SetInputData(original_point_polydata)
glyph_original.SetSourceConnection(sphere_source.GetOutputPort())

# Create a mapper and actor for original points
mapper_original_points = vtk.vtkPolyDataMapper()
mapper_original_points.SetInputConnection(glyph_original.GetOutputPort())

actor_original_points = vtk.vtkActor()
actor_original_points.SetMapper(mapper_original_points)
actor_original_points.GetProperty().SetColor(0.0, 0.0, 1.0)  # Set color to blue

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

# Set up the renderer
renderer = vtk.vtkRenderer()
renderer.AddActor(actor_mesh)
renderer.AddActor(actor_original_points)  # Add original points to the renderer
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
