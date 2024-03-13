# Deterministic-Atlas-Analysis

Code linked to the paper: Comparing landmark-free and manual landmarking methods for macroevolutionary studies using the mammalian crania. DOI: 

# Data 
All data stored here is used in the results section of the paper. These data can be analysed/ visualised using the R scripts wihin the 'R Scripts' Folder. kPCA and eigenvalues from all analyses are stored as followed:

Aligned-Only : Kernel 40.0, Kernel 20.0, Kernel 10.0 

Poisson_Meshes: Kernel 40.0, Kernel 20.0, Kernel 10.0 

# R Scripts 
All the R scripts here are used for data processing, visualisation and statisitcal analysis of the data. These are as follows: 

1.  `Aligned_Only_Interactive_Principal_Component_Plots.R`
2.  `Aligned_Only_Principal_Component_Comparison.R`
3.  `Batch_Mesh_Alignment.R`
4.  `Landmark_Free_Heatmaps.R`
5.  `Manual_Landmarking_Heatmaps.R`
6.  `Manual_Landmarking_Interactive_Principal_Component_Plots.R`
7.  `Manual_Landmarking_Principal_Component_Plots.R`
8.  `Partial_Least_Squares_Aligned_Only_Comparison.R`
9.  `Partial_Least_Squares_Poisson_Meshes_Comparison.R`
10. `Poisson_Mesh_Principal_Component_Comparison.R`
11. `Poisson_Meshes_Interactive_Principal_Component_Plots.R`
12. `Single_Mesh_Alignment.R`

# Python Scripts 
Python scripts used in the analysis, mainly in the processing and generating data. These are as follows: 

1. `Batch_Mesh_to_ASCII.py`: Used to batch convert binary .ply meshes into ASCII format. 
2. `Display_Control_Points_Final.py`: Used to display control points on atlas 
3. `Landmark-Free_Analysis_Mammals.py`: Used to perform shape statistics (kPCA) 
4. `Mammal_Dataset_XML_Generation.py`: Used to generate the data_set.xml file used in the Deterministic Atlas Analysis. 
5. `Mesh_Decimation_Smoothing.py`: Used to decimate & smooth .ply meshes. 
6. `New_Folder_Batch_Ply_to_VTK_Convert.py`: Used to batch convert .ply meshes to .vtk format for use in Deterministic Atlas Analysis.
7. `Variable_Batch_Mesh_to_Label_File_Convertor_final.py` : Used to voxelise mesh 
