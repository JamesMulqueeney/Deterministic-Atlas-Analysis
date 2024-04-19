# Deterministic-Atlas-Analysis

Code linked to the paper: Comparing landmark-free and manual landmarking methods for macroevolutionary studies using the mammalian crania. DOI: 

# Data 
All data stored here is used in the results section of the paper. These data can be analysed/ visualised using the R scripts within the 'R Scripts' Folder. kPCA and eigenvalues from all analyses are stored as followed:

1. `Data_S1-Specimen_Details.csv`: Specimen details and species trait data.
2. `Data_S2-Mirrored_322.csv`: Original landmark data for each species (mirrored for both sides of the cranium). These are used in the alignment of the meshes. 
3. `Data_S3-Shape_Data_322.csv`: Base shape data for the manual landmarking approach (not placed into a principal coordinate analysis). 
4. `Data_S4-PCA_Shape_Data_322.csv`:  Principal coordinates for shape data for the manual landmarking approach. 
5. `Data_S5-Aligned_Only_k40_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the aligned-only data using a kernel width of 40.0mm (45 control points).
6. `Data_S6-Aligned_Only_k20_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the aligned-only data using a kernel width of 20.0mm (270 control points).
7. `Data_S7-Aligned_Only_k10_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the aligned-only data using a kernel width of 10.0mm (1782 control points).
8. `Data_S8-Possion_k40_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 40.0mm (45 control points).
9. `Data_S9-Poisson_k20_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 20.0mm (270 control points).
10. `Data_S10-Poisson_k10_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 10.0mm (1782 control points).
11. `Data_S11-Aligned_Only_k40_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the aligned-only data using a kernel width of 40.0mm (45 control points).
12. `Data_S12-Aligned_Only_k20_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the aligned-only data using a kernel width of 20.0mm (270 control points).
13. `Data_S13-Aligned_Only_k10_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the aligned-only data using a kernel width of 10.0mm (1782 control points).
14. `Data_S14-Poisson_k40_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 40.0mm (45 control points).
15. `Data_S15-Poisson_k20_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 20.0mm (270 control points).
16. `Data_S16-Poisson_k10_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 10.0mm (1782 control points).

# R Scripts 
All the R scripts here are used for data processing, visualisation and statistical analysis of the data. These are as follows: 

1.  `Aligned_Only_Interactive_Principal_Component_Plots.R`: Used to make interactive PC plots for Aligned-only data. 
2.  `Aligned_Only_Principal_Component_Comparison.R`: Used to make ggplot PC plots for Aligned-only data.
3.  `Batch_Mesh_Alignment.R`: Used to batch align the meshes. 
4.  `Landmark_Free_Heatmaps.R`: Used to create heatmaps for landmark-free data. 
5.  `Manual_Landmarking_Heatmaps.R`: Used to create heatmaps for landmark data. 
6.  `Manual_Landmarking_Interactive_Principal_Component_Plots.R`: Used to make interactive PC plots for manual landmarking data. 
7.  `Manual_Landmarking_Principal_Component_Plots.R`: Used to make ggplot PC plots for manual landmarking data. 
8.  `Partial_Least_Squares_Aligned_Only_Comparison.R`: Used to perform PLS on Aligned-only data. 
9.  `Partial_Least_Squares_Poisson_Meshes_Comparison.R`: Used to perform PLS on the Poisson Mesh data. 
10. `Poisson_Mesh_Principal_Component_Comparison.R`: Used to make ggplot PC plots for Poisson Mesh data. 
11. `Poisson_Meshes_Interactive_Principal_Component_Plots.R`: Used to make interactive PC plots for Poisson mesh data. 
12. `Single_Mesh_Alignment.R`: Used to align one mesh to another (single). 

# Python Scripts 
Python scripts are used in the analysis, mainly in the processing and generating of data. These are as follows: 

1. `Batch_Mesh_to_ASCII.py`: Used to batch convert binary .ply meshes into ASCII format. 
2. `Display_Control_Points_Final.py`: Used to display control points on atlas 
3. `Landmark-Free_Analysis_Mammals.py`: Used to perform shape statistics (kPCA) 
4. `Mammal_Dataset_XML_Generation.py`: Used to generate the data_set.xml file used in the Deterministic Atlas Analysis. 
5. `Mesh_Decimation_Smoothing.py`: Used to decimate & smooth .ply meshes. 
6. `New_Folder_Batch_Ply_to_VTK_Convert.py`: Used to batch convert .ply meshes to .vtk format for use in Deterministic Atlas Analysis.
7. `Variable_Batch_Mesh_to_Label_File_Convertor_final.py` : Used to voxelise mesh 
