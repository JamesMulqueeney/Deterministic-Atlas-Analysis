# Deterministic-Atlas-Analysis

Code linked to the paper: Assessing the application of landmark-free morphometrics to macroevolutionary analyses. DOI: https://doi.org/10.1186/s12862-025-02377-9

If you use this code, please cite: Mulqueeney, J.M., Ezard, T.H. and Goswami, A., 2025. Assessing the application of landmark-free morphometrics to macroevolutionary analyses. BMC Ecology and Evolution, 25(1), p.38.

# Data 
All data stored here is used in the results section of the paper. These data can be analysed/ visualised using the R scripts within the 'R Scripts' Folder. kPCA and eigenvalues from all analyses are stored as followed:

1. `Data_A1_A_binturong_Atlas_Results.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 20.0mm (270 control points) using an Arctictis binturong initial template.
2. `Data_A2_C_calvus_Atlas_Results.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 20.0mm (420 control points) using a Cacajao calvus initial template.
3. `Data_A3_S_morckhoviensis_Atals_Results.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 20.0mm (32 control points) using a Schizodelphis morckhoviensis initial template.
4. `Data_A4-Aligned_Only_k40_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the aligned-only data using a kernel width of 40.0mm (45 control points).
5. `Data_A5-Aligned_Only_k20_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the aligned-only data using a kernel width of 20.0mm (270 control points).
6. `Data_A6-Aligned_Only_k10_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the aligned-only data using a kernel width of 10.0mm (1782 control points).
7. `Data_A7-Poisson_k40_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 40.0mm (45 control points).
8. `Data_A8-Poisson_k20_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 20.0mm (270 control points).
9. `Data_A9-Poisson_k10_kpca.csv`: Kernel principal coordinates for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 10.0mm (1782 control points).
10. `Data_A10-Order_Euclidean_Distance_Correlations.csv`: Correlation measures of Euclidean distance values within each major order (>10 specimens) across both the Aligned-only and Poisson mesh analyses.
11. `Data_A11-PCA_Shape_Data_322.csv`: Principal coordinates for shape data for the manual landmarking approach.
12. `Data_A12-Aligned_Only_k40_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the aligned-only data using a kernel width of 40.0mm (45 control points).
13. `Data_A13-Aligned_Only_k20_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the aligned-only data using a kernel width of 20.0mm (270 control points).
14. `Data_A14-Aligned_Only_k10_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the aligned-only data using a kernel width of 10.0mm (1782 control points).
15. `Data_A15-Poisson_k40_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 40.0mm (45 control points).
16. `Data_A16-Poisson_k20_eigenvalues.csv`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 20.0mm (270 control points).
17. `Data_A17-Poisson_k10_eigenvalues`: Eigenvalues generated for the kernel principal coordinate analysis for deterministic atlas analysis measured using the Poisson mesh data using a kernel width of 10.0mm (1782 control points).
18. `Data_A18-Diet_Disparity_Rates_v1.csv`: Estimated values of morphological disparity and evolutionary rates estimated for each different class of diet.
19. `Data_A19-Locomotion_Disparity_Rates_v1.csv`: Estimated values of morphological disparity and evolutionary rates estimated for each different class of locomotion.
20. `Data_A20-Specimen_Details.csv`: Specimen details and species trait data.
21. `Data_A21-Shape_Data_322.csv`: Generalised Procrustes analysis (GPA) shape data for the manual landmarking approach (not placed into a principal coordinate analysis).
22. `Data_A22-Mirrored_322.csv`: Original landmark data for each species (mirrored for both sides of the cranium). These are used in the alignment of the meshes.
23. `Data_A23-New_Order_Colors.csv`: Colours used for each order in the plots.
24. `Data_A24-Combined_Data_322.csv`: Combined data of both the specimen details and shape measurements.
25. `trees90_95subset.tre`: Tree file used in macroevolutionary analyses. 

# R Scripts 
All the R scripts here are used for data processing, visualisation and statistical analysis of the data. These are as follows: 

1.  `Aligned_Only_Interactive_Principal_Component_Plots.R`: Used to make interactive PC plots for Aligned-only data. 
2.  `Aligned_Only_Principal_Component_Comparison.R`: Used to make ggplot PC plots for Aligned-only data.
3.  `All_Correlation_Measure.R`: Used to compute all different correlations. Note: Need to rune Euclidean distances first.
4.  `Alternative_Atlas__Principal_Component_Plots.R`: Used to make ggplot PC plots for different atlas results.
5.  `Alternative_Atlas_Interactive_Principal_Component_Plots.R`: Used to make interactive PC plots for different atlas results. 
6.  `Batch_Mesh_Alignment.R`: Used to batch align the meshes.
7.  `Batch_Mesh_Decimation.R`: Used to batch decimate the meshes.
8.  `Eigenvalue_Scores_Comparison.R`: Used to compare the eigenvalues across all of the analyses.
9.  `Euclidean_Correlation_Within_Orders.R`: Used to compute statistics for within order correlations.
10.  `Euclidean_Distance_Measures.R`: Used to compute Euclidean distances between specimens and Atlas and used to compare correlations between methods. 
11.  `Landmark_Free_Heatmaps_k10.R`: Used to create heatmaps for landmark-free data (Kernel 10).
12.  `Landmark_Free_Heatmaps_k20.R`: Used to create heatmaps for landmark-free data (Kernel 20).
13.  `Landmark_Free_Heatmaps_k40.R`: Used to create heatmaps for landmark-free data (Kernel 40).
14.  `Landmark_Scheme_Visualisation.R`: Used visualise the landmark scheme.
15.  `Macroevo_Analyses.R`: Used to compute phylogenetic signal, disparity and evolutionary rates.
16.  `Macroevo_Analyses_Momentas.R`: Used to compute phylogenetic signal, disparity and evolutionary rates but using Momentas.
17.  `Mantel_Test.R`: Used to perform Mantel test comparisons between methods.
19.  `Manual_Landmarking_Heatmaps.R`: Used to create heatmaps for landmark data. 
20.  `Manual_Landmarking_Interactive_Principal_Component_Plots.R`: Used to make interactive PC plots for manual landmarking data. 
21.  `Manual_Landmarking_Principal_Component_Plots.R`: Used to make ggplot PC plots for manual landmarking data. 
23.  `Partial_Least_Squares_All_Comparison.R`: Used to perform PLS on all of the analyses. 
26. `Poisson_Mesh_Principal_Component_Comparison.R`: Used to make ggplot PC plots for Poisson Mesh data. 
27. `Poisson_Meshes_Interactive_Principal_Component_Plots.R`: Used to make interactive PC plots for Poisson mesh data.
28. `Protest_Distance_Measures.R`: Measure PROTEST correlation between matrices.
29. `Rates_Disparity_Plots.R`: Used to visualise the rates and disparity results.
30. `Scaled_Landmark_Centroid_Measure.R`: Used to measure centroid sizes. 
31. `Single_Mesh_Alignment.R`: Used to align one mesh to another (single). 

# Python Scripts 
Python scripts are used in the analysis, mainly in the processing and generating of data. These are as follows: 

1. `Batch_Mesh_Centroid_Measurement_v1.py`: Use to measure centroid size based on mesh vertices. 
2. `Batch_Mesh_to_ASCII.py`: Used to batch convert binary .ply meshes into ASCII format. 
3. `Display_Control_Points_Final.py`: Used to display control points on the atlas.
4. `Display_Unmapped_Control_Points_Final.py`: Used to display unmapped control points on the atlas. 
5. `Landmark-Free_Analysis_Mammals.py`: Used to perform shape statistics (kPCA). 
6. `Mammal_Dataset_XML_Generation.py`: Used to generate the data_set.xml file used in the Deterministic Atlas Analysis. 
7. `Mesh_Decimation_Smoothing.py`: Used to decimate & smooth .ply meshes. 
8. `New_Folder_Batch_Ply_to_VTK_Convert.py`: Used to batch convert .ply meshes to .vtk format for use in Deterministic Atlas Analysis.
9. `Variable_Batch_Mesh_to_Label_File_Convertor_final.py`: Used to voxelise mesh.
10. `Write_Data_CSV`: Used to write the data.csv file. 
