# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies using the mammalian crania 

# Code for generating kPCA and eigenvalues .csv files

# Author: James M. Mulqueeney

# Date Last Modified: 13/03/2024

%matplotlib inline

### General Imports
import matplotlib.pyplot as plt
from glob import glob as glob
import numpy as np
import pandas as pd
import os
import matplotlib.cm as cm

### VTK imports
import vtk
from vtk import vtkPolyDataReader
from vtk import vtkPolyDataWriter
from vtk import vtkUnstructuredGridReader
from vtk import vtkUnstructuredGridWriter
from vtk import vtkPolyData
from vtk import vtkPoints
from vtk import vtkCellArray
from vtk import vtkProcrustesAlignmentFilter as Procrustes
from vtk import vtkMultiBlockDataGroupFilter as GroupFilter
from vtk import vtkLandmarkTransform
from vtk import vtkTransformPolyDataFilter as TransformFilter
from vtk.util import vtkConstants
from vtk import vtkIdList
from vtk import vtkIdTypeArray
from vtk import vtkTriangle
from vtk import vtkFloatArray
from vtk import vtkTetra
from vtk import vtkMath

# Load Data
working_directory = os.path.join(os.getcwd())

controlpoints = np.loadtxt(os.path.join(working_directory, 'DeterministicAtlas__EstimatedParameters__ControlPoints.txt'))
f = open(os.path.join(working_directory, 'DeterministicAtlas__EstimatedParameters__Momenta.txt'))
first_line = f.readline().split(' ')
number_of_subjects = int(first_line[0])
number_of_controlpoints = int(first_line[1])
dimension = int(first_line[2])
f.close()
momenta = np.loadtxt(os.path.join(working_directory, os.path.join(working_directory, 'DeterministicAtlas__EstimatedParameters__Momenta.txt')), 
                     skiprows=2)
momenta_linearised = momenta.reshape([number_of_subjects, dimension*number_of_controlpoints])

print('Control Points: {}'.format(number_of_controlpoints))
print('Subjects: {}'.format(number_of_subjects))
print('Dimension: {}'.format(dimension))

# Define populations
df = pd.read_csv (os.path.join(working_directory, 'data.csv'))
df

# Kernel PCA
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA, KernelPCA
from sklearn.svm import SVC
from sklearn.model_selection import permutation_test_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
fig = plt.figure(figsize=(7,5))
idx = [0,1]

n_permutations=1000

kpca = KernelPCA(kernel="rbf", fit_inverse_transform=True, n_components=321, gamma=.0000025)
X_kpca = kpca.fit_transform(momenta_linearised)

for i in range(X_kpca.shape[1]):
    df['PC{}'.format(i+1)] = X_kpca[:, i]

# Eigenvalues
eigenvalues = kpca.lambdas_
eigenvectors = kpca.alphas_ / np.linalg.norm(kpca.alphas_, axis=0)

eig = pd.DataFrame()
eig['PCA dimension'] = ['PC{}'.format(idx+1) for idx in range(len(eigenvalues))]
eig['cum. variability (in %)'] = 100 * np.cumsum(eigenvalues) / np.sum(eigenvalues)
pd.set_option('precision', 2)
eig['lambda'] = eigenvalues
eig['alpha'] = list(np.transpose(eigenvectors))

eig
eig.to_csv(os.path.join(working_directory, 'eigenvalues.csv'))

# Write PCA points on disk

for idx in range(X_kpca.shape[1]):
    df['PC{}'.format(idx+1)] = X_kpca[:, idx]
    
df.to_csv(os.path.join(working_directory, 'kpca.csv'))
df
