# Paper - Comparing landmark-free and manual landmarking methods for macroevolutionary studies using the mammalian crania 

# Code for generating the data_set.xml file used in Deformetrica 

# Author: James M. Mulqueeney

# Date Last Modified: 13/03/2024

# Load in Libraries 
import os
from xml.etree.ElementTree import Element, SubElement, tostring

# Define the directory containing the files
directory = r'path/to/directory'

# Create a list to store the file names
file_list = []

# Iterate through the directory and add VTK files to the list
for filename in os.listdir(directory):
    if filename.endswith('.vtk'):
        file_list.append(filename)

# Sort the file names in alphabetical order
sorted_files = sorted(file_list)

# Create a string for the XML data
xml_data = '<?xml version="1.0"?>\n<data-set>\n'

# Iterate through the sorted file names and add subject elements
for filename in sorted_files:
    xml_data += '    <subject id="' + filename + '">\n'
    xml_data += '        <visit id="cranium">\n'
    xml_data += '            <filename object_id="cranium">' + filename + '</filename>\n'
    xml_data += '        </visit>\n'
    xml_data += '    </subject>\n'

# Close the root element
xml_data += '</data-set>\n'

# Write the XML data to a file
with open('data_set.xml', 'w') as f:
    f.write(xml_data)
