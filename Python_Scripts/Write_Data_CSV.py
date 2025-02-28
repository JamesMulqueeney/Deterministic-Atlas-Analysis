import xml.etree.ElementTree as ET
import csv

# Load the XML file
file_path = 'F:/aligned-binary-meshes/Yichen Meshes/VTK Files/data_set.xml'
tree = ET.parse(file_path)
root = tree.getroot()

# Extract and list <filename> elements in the order they appear
filenames = []
for filename in root.findall('.//filename'):
    filenames.append(filename.text)

# Specify the path for the output CSV file
output_file = 'F:/aligned-binary-meshes/Yichen Meshes/VTK Files/data.csv'

# Write the filenames to the CSV file
with open(output_file, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    # Write header
    csvwriter.writerow(['Filename'])
    # Write each filename
    for name in filenames:
        csvwriter.writerow([name])

print(f'Filenames have been written to {output_file}')
