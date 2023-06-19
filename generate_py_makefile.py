

import os
import os.path

# Set the directory where the notebooks are located
NOTEBOOK_DIR = 'numcosmo_py'

# Get a list of all subdirectories in the notebook directory
subdirs = []
for root, dirs, files in os.walk(NOTEBOOK_DIR):
    if "__pycache__" not in root:
        normalized_path = os.path.normpath(root)
        if normalized_path.startswith("./"):
            normalized_path = normalized_path[2:]
        subdirs.append(normalized_path)

subdirs.sort()

# Generate the Makefile.am contents
makefile_contents = "## Generated automatically with generate_py_makefile.py\n\n\n"
python_files = []

for subdir in subdirs:
    # Get a list of all notebook files in the subdirectory
    for root, dirs, files in os.walk(subdir):
        if subdir == root:
            for file in files:
                if file.endswith(('.py', '.pyi')):
                    file = os.path.join(root, file)
                    python_files.append(file)
                

max_length = max(len(extra_dist_file) for extra_dist_file in python_files)
python_files = [extra_dist_file.ljust(max_length) for extra_dist_file in python_files]

# Generate the EXTRA_DIST variable
makefile_contents += "nobase_python_PYTHON = \\\n"
for file in python_files[:-1]:
    makefile_contents += f"\t{file} \\\n"
makefile_contents += f"\t{python_files[-1]}\n"

# Write the Makefile.am file
with open('numcosmo_py.mk', 'w') as f:
    f.write(makefile_contents)

