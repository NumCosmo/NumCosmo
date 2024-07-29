

import os
import os.path

# Set the directory where the notebooks are located
NOTEBOOK_DIR = '.'

# Get a list of all subdirectories in the notebook directory
subdirs = []
for root, dirs, files in os.walk(NOTEBOOK_DIR):
    if root != NOTEBOOK_DIR:
        if ".ipynb_checkpoints" not in root:
            normalized_path = os.path.normpath(root)
            if normalized_path.startswith("./"):
                normalized_path = normalized_path[2:]
            subdirs.append(normalized_path)
    else:
        subdirs.append(".")

subdirs.sort()

# Generate the Makefile.am contents
makefile_contents = "## Process this file with automake to produce Makefile.in\n\n"
makefile_contents += "SUBDIRS = .\n\n"
makefile_contents += "nbs_datadir = $(pkgdatadir)-$(VERSION)/notebooks\n\n"
extra_dist_files = []

for subdir in subdirs:
    # Get a list of all notebook files in the subdirectory
    notebook_files = []
    for root, dirs, files in os.walk(subdir):
        if subdir == root:
            for file in files:
                if file.endswith('.ipynb'):
                    if subdir != ".":
                        file = os.path.join(root, file)
                    notebook_files.append(file)
                
    if len(notebook_files) == 0:
        continue
    max_length = max(len(notebook_file) for notebook_file in notebook_files)
    notebook_files = [notebook_file.ljust(max_length) for notebook_file in notebook_files]
    notebook_files.sort()

    # Generate the _DATA variable for the subdirectory
    if subdir == ".":
        subdir_varname = "nbs"
        makefile_contents += f"{subdir_varname}dir = $(nbs_datadir)\n"
        makefile_contents += f"{subdir_varname}_DATA = \\\n"
    else:
        subdir_varname = subdir.replace('/', '_').replace('.', '').lower()
        makefile_contents += f"{subdir_varname}dir = $(nbs_datadir)/{subdir}\n"
        makefile_contents += f"{subdir_varname}_DATA = \\\n"

    extra_dist_files.append(f"$({subdir_varname}_DATA)")
    
    for file in notebook_files[:-1]:
        makefile_contents += f"\t{file} \\\n"
    makefile_contents += f"\t{notebook_files[-1]}\n"
    makefile_contents += "\n\n"


max_length = max(len(extra_dist_file) for extra_dist_file in extra_dist_files)
extra_dist_files = [extra_dist_file.ljust(max_length) for extra_dist_file in extra_dist_files]

# Generate the EXTRA_DIST variable
makefile_contents += "EXTRA_DIST = \\\n"
for file in extra_dist_files[:-1]:
    makefile_contents += f"\t{file} \\\n"
makefile_contents += f"\t{extra_dist_files[-1]}\n"

# Write the Makefile.am file
with open('Makefile.am', 'w') as f:
    f.write(makefile_contents)

