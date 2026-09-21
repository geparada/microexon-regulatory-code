#!/bin/bash

# Set the maximum file size limit (in KB)
max_size=50000  # 50MB

# Find all .ipynb files in the current directory and its subdirectories
for file in $(find . -name "*.ipynb")
do
    # Get the file size in KB
    file_size=$(du -k "$file" | cut -f1)

    # If the file size exceeds the limit, track it with git lfs
    if [ $file_size -gt $max_size ]; then
        git lfs track "$file"
    fi
done

# Add the .gitattributes file to the staging area
git add .gitattributes

# Commit the changes
git commit -m "Track large .ipynb files with git lfs"

