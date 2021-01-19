#!/bin/sh

# Make sure globstar is enabled
shopt -s globstar
for i in **/*.ipynb; do # Whitespace-safe and recursive
    echo "======================"
    echo "Running notebook $i..."

    papermill --cwd "$(dirname "$i")" "$i" -

    echo "======================"
done


