#!/bin/bash

exit_code=0

# Make sure globstar is enabled
shopt -s globstar
for notebook in **/*.ipynb; do # Whitespace-safe and recursive
    echo "======================"
    echo "Running $notebook:"
    echo ""

    papermill --cwd "$(dirname "$notebook")" "$notebook" -

    papermill_exit=$?

    if [ $papermill_exit -eq 0 ]
    then
        echo "Notebook $notebook successful."
    else
        echo "Notebook $notebook FAILED!"
    fi

    exit_code=$(($exit_code + $papermill_exit)) # Simply add exit codes; will be zero if all tests successful
done

exit $exit_code
