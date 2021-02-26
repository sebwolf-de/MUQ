\page muqinstall Installation

# Getting MUQ
MUQ currently supports Linux and OSX systems.  We do not currently support Windows, but let us know if that is important to you and we may consider it in the future.

On Linux or OSX, there are three primary ways of installing MUQ:
 1. Using the \c muq conda package from the conda-forge channel.  (See instructions below.)
 2. Using the \c muq docker image. (See \subpage docker )
 3. Installing MUQ from source. (See \subpage source_install )


# Conda
A conda package containing the c++ and python components of MUQ is available via conda-forge.   To install and use MUQ using conda, use
```
conda install -c conda-forge muq
```

You may run into conflicts when installing MUQ and Fenics with conda if they are installed one at a time.  To avoid this, we recommend installing MUQ and Fenics at the same time:
```
conda install -c conda-forge muq fenics
```
