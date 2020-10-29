\page muqinstall MUQ Installation Guide

# Getting MUQ2
MUQ currently supports Linux and OSX systems.  We do not currently support Windows, but let us know if that is important to you and we may consider it in the future.

On Linux or OSX, there are three primary ways of installing MUQ:
 1. Using the \c muq docker image
 2. Using the \c muq conda package from the conda-forge channel
 3. Installing MUQ from source.

## Docker
First, make sure you have docker installed.  Follow [these](https://docs.docker.com/get-started/) instructions if you don't.  We provide several MUQ docker images:

| Image Name | Description  |
|---|---|
| `mparno/muq`  | A simple image containing MUQ (c++ and python) installed in a debian environment with minimal dependencies.  |
| `mparno/muq-jupyter` | Same as the muq image, but with jupyter lab and several other python packages installed. |
| `mparno/muq-hippylib`  | An image building on the [FEniCS](https://fenicsproject.org/) `quay.io/fenicsproject/stable` image with an additional installation of [hIPPylib](https://hippylib.github.io/).  |
| `mparno/muq-build` | An image containing necessary dependencies for MUQ, but not MUQ itself.  This image is useful for continuous integration purposes. |

The following command launches a container using the `muq` image and sharing the current directory with the `/home/muq-user` folder in the container:
```
docker pull muq
docker run -it -v ${PWD}:/home/muq-user muq bash
```

To launch a container using the `muq-jupyter` image, we should also share port 8888 so we can launch jupyter notebooks inside the container but view them in a host browser.   An example doing this is
```
docker run -it -p 8888:8888 -v ${PWD}:/home/muq-user muq-jupyter bash
```

## Conda
A conda package containing the c++ and python components of MUQ will soon be available on conda-forge.   Once accepted, this package will enable you to install MUQ with
```
conda install -c conda-forge muq
```

## Building from source
 MUQ2 is hosted on<a href=https://bitbucket.org/mituq/muq2>BitBucket</a>.  To get the source code, you can clone our git repository by running
 ```
 git clone https://bitbucket.org/mituq/muq2
 ```
 While it is usually better to clone MUQ so you can easily obtain updates and bugfixes, you can also download a version of MUQ <a href=https://bitbucket.org/mituq/muq2/downloads>here</a>.

### Configuring the Build
MUQ is written in c++ with a python wrapper provided through <a href=https://pybind11.readthedocs.io>pybind11</a>.   CMake is used to configure the MUQ build and generate makefiles.   Like most projects based on CMake, a basic configuration of MUQ can be accomplished by running CMake and only specifying the installation path.  For example, from the `muq2/build` directory, run
```
cmake -DCMAKE_INSTALL_PREFIX=/my/install/path ..
```

During this call, CMake will search for all necessary dependencies, test those dependencies, summarize the configuration in a file called `summary.log`, and create a makefile capable of compiling the MUQ libraries.

Additional configuration options can also be added to the CMake command to control the build process; the most important of these are described in the subsections below.

#### Handling Dependencies
Various components of MUQ use thee following external libraries
- <a href=http://www.boost.org/>Boost</a>
- <a href=http://eigen.tuxfamily.org/>Eigen</a>
- <a href=http://www.hdfgroup.org/HDF5/>HDF5</a>
- <a href=http://mc-stan.org/>Stan Math</a>
- <a href=https://computation.llnl.gov/projects/sundials>Sundials</a>
- <a href=https://github.com/jlblancoc/nanoflann>nanoflann</a>

#### Compile Groups
On a high level, MUQ is organized into separate libraries (e.g., muqModeling, muqApproximation). However, users have more granular control over what components of MUQ are compiled through our "compile group" concept.   A compile group is simply a set of related c++ source files with common dependencies that can be enabled or disabled during the CMake configuration.  If a group is disabled, it will not be compiled and its dependencies will not used (or even searched for).

### Compiling

### Testing
