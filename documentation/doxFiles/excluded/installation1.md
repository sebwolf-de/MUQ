\page muqinstall MUQ Installation Guide

# Getting MUQ
MUQ currently supports Linux and OSX systems.  We do not currently support Windows, but let us know if that is important to you and we may consider it in the future.

On Linux or OSX, there are three primary ways of installing MUQ:
 1. Using the \c muq docker image
 2. Using the \c muq conda package from the conda-forge channel
 3. Installing MUQ from source.

# Conda
A conda package containing the c++ and python components of MUQ is available via conda-forge.   To install and use MUQ using conda, use
```
conda install -c conda-forge muq
```

You may run into conflicts when installing MUQ and Fenics with conda if they are installed one at a time.  To avoid this, we recommend installing MUQ and Fenics at the same time:
```
conda install -c conda-forge muq fenics
```

# Docker
First, make sure you have docker installed.  Follow [these](https://docs.docker.com/get-started/) instructions if you don't.  We provide several MUQ docker images:

| Image Name | Description  |
|---|---|
| `mparno/muq`  | A simple image containing MUQ (c++ and python) installed in a debian environment with minimal dependencies.  |
| `mparno/muq-jupyter` | Same as the muq image, but with jupyter lab and several other python packages installed. |
| `mparno/muq-hippylib`  | An image building on the [FEniCS](https://fenicsproject.org/) `quay.io/fenicsproject/stable` image with additional installations of [hIPPylib](https://hippylib.github.io/), [hIPPylib2MUQ](https://github.com/hippylib/hippylib2muq), andd MUQ.  |
| `mparno/muq-build` | An image containing necessary dependencies for MUQ, but not MUQ itself.  This image is useful for continuous integration purposes and serves as the base for the `mparno/muq` and `mparno/muq-jupyter` images. |

The following command launches a container using the `muq` image and sharing the current directory with the `/home/muq-user` folder in the container:
```
docker pull muq
docker run -it -v ${PWD}:/home/muq-user muq bash
```

To launch a container using the `muq-jupyter` image, we should also share port 8888 so we can launch jupyter notebooks inside the container but view them in a host browser.   An example doing this is
```
docker run -it -p 8888:8888 -v ${PWD}:/home/muq-user muq-jupyter bash
```


# Building from source
 MUQ2 is hosted on<a href=https://bitbucket.org/mituq/muq2>BitBucket</a>.  To get the source code, you can clone our git repository by running
 ```
 git clone https://bitbucket.org/mituq/muq2
 ```
 While it is usually better to clone MUQ so you can easily obtain updates and bugfixes, you can also download a version of MUQ <a href=https://bitbucket.org/mituq/muq2/downloads>here</a>.

## QuickStart
The following commands can be used on Ubuntu to compile MUQ and its examples.
```
