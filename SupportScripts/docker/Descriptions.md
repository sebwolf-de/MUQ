## Descriptions

| Image Name | Description  |
|---|---|
| `mparno/muq`  | A simple image containing MUQ (c++ and python) installed in a debian environment with minimal dependencies.  |
| `mparno/muq-jupyter` | Same as the muq image, but with jupyter lab and several other python packages installed. |
| `mparno/muq-hippylib`  | An image building on the [FEniCS](https://fenicsproject.org/) `quay.io/fenicsproject/stable` image with additional installations of [hIPPylib](https://hippylib.github.io/), [hIPPylib2MUQ](https://github.com/hippylib/hippylib2muq), andd MUQ.  |
| `mparno/muq-build` | An image containing necessary dependencies for MUQ, but not MUQ itself.  This image is useful for continuous integration purposes and serves as the base for the `mparno/muq` and `mparno/muq-jupyter` images. |
