# MUQ: MIT Uncertainty Quantification Library

Welcome to MUQ (pronounced “muck”).  A collection of tools for defining and solving forward and inverse Bayesian uncertainty quantification problems.

## Purpose
Uncertainty quantification UQ is important in many different applications and MUQ aims to make advanced probabilistic UQ tools easy to use in either c++ or python.   MUQ's primary emphasis is on setting up and solving Bayesian inference problems, but has some tools for forward UQ as well.

MUQ has a variety of capabilities, including
- Markov chain Monte Carlo
- Graphical modeling with a mix of statistical and physical components.
- Gaussian processes
- Karhunen Lo&egrave;ve expansions.
- Transport maps
- Nonlinear Optimization
- Generalized Polynomial Chaos Expansions

## Installation:
MUQ is available on Linux and OSX as a conda package, docker image, or from source.   For many users, getting started can be as easy as running
```
conda install -c conda-forge muq
```
For more details on other installation options, check out the [installation guide](@ref muqinstall).

## Getting Started
MUQ has many features but most of them rely on defining your model in a way that MUQ can understand (i.e., as a child of the muq::Modeling::ModPiece class).  Check out the `CustomModPiece` example to get started creating your own ModPiece.  The other examples described [here](http://muq.mit.edu/examples/) can also serve as a jumping-off point for your specific application.

More detailed [API documentation](http://muq.mit.edu/master-muq2-docs/) is generated with doxygen.

## Citing
Parno, M., Davis, A., Seelinger L., & Marzouk, Y. (2014). MIT Uncertainty Quantification (MUQ) library.

```
@misc{MUQ,
  title={MIT uncertainty quantification (MUQ) library},
  author={Parno, Matthew and Davis, Andrew and Linus, Seelinger and Marzouk, Youssef},
  year={2014}
}
```

## Contributing

#### Want to help develop MUQ?
Yes, please! Fork the [muq2](https://bitbucket.org/mituq/muq2/src/master/) repository and submit a pull request when ready.  Also check out our [style guide](@ref muqstyle) for details on what we expect in the submission. </br>

#### Find a bug?
[Submit the issue on bitbucket](https://bitbucket.org/mituq/muq2/issues/new).  Make sure to label the issue as a bug.</br>

#### Want a new feature?
[Submit a request on bitbucket](https://bitbucket.org/mituq/muq2/issues/new).  Label the issue as and enhancement or proposal.
