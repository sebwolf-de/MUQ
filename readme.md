# MUQ: MIT Uncertainty Quantification Library

Welcome to MUQ (pronounced “muck”), a modular software framework for defining and solving forward and inverse uncertainty quantification problems.

## Purpose
Uncertainty quantification (UQ) is important in many different applications.
MUQ aims to make advanced probabilistic UQ tools easy to use in either c++ or python,
and enable cutting-edge method development through its modular structure.

MUQ has a variety of capabilities, including:

*  Various Markov chain Monte Carlo methods
*  Graphical modeling with a mix of statistical and physical components.
*  Gaussian processes
*  Karhunen Loève expansions.
*  Transport maps
*  Nonlinear Optimization
*  Generalized Polynomial Chaos Expansions

## Installation:
MUQ is available on Linux and OSX as a conda package, docker image, or from source. For many users, getting started can be as easy as running

```
conda install -c conda-forge muq
```

For more installation options, check out the [installation guide](https://mituq.bitbucket.io/source/_site/latest/muqinstall.html).

## Getting Started

You can find an extensive amount of easy to run [examples](https://mituq.bitbucket.io/source/_site/examples.html)
in both c++ and Python as part of our repository.

For more detailed information about MUQ, refer to the [API documentation](https://mituq.bitbucket.io/source/_site/latest/index.html).

If you are mostly interested in modelling, it is a good idea to focus on MUQ's model graphs.
They allow you to construct complex models while keeping a clean code structure.

For method developers, examples showing various existing methods and the API documentation will be a good start to see what components you can build on.

Also, join our Slack channel via our [website](http://muq.mit.edu/) to get in touch with other developers. We are always happy to help!

## Citing
Parno, M., Davis, A., Seelinger, L., and Marzouk, Y. (2014). MIT Uncertainty Quantification (MUQ) library.

```
@misc{MUQ,
  title={MIT uncertainty quantification (MUQ) library},
  author={Parno, Matthew and Davis, Andrew and Seelinger, Linus and Marzouk, Youssef},
  year={2014}
}
```

## Contributing

#### Want to help develop MUQ?
Yes, please! We frequently discuss future developments on Slack ([join via our website](http://muq.mit.edu/)), so feel free to drop by!
Then fork the [muq2 repository](https://bitbucket.org/mituq/muq2/src/master/) and submit a pull request when ready.
Also check out our [style guide](https://mituq.bitbucket.io/source/_site/latest/muqstyle.html).

#### Find a bug?
[Submit the issue on bitbucket](https://bitbucket.org/mituq/muq2/issues/new).  Make sure to label the issue as a bug.

#### Want a new feature?
[Submit a request on bitbucket](https://bitbucket.org/mituq/muq2/issues/new).  Label the issue as an enhancement or proposal.
