#ifndef MUQ_SAMPLING_H
#define MUQ_SAMPLING_H

/** @defgroup sampling Sampling

## Background and Motivation
Uncertainty quantification problems often require computing expectations of the form
\f\[
\bar{h} = \int_\Omega h(x) p(x) dx,
\f\]
where \f$h(x)\f$ is some utility function and \f$p(x)\f$ is the probability density function for the random variable \f$x\f$.  Monte Carlo approximations to \f$\bar{h}\f$ use random realizations \f$x^{(i)}\sim p(x)\f$ to approximate \f$\bar{h}\f$ with an estimator \f$\hat{h}\f$ defined by
\f\[
\hat{h} = \sum_{i=1}^N w_i h(x^{(i)}),
\f\]
where \f$w_i\f$ are appropriately defined weights.  Typically, \f$w_i = N^{-1}\f$.   MUQ's sampling module provides tools for constructing Monte Carlo estimates like \f$\hat{h}\f$.  In particular, MUQ provides a suite of [Markov chain Monte Carlo](\ref mcmc) (MCMC) algorithms for generating samples \f$x^{(i)}\f$ that can be used in Monte Carlo.


- The basics of **defining a sampling problem** in MUQ can be found [here](\ref mcmcprob).
- An introduction to **constructing MCMC algorithms** is provided [here](\ref mcmc).
- An introduction to **dimension independent MCMC algorithms** can also be found [here](\ref disamp).
- **Multi-fidelity and Multi-index MCMC algorithms** are described [here](\ref MIMCMC).

*/

/**
@defgroup mcmc Markov chain Monte Carlo
@ingroup sampling

## Markov chain Monte Carlo (MCMC)
*/

/** @defgroup mcmcprob Getting Started 1: Creating a Sampling Problem
@ingroup mcmc

## Sampling Problems

## Posterior Sampling with MCMC
*/

/** @defgroup mcmcalg Getting Started 2: Defining an MCMC Algorithm
@ingroup mcmc

## Overview

## Chains

## Kernels

## Proposals
*/


/**
@defgroup disamp Dimension-Independent MCMC
@ingroup mcmc

## Dimension-Independent MCMC
Coming soon...
*/

/**
@defgroup MIMCMC Multi-Index MCMC
@ingroup mcmc
*/

#endif // #ifndef MUQ_SAMPLING_H
