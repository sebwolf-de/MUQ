#!/bin/bash

cmake \
  -DCMAKE_INSTALL_PREFIX=~/Installations/MUQ_INSTALL/ \
  -DCMAKE_BUILD_TYPE=Release \
  -DMUQ_USE_NLOPT=ON \
  -DMUQ_NLOPT_DIR=~/Installations/nlopt \
  -DMUQ_USE_PYTHON=OFF \
  -DMUQ_USE_DOLFIN=OFF \
  -DMUQ_USE_GTEST=ON \
  -DMUQ_GTEST_DIR=~/Installations/gtest_install \
  -DMUQ_ENABLEGROUP_DEFAULT=ON \
  -DMUQ_ENABLEGROUP_MODELING_CORE_PYTHON=OFF \
  -DMUQ_ENABLEGROUP_MODELING_CORE=ON \
  -DMUQ_ENABLEGROUP_MODELING_DOLFIN=OFF \
  -DMUQ_ENABLEGROUP_APPROXIMATION_GP_Kernels=ON \
../