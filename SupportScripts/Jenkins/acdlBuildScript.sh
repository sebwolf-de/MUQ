#!/bin/bash

shopt -s nocasematch

#######################################
##### DECIDE TO BUILD RELEASE OR DEBUG
#######################################
dir=`pwd`
bin_dir=$(echo "$dir" | tr '[:lower:]' '[:upper:]')

build_type="Debug"
if echo $bin_dir| grep -q "RELEASE"; then
  build_type="Release"
else
  build_type="Debug"
fi

#######################################
##### DECIDE TO BUILD WITH PYTHON OR NOT
#######################################
if echo $bin_dir| grep -q "PYTHON"; then
  with_python="ON"
else
  with_python="OFF"
fi

#######################################
##### DECIDE TO BUILD WITH NLOPT OR NOT
#######################################
if echo $bin_dir| grep -q "NLOPT"; then
  with_nlopt="ON"
else
  with_nlopt="OFF"
fi

####################################
##### SET MACHINE SPECIFIC PATHS
####################################

if [[ `hostname` == "reynolds" ]]; then

  BOOST_SOURCE=/home/mparno/util/boost_1_63_0.tar.gz
    
  if echo $bin_dir| grep -q "CLANG"; then
    GTEST_DIR=/home/mparno/util/gtest/clang
    echo "Using gtest compiled with clang in ${GTEST_DIR}"
  else
    GTEST_DIR=/home/mparno/util/gtest/gnu
    echo "Using gtest compiled with gcc in ${GTEST_DIR}"
  fi

  NLOPT_DIR=/home1/jenkins/util/NLOPT/nlopt-2.4.2/install

elif [[ `hostname` == "macys.mit.edu" ]]; then

  BOOST_SOURCE=/Users/mparno/util/boost_1_63_0.tar.gz
    
  # export the path so we can get cmake
  export PATH=/usr/local/bin:/usr/local/sbin:$PATH

  GTEST_DIR=/Users/jenkins/util/gtest/
  NLOPT_DIR=/usr/local/Cellar/nlopt/2.4.2/

else

  # if not on macys or reynolds, assume I'm on Matt's macbook for testing
  NLOPT_DIR=/Users/mparno/Desktop/NLOPT_CLANG
  GTEST_DIR=/Users/mparno/Documents/Repositories/gtest_clang/

fi

#######################################
##### EXTRACT COMPILER FROM WORKSPACE
#######################################

if echo $bin_dir| grep -q "CLANG35"; then

  if [[ `hostname` == "reynolds" ]]; then
    my_cc_compiler="clang-3.5"
    my_cxx_compiler="clang++-3.5"
  else
    my_cc_compiler="clang"
    my_cxx_compiler="clang++"
  fi

elif echo $bin_dir | grep -q "CLANG"; then

  if [[ `hostname` == "macys.mit.edu" ]]; then
    my_cc_compiler="clang"
    my_cxx_compiler="clang++"
  else
    my_cc_compiler="clang-3.5"
    my_cxx_compiler="clang++-3.5"
  fi

elif echo $bin_dir | grep -q "INTEL"; then
  my_cc_compiler="icc"
  my_cxx_compiler="icpc"
elif echo $bin_dir | grep -q "GNU49"; then
  my_cc_compiler="gcc-4.9"
  my_cxx_compiler="g++-4.9"
elif echo $bin_dir | grep -q "GNU48"; then
  my_cc_compiler="gcc-4.8"
  my_cxx_compiler="g++-4.8"
elif echo $bin_dir | grep -q "GNU47"; then
  my_cc_compiler="gcc-4.7"
  my_cxx_compiler="g++-4.7"
elif echo $bin_dir | grep -q "GNU5"; then
  my_cc_compiler="gcc-5"
  my_cxx_compiler="g++-5"
elif echo $bin_dir | grep -q "GNU"; then
  my_cc_compiler="gcc"
  my_cxx_compiler="g++"
fi

echo "C Compiler = $my_cc_compiler"
echo "CXX Compiler = $my_cxx_compiler"


#######################################
##### CHANGE INTO BUILD DIR
#######################################
BUILD_DIR="$dir/build"
INSTALL_DIR="$dir/install"
echo "BUILD_DIR = $BUILD_DIR"
echo "SOURCE = $dir"

if [ -d "$BUILD_DIR" ]; then
  echo "Build directory, $BUILD_DIR already exists."
else
  echo "Making build directory, $BUILD_DIR"
  mkdir "$BUILD_DIR"
fi

# cd into build directory and remove all previous files
cd "$BUILD_DIR"
rm -rf *

#######################################
##### RUN CMAKE
#######################################
cmake \
-DCMAKE_BUILD_TYPE=$build_type \
-DCMAKE_INSTALL_PREFIX=$INSTALL_DIR \
-DCMAKE_CXX_COMPILER=$my_cxx_compiler \
-DCMAKE_C_COMPILER=$my_cc_compiler \
-DMUQ_USE_GTEST=ON \
-DMUQ_GTEST_DIR=$GTEST_DIR \
-DMUQ_USE_PYTHON=$with_python \
-DMUQ_USE_NLOPT=$with_nlopt \
-DMUQ_NLOPT_DIR=$NLOPT_DIR \
-DMUQ_USE_OPENMPI=OFF \
-DBOOST_EXTERNAL_SOURCE=$BOOST_SOURCE \
$dir

#######################################
##### BUILD MUQ
#######################################
make install

cd $dir

exit 0
