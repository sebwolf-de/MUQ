#include "AllClassWrappers.h"

#include "AnyAlgebra.h"
#include "BlockDiagonalOperator.h"
#include "BlockRowOperator.h"
#include "CompanionMatrix.h"
#include "ConcatenateOperator.h"
#include "DiagonalOperator.h"
#include "EigenMatrixAlgebra.h"
#include "EigenVectorAlgebra.h"
#include "KroneckerProductOperator.h"
#include "LinearOperator.h"
#include "ProductOperator.h"
#include "ScalarAlgebra.h"
#include "SumOperator.h"
#include "SundialsAlgebra.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

void muq::SamplingAlgorithms::PythonBindings::LinearAlgebraWrapper(py::module &m)
{
    
}