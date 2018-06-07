#include "AllClassWrappers.h"

#include "MUQ/Utilities/LinearAlgebra/AnyAlgebra.h"
#include "MUQ/Utilities/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Utilities/LinearAlgebra/BlockRowOperator.h"
#include "MUQ/Utilities/LinearAlgebra/CompanionMatrix.h"
#include "MUQ/Utilities/LinearAlgebra/ConcatenateOperator.h"
#include "MUQ/Utilities/LinearAlgebra/DiagonalOperator.h"
#include "MUQ/Utilities/LinearAlgebra/EigenMatrixAlgebra.h"
#include "MUQ/Utilities/LinearAlgebra/EigenVectorAlgebra.h"
#include "MUQ/Utilities/LinearAlgebra/KroneckerProductOperator.h"
#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"
#include "MUQ/Utilities/LinearAlgebra/ProductOperator.h"
#include "MUQ/Utilities/LinearAlgebra/ScalarAlgebra.h"
#include "MUQ/Utilities/LinearAlgebra/SumOperator.h"
#include "MUQ/Utilities/LinearAlgebra/SundialsAlgebra.h"

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