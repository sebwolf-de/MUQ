#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Utilities/Exceptions.h"

#include <chrono>

using namespace muq::Modeling;

ModPiece::ModPiece(Eigen::VectorXi const& inputSizesIn,
                   Eigen::VectorXi const& outputSizesIn) : WorkPiece(std::vector<std::string>(inputSizesIn.size(), typeid(Eigen::VectorXd).name()),
                                                                     std::vector<std::string>(outputSizesIn.size(), typeid(Eigen::VectorXd).name())),
                                                           inputSizes(inputSizesIn),
                                                           outputSizes(outputSizesIn){};


std::vector<Eigen::VectorXd> ModPiece::Evaluate(std::vector<Eigen::VectorXd> const& input)
{
  return Evaluate(ToRefVector(input));
}

std::vector<Eigen::VectorXd> ModPiece::Evaluate(ref_vector<Eigen::VectorXd> const& input)
{
  CheckInputs(input);

  numEvalCalls++;
  auto start_time = std::chrono::high_resolution_clock::now();

  EvaluateImpl(input);

  auto end_time = std::chrono::high_resolution_clock::now();
  evalTime += 1e6*static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

  return outputs;
}

Eigen::VectorXd ModPiece::Gradient(unsigned int                 const  outputDimWrt,
                                   unsigned int                 const  inputDimWrt,
                                   std::vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd              const& sensitivity)
{
  return Gradient(outputDimWrt, inputDimWrt, ToRefVector(input), sensitivity);
}


Eigen::VectorXd ModPiece::Gradient(unsigned int                const  outputDimWrt,
                                   unsigned int                const  inputDimWrt,
                                   ref_vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd             const& sensitivity)
{
  CheckInputs(input);

  numGradCalls++;
  auto start_time = std::chrono::high_resolution_clock::now();

  GradientImpl(outputDimWrt, inputDimWrt, input, sensitivity);

  auto end_time = std::chrono::high_resolution_clock::now();
  gradTime += 1e6*static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

  return gradient;
}


Eigen::VectorXd ModPiece::GradientByFD(unsigned int                 const  outputDimWrt,
                                       unsigned int                 const  inputDimWrt,
                                       std::vector<Eigen::VectorXd> const& input,
                                       Eigen::VectorXd              const& sensitivity)
{
  return GradientByFD(outputDimWrt, inputDimWrt, ToRefVector(input), sensitivity);
}

Eigen::VectorXd ModPiece::GradientByFD(unsigned int                const  outputDimWrt,
                                       unsigned int                const  inputDimWrt,
                                       ref_vector<Eigen::VectorXd> const& input,
                                       Eigen::VectorXd             const& sensitivity)
{
  return JacobianByFD(outputDimWrt,inputDimWrt, input).transpose()*sensitivity;
}

Eigen::MatrixXd ModPiece::Jacobian(unsigned int                 const  outputDimWrt,
                                   unsigned int                 const  inputDimWrt,
                                   std::vector<Eigen::VectorXd> const& input)
{
  return Jacobian(outputDimWrt, inputDimWrt, ToRefVector(input));
}

Eigen::MatrixXd ModPiece::Jacobian(unsigned int                const  outputDimWrt,
                                   unsigned int                const  inputDimWrt,
                                   ref_vector<Eigen::VectorXd> const& input)
{
  CheckInputs(input);

  numJacCalls++;
  auto start_time = std::chrono::high_resolution_clock::now();

  JacobianImpl(outputDimWrt, inputDimWrt, input);

  auto end_time = std::chrono::high_resolution_clock::now();
  jacTime += 1e6*static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

  return jacobian;
}


Eigen::VectorXd ModPiece::ApplyJacobian(unsigned int                 const  outputDimWrt,
                                        unsigned int                 const  inputDimWrt,
                                        std::vector<Eigen::VectorXd> const& input,
                                        Eigen::VectorXd              const& vec)
{
  return ApplyJacobian(outputDimWrt, inputDimWrt, ToRefVector(input), vec);
}

Eigen::VectorXd ModPiece::ApplyJacobian(unsigned int                const  outputDimWrt,
                                        unsigned int                const  inputDimWrt,
                                        ref_vector<Eigen::VectorXd> const& input,
                                        Eigen::VectorXd             const& vec)
{
  CheckInputs(input);

  numJacActCalls++;
  auto start_time = std::chrono::high_resolution_clock::now();

  ApplyJacobianImpl(outputDimWrt, inputDimWrt, input, vec);

  auto end_time = std::chrono::high_resolution_clock::now();
  jacActTime += 1e6*static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

  return jacobianAction;
}

void ModPiece::EvaluateImpl(ref_vector<boost::any> const& inputs){

  ref_vector<Eigen::VectorXd> eigenInputs;
  for(int i=0; i<inputs.size(); ++i)
    eigenInputs.push_back( std::cref(boost::any_cast<Eigen::VectorXd const&>(inputs.at(i))) );

  EvaluateImpl(eigenInputs);
}


// Eigen::VectorXd ModPiece::ApplyHessian(int                         const  outputDimWrt,
//                                        int                         const  inputDimWrt1,
//                                        int                         const  inputDimWrt2,
//                                        ref_vector<Eigen::VectorXd> const& input,
//                                        Eigen::VectorXd             const& sensitivity,
//                                        Eigen::VectorXd             const& vec)
// {
//  CheckInputs(input);
//
//  numHessCalls++;
//  auto start_time = std::chrono::high_resolution_clock::now();
//
//  ApplyHessianImpl(outputDimWrt, inputDimWrt1, inputDimWrt2, input, sensitivity, vec);
//
//  auto end_time = std::chrono::high_resolution_clock::now();
//  hessTime += 1e6*static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());
// }


void ModPiece::CheckInputs(ref_vector<Eigen::VectorXd> const& input)
{
  bool errorOccured = false;
  std::string msg = "Error evaluating ModPiece:\n";

  if(input.size() != inputSizes.size())
    msg += "  - Wrong number of input arguments.  Expected " + std::to_string(inputSizes.size()) + " inputs, but " + std::to_string(input.size()) + " were given.\n";

  for(int i=0; i<std::min<int>(input.size(), inputSizes.size()); ++i){
    if(input.at(i).get().size() != inputSizes(i))
      msg += "  - Input " + std::to_string(i) + " has the wrong size.  Expected size " + std::to_string(inputSizes(i)) + ", but given input with size " + std::to_string(input.at(i).get().size()) + "\n";
  }

  throw muq::WrongSizeError(msg);
}


void ModPiece::GradientImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input,
                            Eigen::VectorXd             const& sensitivity)
{
  gradient = Jacobian(outputDimWrt, inputDimWrt, input).transpose() * sensitivity;
}

void ModPiece::JacobianImpl(unsigned int                const  outputDimWrt,
                            unsigned int                const  inputDimWrt,
                            ref_vector<Eigen::VectorXd> const& input)
{
  jacobian = JacobianByFD(outputDimWrt, inputDimWrt, input);
}

void ModPiece::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& vec)
{
  jacobianAction = ApplyJacobianByFD(outputDimWrt, inputDimWrt, input, vec);
}

// void ModPiece::ApplyHessianImpl(int                         const  outputDimWrt,
//                                 int                         const  inputDimWrt1,
//                                 int                         const  inputDimWrt2,
//                                 ref_vector<Eigen::VectorXd> const& input,
//                                 Eigen::VectorXd             const& sensitivity
//                                 Eigen::VectorXd             const& vec)
// {
//   hessianAction = ApplyHessianByFD(outputDimWrt, inputDimWrt1, inputDimWrt2, input, sensitvity, vec);
// }

Eigen::MatrixXd ModPiece::JacobianByFD(unsigned int                 const  outputDimWrt,
                                       unsigned int                 const  inputDimWrt,
                                       std::vector<Eigen::VectorXd> const& input)
{
  return JacobianByFD(outputDimWrt, inputDimWrt, ToRefVector(input));
}

Eigen::MatrixXd ModPiece::JacobianByFD(unsigned int                const  outputDimWrt,
                                       unsigned int                const  inputDimWrt,
                                       ref_vector<Eigen::VectorXd> const& input)
{
  Eigen::VectorXd f0 = Evaluate(input).at(outputDimWrt);
  Eigen::VectorXd f;

  double eps;
  Eigen::VectorXd newInput(input.at(inputDimWrt));
  ref_vector<Eigen::VectorXd> newInputVec = input;

  Eigen::MatrixXd jac(outputSizes(outputDimWrt), inputSizes(inputDimWrt));
  for(int i=0; i<inputSizes(inputDimWrt); ++i){
    eps  = std::max(1e-8, 1e-10*std::abs(input.at(inputDimWrt)(i)));

    newInput(i) = input.at(inputDimWrt)(i) + eps;
    newInputVec.at(inputDimWrt) = std::cref(newInput);

    f = Evaluate(newInputVec).at(outputDimWrt);

    jac.col(i) = (f-f0)/eps;

    newInput(i) = input.at(inputDimWrt)(i);
  }

  return jac;
}

Eigen::VectorXd ModPiece::ApplyJacobianByFD(unsigned int                 const  outputDimWrt,
                                            unsigned int                 const  inputDimWrt,
                                            std::vector<Eigen::VectorXd> const& input,
                                            Eigen::VectorXd              const& vec)
{
  return ApplyJacobianByFD(outputDimWrt, inputDimWrt, ToRefVector(input), vec);
}
Eigen::VectorXd ModPiece::ApplyJacobianByFD(unsigned int                const  outputDimWrt,
                                            unsigned int                const  inputDimWrt,
                                            ref_vector<Eigen::VectorXd> const& input,
                                            Eigen::VectorXd             const& vec)
{
  const double eps = std::max(1e-8, 1e-10*vec.norm());

  ref_vector<Eigen::VectorXd> newInputVec = input;
  Eigen::VectorXd newInput = input.at(inputDimWrt).get() + eps*vec;
  newInputVec.at(inputDimWrt) = std::cref(newInput);

  Eigen::VectorXd f0 = Evaluate(input).at(outputDimWrt);
  Eigen::VectorXd f  = Evaluate(newInputVec).at(outputDimWrt);

  return (f-f0)/eps;
}

// Eigen::VectorXd ModPiece::ApplyHessianByFD(int                         const  outputDimWrt,
//                                            int                         const  inputDimWrt1,
//                                            int                         const  inputDimWrt2,
//                                            ref_vector<Eigen::VectorXd> const& input,
//                                            Eigen::VectorXd             const& sensitivity
//                                            Eigen::VectorXd             const& vec)
// {
//   Eigen::MatrixXd hess(inputSizes(inputDimWrt1), inputSizes(inputDimWrt2));
//
//   const double eps = std::max(1e-8, 1e-10*vec.norm());
//
//   Eigen::VectorXd newInput(input.at(inputDimWrt2));
//   ref_vector<Eigen::VectorXd> newInputVec = input;
//
//   Eigen::VectorXd newInput = input.at(inputDimWrt2) + eps*vec;
//
//   Eigen::VectorXd f0 = Gradient(outputDimWrt, inputDimWrt1, input, sensitvity);
//   Eigen::VectorXd f  = Gradient(outputDimWrt, inputDimWrt1, newInput, sensitvity);
//
//   return (f-f0)/eps;
// }
