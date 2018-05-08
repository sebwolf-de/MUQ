#include "MUQ/Modeling/Distributions/RandomVariable.h"

using namespace muq::Modeling;


RandomVariable::RandomVariable(std::shared_ptr<Distribution> distIn) : WorkPiece(), dist(distIn), input0(Distribution::Mode::SampleDistribution)
{
  assert(dist);
  numInputs = std::max(dist->numInputs-1, -1);

  for(auto& it : dist->inputTypes){
    if(it.first > 1)
      inputTypes[it.first-1] = it.second;
  }
}

ref_vector<boost::any> RandomVariable::CreateInputs(ref_vector<boost::any> const& oldInputs)
{
  ref_vector<boost::any> newInputs(1, std::cref(input0));
  newInputs.insert(newInputs.end(), oldInputs.begin(),oldInputs.end());
  return newInputs;
}

boost::any RandomVariable::Sample()
{
  return dist->Sample();
}

boost::any RandomVariable::Sample(ref_vector<boost::any> const& inputs)
{
  return dist->Sample(inputs);
}


void RandomVariable::EvaluateImpl(ref_vector<boost::any> const& inputs)
{
  outputs = dist->Evaluate(CreateInputs(inputs));
}


void RandomVariable::JacobianImpl(unsigned int           const  wrtIn,
                                  unsigned int           const  wrtOut,
                                  ref_vector<boost::any> const& inputs)
{
  jacobian = dist->Jacobian(wrtIn+1,wrtOut, CreateInputs(inputs));
}

void RandomVariable::JacobianActionImpl(unsigned int           const  wrtIn,
                                        unsigned int           const  wrtOut,
                                        boost::any             const& vec,
                                        ref_vector<boost::any> const& inputs)
{
  jacobianAction = dist->JacobianAction(wrtIn+1,wrtOut, vec, CreateInputs(inputs));
}

void RandomVariable::JacobianTransposeActionImpl(unsigned int           const  wrtIn,
                                                 unsigned int           const  wrtOut,
                                                 boost::any             const& vec,
                                                 ref_vector<boost::any> const& inputs)
{
  jacobianTransposeAction = dist->JacobianTransposeAction(wrtIn+1, wrtOut, vec, CreateInputs(inputs));
}
