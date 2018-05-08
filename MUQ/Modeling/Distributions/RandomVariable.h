#ifndef RANDOMVARIABLE_H
#define RANDOMVARIABLE_H

#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq{
  namespace Modeling{

    class RandomVariable : public WorkPiece{

    public:
      RandomVariable(std::shared_ptr<Distribution> distIn);

      virtual ~RandomVariable() = default;

      boost::any Sample(ref_vector<boost::any> const& inputs);
      boost::any Sample(std::vector<boost::any> const& inputs){return Sample(ToRefVector(inputs));};

      boost::any Sample();


    protected:
      std::shared_ptr<Distribution> dist;

      boost::any input0;

      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

      virtual void JacobianImpl(unsigned int           const  wrtIn,
                                unsigned int           const  wrtOut,
                                ref_vector<boost::any> const& inputs) override;

      virtual void JacobianActionImpl(unsigned int           const  wrtIn,
                                      unsigned int           const  wrtOut,
                                      boost::any             const& vec,
                                      ref_vector<boost::any> const& inputs) override;

      virtual void JacobianTransposeActionImpl(unsigned int           const  wrtIn,
                                               unsigned int           const  wrtOut,
                                               boost::any             const& vec,
                                               ref_vector<boost::any> const& inputs) override;

    private:

      template<typename ith, typename... Args>
      inline double Sample(ref_vector<boost::any>& inputs, ith const& in, Args... args) {
        const int inputNum = inputs.size();
        assert(numInputs<0 || inputNum<numInputs);
        assert(CheckInputType(inputNum, typeid(in).name()));

        const boost::any in_any(in);
        inputs.push_back(std::cref(in_any));
        return Sample(inputs, args...);
      }


      ref_vector<boost::any> CreateInputs(ref_vector<boost::any> const& oldInputs);

    }; // class RandomVariable

  } // namespace Modeling
} // namespace muq



#endif // #ifndef RANDOMVARIABLE_H
