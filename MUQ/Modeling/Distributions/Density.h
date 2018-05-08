#ifndef DENSITY_H
#define DENSITY_H

#include "MUQ/Modeling/Distributions/Distribution.h"

namespace muq{
  namespace Modeling{

    class Density : public WorkPiece{

    public:
      Density(std::shared_ptr<Distribution> distIn);

      virtual ~Density() = default;

      virtual double LogDensity(ref_vector<boost::any> const& inputs);
      virtual double LogDensity(std::vector<boost::any> const& inputs){return LogDensity(ToRefVector(inputs));};

      template<typename... Args>
	    inline double LogDensity(Args... args) {
	      ref_vector<boost::any> inputs;
	      inputs.reserve( (numInputs<0) ? 0 : numInputs-1);
	      return LogDensity(inputs, args...);
      }


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
      inline double LogDensity(ref_vector<boost::any>& inputs, ith const& in, Args... args) {
        const int inputNum = inputs.size();
        assert(numInputs<0 || inputNum<numInputs);
        assert(CheckInputType(inputNum, typeid(in).name()));

        const boost::any in_any(in);
        inputs.push_back(std::cref(in_any));
        return LogDensity(inputs, args...);
      }


      ref_vector<boost::any> CreateInputs(ref_vector<boost::any> const& oldInputs);

    }; // class Density

  } // namespace Modeling
} // namespace muq



#endif // #ifndef DENSITY_H
