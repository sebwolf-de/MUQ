#include<iostream>
#include<vector>
#include<string>
#include<sstream>
#include<cassert>
#include<map>

#include "boost/variant.hpp"

#define WorkParent(name, numIns, inouts...)				\
  namespace muq {							\
    namespace Modelling {						\
      /* c */								\
      class name : public WorkPiece<inouts>  {				\
      public:								\
      name() :								\
	WorkPiece(numIns, #inouts) {}					\
									\
	virtual ~name() {}						\
									\
      private:								\
      };								\
    }									\
  }						

namespace muq {
  namespace Modelling {
    /// Base class for MUQ's modeling envronment
    template<typename... inouts> // template the class on its input/output types
      class WorkPiece {
    public:
      /**
	 @param[in] numIns The number of inputs
      */
    WorkPiece(unsigned int const numIns, std::string const& ios) : 
      types(InputOutputTypes(ios)), numInputs(numIns) {
	// make sure we have enough in/outputs
	assert(numIns<types.size());
      }
      
      virtual ~WorkPiece() {}
      
      template<typename first, typename... Args>			
	void Evaluate(first const& in1, Args... args) {		
	std::cout << in1 << std::endl;				
	/*return std::tuple<inouts> (100);			*/	
      }								
      
      /// Get the number of inputs
      /**
	 \return The number of inputs
       */
      unsigned int NumInputs() const { return numInputs; }

      /// Get the number of outputs
      /**
	 \return The number of outputs
       */
      unsigned int NumOutputs() const { return types.size()-numInputs; }

      /// Get the input types
      /**
	 \return A vector in input types
       */
      std::vector<std::string> InputTypes() const { return std::vector<std::string>(types.begin()+NumOutputs(), types.end()); }

      /// Get the output types
      /**
	 \return A vector in output types
       */
      std::vector<std::string> OutputTypes() const { return std::vector<std::string>(types.begin(), types.begin()+NumOutputs()); }
     
    private:
      
      /// Convert a comma-separated list of inputs into a vector
      /**
	 @param[in] ios A comma-separated list of in/outputs
	 \return a vector of in/outputs
       */
      static std::vector<std::string> InputOutputTypes(std::string ios) {
	// remove any whitespace from the in/output string
	ios.erase(remove_if(ios.begin(), ios.end(), isspace), ios.end());

	// allocate memory for a vector of in/outputs
	std::vector<std::string> iosvec(std::count(ios.begin(), ios.end(), ',')+1);

	// convert to stringstream
	std::stringstream ss(ios);

	// get each in/output time and put it into the vector
	for( auto it=iosvec.begin(); it!=iosvec.end(); ++it ) {
	  // get the next in/output
	  std::string substr;
	  std::getline(ss, substr, ','); 

	  // store it in the vector
	  *it = substr;
	}

	// return the in/outputs
	return iosvec;
      }

      /// A vector of in/output types
      std::vector<std::string> types;

      /// The number of inputs
      const unsigned int numInputs;

      /// A map that stores the inputs 
      std::map<unsigned int, boost::variant<inouts...> > inputs;

    }; // class WorkPiece
  } // namespace Modelling
} // namespace muq
