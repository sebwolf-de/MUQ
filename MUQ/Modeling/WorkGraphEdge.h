#ifndef WORKGRAPHEDGE_H_
#define WORKGRAPHEDGE_H_

namespace muq {
  namespace Modeling { 
    namespace Core {

      /// An edge in a muq::Modeling::Core::WorkGraph
      class WorkGraphEdge {
      public:

	/// Create an edge for muq::Modeling::Core::WorkGraph
	/**
	   @param[in] dim The input number
	 */
	WorkGraphEdge(unsigned int dim);
	
	/// The input number
	const unsigned int dim;

      private:
       	
      };

    } // namespace Core
  } // namespace Modeling
} // namespace muq

#endif
