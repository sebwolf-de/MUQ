#ifndef DISTRIBUTEDCOLLECTION_H_
#define DISTRIBUTEDCOLLECTION_H_

#include "MUQ/config.h"

#if MUQ_HAS_MPI

#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
  namespace SamplingAlgorithms {
    class DistributedCollection : public SampleCollection {
    public:

      DistributedCollection() = default;

      virtual ~DistributedCollection() = default;
      
    private:
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif // end MUQ_HAS_MPI
#endif
