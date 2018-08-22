#ifndef DISTRIBUTEDCOLLECTION_H_
#define DISTRIBUTEDCOLLECTION_H_

#include "MUQ/config.h"

#if MUQ_HAS_MPI

#if !MUQ_HAS_PARCER
#error
#endif

#include <parcer/Communicator.h>

#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
  namespace SamplingAlgorithms {
    class DistributedCollection : public SampleCollection {
    public:

      DistributedCollection(std::shared_ptr<SampleCollection> collection, std::shared_ptr<parcer::Communicator> comm);

      virtual ~DistributedCollection() = default;

      /// Add a sample that is stored on this processor
      /**
	 @param[in] newSamp The sample to be added
       */
      void Add(std::shared_ptr<SamplingState> newSamp);

            virtual std::shared_ptr<SamplingState> at(unsigned i);
      virtual const std::shared_ptr<SamplingState> at(unsigned i) const;


    private:

      /// The local sample collection (stored on this processor)
      std::shared_ptr<SampleCollection> collection;

      /// The communicator for this collection
      std::shared_ptr<parcer::Communicator> comm;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif // end MUQ_HAS_MPI
#endif
