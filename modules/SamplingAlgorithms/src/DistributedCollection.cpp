#include "MUQ/SamplingAlgorithms/DistributedCollection.h"

#if MUQ_HAS_MPI

using namespace muq::SamplingAlgorithms;

DistributedCollection::DistributedCollection(std::shared_ptr<SampleCollection> collection, std::shared_ptr<parcer::Communicator> comm) : SampleCollection(), collection(collection), comm(comm) {}

void DistributedCollection::Add(std::shared_ptr<SamplingState> newSamp) {
  collection->Add(newSamp);
}

#endif // end MUQ_HAS_MPI
