#include "MUQ/SamplingAlgorithms/DistributedCollection.h"

#include "parcer/Eigen.h"

#if MUQ_HAS_MPI

using namespace muq::SamplingAlgorithms;

DistributedCollection::DistributedCollection(std::shared_ptr<SampleCollection> collection, std::shared_ptr<parcer::Communicator> comm) : SampleCollection(), collection(collection), comm(comm) {}

void DistributedCollection::Add(std::shared_ptr<SamplingState> newSamp) {
  collection->Add(newSamp);
}

std::shared_ptr<SamplingState> DistributedCollection::at(unsigned i) { return collection->at(i); }
	    
const std::shared_ptr<SamplingState> DistributedCollection::at(unsigned i) const { return collection->at(i); }

unsigned int DistributedCollection::LocalSize() const { return collection->size(); }

unsigned int DistributedCollection::GlobalSize() const {
  int size = 0;
  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    int localSize = i==comm->GetRank() ? LocalSize() : 0;
    comm->Bcast(localSize, i);

    size += localSize;
  }

  return size;
}

unsigned int DistributedCollection::size() const { return GlobalSize(); }

Eigen::VectorXd DistributedCollection::LocalCentralMoment(unsigned order, int blockDim) const { return collection->CentralMoment(order, blockDim); }

Eigen::VectorXd DistributedCollection::GlobalCentralMoment(unsigned order, int blockDim) const { return GlobalEigenMean(LocalCentralMoment(order, blockDim)); }

Eigen::VectorXd DistributedCollection::CentralMoment(unsigned order, int blockDim) const { return GlobalCentralMoment(order, blockDim); }

Eigen::VectorXd DistributedCollection::LocalMean(int blockDim) const { return collection->Mean(blockDim); }

Eigen::VectorXd DistributedCollection::GlobalMean(int blockDim) const { return GlobalEigenMean(LocalMean(blockDim)); }

Eigen::VectorXd DistributedCollection::Mean(int blockDim) const { return GlobalMean(blockDim); }

Eigen::VectorXd DistributedCollection::LocalVariance(int blockDim) const { return collection->Variance(blockDim); }

Eigen::VectorXd DistributedCollection::GlobalVariance(int blockDim) const { return GlobalEigenMean(LocalVariance(blockDim)); }

Eigen::VectorXd DistributedCollection::Variance(int blockDim) const { return GlobalVariance(blockDim); }

Eigen::MatrixXd DistributedCollection::LocalCovariance(int blockDim) const { return collection->Covariance(blockDim); }

Eigen::MatrixXd DistributedCollection::GlobalCovariance(int blockDim) const { return GlobalEigenMean(LocalCovariance(blockDim)); }

Eigen::MatrixXd DistributedCollection::Covariance(int blockDim) const { return GlobalCovariance(blockDim); }

Eigen::VectorXd DistributedCollection::LocalESS(int blockDim) const { return collection->ESS(blockDim); }

Eigen::VectorXd DistributedCollection::GlobalESS(int blockDim) const { return GlobalEigenMean(LocalESS(blockDim)); }

Eigen::VectorXd DistributedCollection::ESS(int blockDim) const { return GlobalESS(blockDim); }

Eigen::MatrixXd DistributedCollection::AsLocalMatrix(int blockDim) const { return collection->AsMatrix(blockDim); }

Eigen::MatrixXd DistributedCollection::AsGlobalMatrix(int blockDim) const {
  const Eigen::MatrixXd& local = AsLocalMatrix(blockDim);
  Eigen::MatrixXd global(local.rows(), GlobalSize());

  int numSamps = 0;
  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    int localSize = i==comm->GetRank() ? LocalSize() : 0;
    comm->Bcast(localSize, i);

    Eigen::MatrixXd l(local.rows(), localSize);
    if( comm->GetRank()==i ) { l = local; }

    comm->Bcast(l, i);
    global.block(0, numSamps, local.rows(), localSize) = l;

    numSamps += localSize;
  }

  return global;
}

Eigen::MatrixXd DistributedCollection::AsMatrix(int blockDim) const { return AsGlobalMatrix(blockDim); }

Eigen::VectorXd DistributedCollection::LocalWeights() const { return collection->Weights(); }

Eigen::VectorXd DistributedCollection::GlobalWeights() const {
  const Eigen::VectorXd& local = LocalWeights();
  Eigen::VectorXd global = Eigen::VectorXd::Constant(GlobalSize(), std::numeric_limits<double>::quiet_NaN());

  int numSamps = 0;
  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    int localSize = i==comm->GetRank() ? LocalSize() : 0;
    comm->Bcast(localSize, i);

    Eigen::VectorXd l(localSize);
    if( comm->GetRank()==i ) { l = local; }

    comm->Bcast(l, i);
    global.segment(numSamps, localSize) = l;

    numSamps += localSize;
  }

  return global;
}

Eigen::VectorXd DistributedCollection::Weights() const { return GlobalWeights(); }

void DistributedCollection::WriteToFile(std::string const& filename, std::string const& dataset) const {
  const unsigned int size = GlobalSize();
  int numSamps = 0;
  for( unsigned int i=0; i<comm->GetSize(); ++i ) {
    int localSize = i==comm->GetRank() ? LocalSize() : 0;
    comm->Bcast(localSize, i);

    if( comm->GetRank()==i ) { collection->WriteToFile(numSamps, filename, dataset, size); }
    numSamps += localSize;

    comm->Barrier();
  }
}

#endif // end MUQ_HAS_MPI
