#include "MUQ/config.h"

#if MUQ_HAS_MPI

#include <gtest/gtest.h>

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/SamplingAlgorithms/DistributedCollection.h"

using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;

class DistributedCollectionTest : public::testing::Test {
protected:

  inline virtual void SetUp() override {
    L.resize(2,2);
    L << 1.0, 0.0,
      1.0, 2.0;

    {
      // create the local collection
      auto local = std::make_shared<SampleCollection>();

      // the default communicator is MPI_COMM_WORLD
      auto comm = std::make_shared<parcer::Communicator>();
      
      // create the global collection
      collection = std::make_shared<DistributedCollection>(local, comm);
    }
    
    samps = L * RandomGenerator::GetNormal(2,numSamps);
    weights = RandomGenerator::GetUniform(numSamps);
    weights /= weights.sum();
    
    for(int i=0; i<numSamps; ++i) {
      auto state = std::make_shared<SamplingState>(Eigen::VectorXd(samps.col(i)), weights(i));
      state->meta["id"] = i;
      state->meta["x norm"] = samps.col(i).norm();
      state->meta["vec2"] = (Eigen::Vector2d)(i*Eigen::Vector2d::Ones());
      state->meta["vec3"] = (Eigen::Vector3d)(i*Eigen::Vector3d::Ones());
      state->meta["vec4"] = (Eigen::Vector4d)(i*Eigen::Vector4d::Ones());
      state->meta["vecX"] = (Eigen::VectorXd)(i*Eigen::VectorXd::Ones(5));
      
      collection->Add(state);
    }
  }
  
  inline virtual void TearDown() override {}
  
  const int numSamps = 1e4;
  Eigen::MatrixXd L;
  
  Eigen::MatrixXd samps;
  Eigen::VectorXd weights;
  std::shared_ptr<DistributedCollection> collection;
};

TEST_F(DistributedCollectionTest, Basic) {
}

#endif // end MUQ_HAS_MPI
