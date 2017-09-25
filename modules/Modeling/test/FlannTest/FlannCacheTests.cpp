#include <gtest/gtest.h>

#include "MUQ/Modeling/Flann/FlannCache.h"

using namespace muq::Modeling;

/// A model to whose input/output pairs we can store (3 input dimensions, 2 output dimensions)
class func : public WorkPiece {
public:

  inline func() :
    WorkPiece(std::vector<std::string>({typeid(Eigen::Vector3d).name()}), // inputs: Eigen::Vector3d
	      std::vector<std::string>({typeid(Eigen::Vector2d).name()})) // outputs: Eigen::Vector2d
  {}

  inline virtual ~func() {}

private:

  inline virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    const Eigen::Vector3d& in = boost::any_cast<Eigen::Vector3d const&>(inputs[0]);

    outputs.resize(1);
    outputs[0] = Eigen::Vector2d(in(0)*in(1), in(2));
  }
};


TEST(FlannCache, CreateCache) {
  // create the model whose input/output pairs we want to store
  auto f = std::make_shared<func>();

  // create a cache
  auto cache = std::make_shared<FlannCache>(f);

  // generate some random inputs
  //const 
}
