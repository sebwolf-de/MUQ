#include "MUQ/Modeling/Flann/FlannCache.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

FlannCache::FlannCache(std::shared_ptr<WorkPiece> function,
                       int                        inputDim) : WorkPiece(1, 1), // can only have one input and output
							                                                function(function),
                                                              kdTree(std::make_shared<DynamicKDTreeAdaptor<>>(inputDim)) {

  // the target function can only have one input/output
  assert(function->numInputs==1);
  assert(function->numOutputs==1);
}

FlannCache::~FlannCache() {}

void FlannCache::EvaluateImpl(ref_vector<boost::any> const& inputs)
{
    int cacheId = InCache(inputs.at(0));
    if(cacheId < 0){
      Add(*boost::any_cast<Eigen::VectorXd>(&inputs.at(0).get()));
      outputs.resize(1);
      outputs.at(0) = outputCache.at(outputCache.size()-1);
    }else{
      outputs.resize(1);
      outputs.at(0) = outputCache.at(cacheId);
    }
}

int FlannCache::InCache(boost::any const& input) const {
  if( Size()>0 ) { // if there are points in the cache

    Eigen::VectorXd const& testPt = *boost::any_cast<Eigen::VectorXd>(&input);

    std::vector<size_t> indices;
    std::vector<double> squaredDists;
    std::tie(indices, squaredDists) = kdTree->query(testPt,1);

    if(squaredDists.at(0)<1e-10){
      return indices.at(0);
    }
  }

  // the cache is either empty or none of the points in a small radius are exactly the point we care about
  return -1;
}

void FlannCache::Add(Eigen::VectorXd const& newPt) {
  kdTree->add(newPt);

  Eigen::VectorXd newOutput = boost::any_cast<Eigen::VectorXd>(function->Evaluate(newPt).at(0));
  outputCache.push_back(newOutput);
}

void FlannCache::Remove(boost::any const& input) {

  // get the index of the point
  const int id = InCache(input);

  // the point is not in the cache ... nothing to remove
  if( id<0 ) { return; }

  kdTree->m_data.erase(kdTree->m_data.begin()+id);
  kdTree = std::make_shared<DynamicKDTreeAdaptor<>>(kdTree->m_data);
}

void FlannCache::NearestNeighbors(boost::any const& point,
                                  unsigned int const k,
                                  std::vector<Eigen::VectorXd>& neighbors,
                                  std::vector<Eigen::VectorXd>& result) const {


  // make sure we have enough
  assert(k<=Size());

  Eigen::VectorXd const& testPt = *boost::any_cast<Eigen::VectorXd>(&point);

  std::vector<size_t> indices;
  std::vector<double> squaredDists;
  std::tie(indices, squaredDists) = kdTree->query(testPt,k);

  neighbors.resize(k);
  result.resize(k);
  for(int i=0; i<k; ++i){
    neighbors.at(i) = kdTree->m_data.at(i);
    result.at(i) = outputCache.at(i);
  }
}

unsigned int FlannCache::Size() const { return kdTree->m_data.size(); }
