#include "MUQ/Modeling/Flann/FlannCache.h"

using namespace muq::Modeling;

FlannCache::FlannCache(std::shared_ptr<ModPiece> function) : ModPiece(function->inputSizes, function->outputSizes), // can only have one input and output
							     function(function),
							     kdTree(std::make_shared<DynamicKDTreeAdaptor<>>(function->inputSizes(0))) {

  // the target function can only have one input/output
  assert(function->numInputs==1);
  assert(function->numOutputs==1);
}

FlannCache::~FlannCache() {}

void FlannCache::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
    int cacheId = InCache(inputs.at(0));
    outputs.resize(1);
    if(cacheId < 0){
      Add(inputs.at(0));
      outputs.at(0) = outputCache.at(outputCache.size()-1);
    }else{
      outputs.at(0) = outputCache.at(cacheId);
    }
}

int FlannCache::InCache(Eigen::VectorXd const& input) const {
  if( Size()>0 ) { // if there are points in the cache
    std::vector<size_t> indices;
    std::vector<double> squaredDists;
    std::tie(indices, squaredDists) = kdTree->query(input, 1);
    
    if(squaredDists.at(0)<std::numeric_limits<double>::epsilon()){
      return indices.at(0);
    }
  }
  
  // the cache is either empty or none of the points in a small radius are exactly the point we care about
  return -1;
}

void FlannCache::Add(Eigen::VectorXd const& newPt) {
  kdTree->add(newPt);

  Eigen::VectorXd newOutput = function->Evaluate(newPt).at(0);
  outputCache.push_back(newOutput);
}

void FlannCache::Remove(Eigen::VectorXd const& input) {
  // get the index of the point
  const int id = InCache(input);

  // the point is not in the cache ... nothing to remove
  if( id<0 ) { return; }

  kdTree->m_data.erase(kdTree->m_data.begin()+id);
  kdTree->UpdateIndex();
}

void FlannCache::NearestNeighbors(Eigen::VectorXd const& point,
                                  unsigned int const k,
                                  std::vector<Eigen::VectorXd>& neighbors,
                                  std::vector<Eigen::VectorXd>& result) const {


  // make sure we have enough
  assert(k<=Size());

  std::vector<size_t> indices;
  std::vector<double> squaredDists;
  std::tie(indices, squaredDists) = kdTree->query(point, k);

  neighbors.resize(k);
  result.resize(k);
  for(int i=0; i<k; ++i){
    neighbors.at(i) = kdTree->m_data.at(i);
    result.at(i) = outputCache.at(i);
  }
}

unsigned int FlannCache::Size() const { return kdTree->m_data.size(); }

void FlannCache::Add(std::vector<Eigen::VectorXd> const& inputs) {
  for( auto it : inputs ) {
    // add the point if is not already there
    if( InCache(it)<0 ) { Add(it); }
    
    // make sure it got added
    assert(InCache(it)>=0);
  }
}
