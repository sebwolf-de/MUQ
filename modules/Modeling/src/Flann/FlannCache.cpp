#include "MUQ/Modeling/Flann/FlannCache.h"

using namespace muq::Modeling;

FlannCache::Entry::Entry() {}

FlannCache::FlannCache(std::shared_ptr<WorkPiece> function) : function(function) {
  // the target function can only have one input/output
  assert(function->numInputs==1);
  assert(function->numOutputs==1);
  
  algebra = std::make_shared<AnyAlgebra>();
}

FlannCache::~FlannCache() {}

void FlannCache::EvaluateImpl(ref_vector<boost::any> const& inputs) {}

void FlannCache::DeepVectorCopy(boost::any const& vec, flann::Matrix<double>& fvec) const {
  // get the dimension of the vector
  const unsigned int dim = algebra->VectorDimensionBase(vec);

  // create a flann vector of the same size
  fvec = flann::Matrix<double>(new double[dim], 1, dim);
  assert(fvec.cols==dim); assert(fvec.rows==1);

  // copy all of the elements into the flann vector
  for( unsigned int i=0; i<dim; ++i ) {
    fvec[0][i] = boost::any_cast<double const>(algebra->AccessElementBase(i, vec));
  }
}

int FlannCache::InCache(boost::any const& input) const {
  if( Size()>0 ) { // if there are points in the cache
    assert(nnIndex);
    
    // copy the vector into a flann type
    flann::Matrix<double> in;
    DeepVectorCopy(input, in);

    // find all of the neighbors in a small ball around the point
    std::vector<std::vector<int> > indices;
    std::vector<std::vector<double> > dists;
    int num = nnIndex->radiusSearch(in, indices, dists, 1.0e-10, flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
    assert(indices.size()==1); assert(dists.size()==1);
    assert(indices[0].size()==num); assert(dists[0].size()==num);

    // if one of the distances is zero---we have this point
    const int d = std::distance(dists[0].begin(), std::find(dists[0].begin(), dists[0].end(), 0.0));
    if( d<dists[0].size() ) { // there is a "distance zero" point ...
      // ... return the index 
      return indices[0][d];
    }
  }
  
  // the cache is either empty or none of the points in a small radius are exactly the point we care about
  return -1;
}

void FlannCache::Add(boost::any const& input) {
  // the input dimension
  const unsigned int dimIn = algebra->VectorDimensionBase(input);
  
  if( !nnIndex ) { // if the nearest neighbor index does not exist ...
    // ... create the nearest neighbor index
    auto nnDataset = flann::Matrix<double>(nullptr, 0, dimIn);
    nnIndex = std::make_shared<flann::Index<flann::L2<double> > >(nnDataset, flann::KDTreeSingleIndexParams());
  }
  // make sure we have a nearest neighbor index
  assert(nnIndex);

  // copy the input point
  flann::Matrix<double> in;
  DeepVectorCopy(input, in);

  // add the point to the nearest neighbor index
  nnIndex->addPoints(in);

  // evaluate the point
  const std::vector<boost::any>& result = function->Evaluate(ref_vector<boost::any>(1, input));
  assert(result.size()==1);
  
  // create a new entry in the cache
  auto entry = std::make_shared<Entry>();
  
  // copy the point as a flann type (copy the input and output into the entry)
  DeepVectorCopy(result[0], entry->output);
  
  // add the entry to the cache
  cache[nextID++] = entry;
}

void FlannCache::Remove(boost::any const& input) {
  // get the index of the point
  const int id = InCache(input);

  // the point is not in the cache ... nothing to remove
  if( id<0 ) { return; }

  // remove it from the nearest neighbor index
  nnIndex->removePoint(id);

  // remove it from the cache
  cache.erase(id);
}

void FlannCache::NearestNeighbors(boost::any const& point, unsigned int const k) const {
  // make sure we have enough
  assert(k<=Size());

  // conver the input to a flann matrix
  flann::Matrix<double> input;
  DeepVectorCopy(point, input);

  std::vector<std::vector<int> > indices;
  std::vector<std::vector<double> > dists;

  nnIndex->knnSearch(input, indices, dists, k, flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));

  std::cout << indices.size() << std::endl;
  std::cout << dists.size() << std::endl;

  std::cout << indices[0].size() << std::endl;
  std::cout << dists[0].size() << std::endl;

  for( auto ind : indices[0] ) {
    std::cout << ind << std::endl;
  }

  std::cout << std::endl;
}

unsigned int FlannCache::Size() const { return cache.size(); }
