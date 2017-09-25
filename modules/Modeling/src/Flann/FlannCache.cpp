#include "MUQ/Modeling/Flann/FlannCache.h"

#include <flann/flann.hpp>

using namespace muq::Modeling;

FlannCache::FlannCache(std::shared_ptr<WorkPiece> function) {}

FlannCache::~FlannCache() {}

void FlannCache::EvaluateImpl(ref_vector<boost::any> const& inputs) {}
