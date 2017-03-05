#include "MUQ/Modeling/Core/WorkGraphNode.h"

using namespace muq::Modeling::Core;

WorkGraphNode::WorkGraphNode(std::shared_ptr<WorkPiece> piece, std::string const& name) : piece(piece), name(name) {
  assert(piece);
}
