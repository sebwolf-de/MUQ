#include "MUQ/Modeling/NodeNameFinder.h"

using namespace muq::Modeling;

NodeNameFinder::NodeNameFinder(std::string const& name, std::shared_ptr<const Graph> graph) : name(name), graph(graph) {}

bool NodeNameFinder::operator()(boost::graph_traits<Graph>::vertex_descriptor vertex) const {
  // check the names
  return name.compare(graph->operator[](vertex)->name)==0;
}
