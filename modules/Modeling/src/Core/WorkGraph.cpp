#include "MUQ/Modeling/Core/WorkGraph.h"

using namespace muq::Modeling::Core;

/// A helper class that determines if a node in the graph has a given name
class NodeNameFinder {
public:

  /**
     @param[in] name We are looking for nodes with this name
     @param[in] graph A pointer to the graph that stores the nodes
   */
  NodeNameFinder(std::string const& name, std::shared_ptr<Graph> graph) : name(name), graph(graph) {}

  /// Does a given vertex have the same name as the given name
  /**
     param[in] vertex The vertex of the graph
     \return true if the names are the same, false if not
   */
  bool operator()(boost::graph_traits<Graph>::vertex_descriptor vertex) const {
    // check the names
    return name.compare(graph->operator[](vertex)->name)==0;
  }

  /// We are looking for vertices with this name
  const std::string name;

  /// This graph stores the vertices
  std::shared_ptr<Graph> graph;
};

WorkGraph::WorkGraph() {
  // create an empty graph
  graph = std::make_shared<Graph>();
}

WorkGraph::~WorkGraph() {}

unsigned int WorkGraph::NumNodes() const {
  assert(graph);

  // return the number of vertices
  return boost::num_vertices(*graph);
}

unsigned int WorkGraph::NumEdges() const {
  assert(graph);

  // return the number of edges
  return boost::num_edges(*graph);
}

bool WorkGraph::HasNode(std::string const& name) const {
  assert(graph);

  // get iterators to the begining and end of the graph
  boost::graph_traits<Graph>::vertex_iterator v, v_end;
  boost::tie(v, v_end) = vertices(*graph);

  // look for a node with this name
  const boost::graph_traits<Graph>::vertex_iterator iter = std::find_if(v, v_end, NodeNameFinder(name, graph));

  // the node exists if the iterator is not the end
  return iter!=v_end;
}

void WorkGraph::AddNode(std::shared_ptr<WorkPiece> input, std::string const& name) {
  // make sure this node does not already exist
  assert(!HasNode(name));

  // add a node to the graph
  auto node = add_vertex(*graph);

  graph->operator[](node) = std::make_shared<WorkGraphNode>(input, name);
}
