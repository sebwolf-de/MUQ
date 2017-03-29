#include "MUQ/Modeling/Core/WorkGraph.h"

using namespace muq::Modeling::Core;

/// A helper struct that determines if a node in the graph has a given name
struct NodeNameFinder {
public:

  /**
     @param[in] name We are looking for nodes with this name
     @param[in] graph A pointer to the graph that stores the nodes
   */
  inline NodeNameFinder(std::string const& name, std::shared_ptr<Graph> graph) : name(name), graph(graph) {}

  /// Does a given vertex have the same name as the given name
  /**
     param[in] vertex The vertex of the graph
     \return true if the names are the same, false if not
   */
  inline bool operator()(boost::graph_traits<Graph>::vertex_descriptor vertex) const {
    // check the names
    return name.compare(graph->operator[](vertex)->name)==0;
  }

  /// We are looking for vertices with this name
  const std::string name;

  /// This graph stores the vertices
  std::shared_ptr<Graph> graph;
};

/// A helper struct that determines if an edge has the same input number 
struct SameDim {
  /**
     @param[in] dim The input number of the edge
     @param[in] graph A pointer to the graph that stores the edges
   */
  SameDim(int dim, std::shared_ptr<Graph> graph) : dim(dim), graph(graph) {}
  
  /**
     @param[in] edge The edge we want to compate the input number to
     \return true if the input number is the same, false otherwise
  */
  bool operator()(boost::graph_traits<Graph>::edge_descriptor edge) {
    return graph->operator[](edge)->dim==dim;
  }

  /// The input number of the edge
  int dim;

  /// This graph stores the edges
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

  // try to find the node with this name
  const boost::graph_traits<Graph>::vertex_iterator iter = GetNodeIterator(name);

  // the node exists if the iterator is not the end
  return iter!=vertices(*graph).second;
}

void WorkGraph::AddNode(std::shared_ptr<WorkPiece> input, std::string const& name) {
  // make sure this node does not already exist
  assert(!HasNode(name));

  // add a node to the graph
  auto node = add_vertex(*graph);

  graph->operator[](node) = std::make_shared<WorkGraphNode>(input, name);
}

void WorkGraph::AddEdge(std::string const& nameFrom, std::string const& nameTo, unsigned int const inputDim) {
  // get iterators to the upstream and downstream nodes (make sure they exist)
  auto itFrom = GetNodeIterator(nameFrom);
  assert(itFrom!=vertices(*graph).second);
  auto itTo = GetNodeIterator(nameTo);
  assert(itTo!=vertices(*graph).second);

  // remove any other edge going into dimension inputDim of the nameTo node
  boost::remove_in_edge_if(*itTo, SameDim(inputDim, graph), *graph);

  // try to add the new edge, if an edge already exists, notFound will be false and we need to delete the current edge first
  auto temp = boost::add_edge(*itFrom, *itTo, *graph);

  // set the edge to have the current dimension
  graph->operator[](temp.first) = std::make_shared<WorkGraphEdge>(inputDim);
}

boost::graph_traits<Graph>::vertex_iterator WorkGraph::GetNodeIterator(std::string const& name) const {
  assert(graph);

  // get iterators to the begining and end of the graph
  boost::graph_traits<Graph>::vertex_iterator v, v_end;
  boost::tie(v, v_end) = vertices(*graph);

  // return the iterator with this name (it is end if that node does not exist)
  return std::find_if(v, v_end, NodeNameFinder(name, graph));
}
