#include "MUQ/Modeling/WorkGraph.h"

#include <fstream>
#include <algorithm>

#include <boost/algorithm/string.hpp>

// boost graph library includes
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>

using namespace muq::Modeling;

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
struct SameInputDim {
  /**
     @param[in] dim The input number of the edge
     @param[in] graph A pointer to the graph that stores the edges
   */
  SameInputDim(unsigned int const dim, std::shared_ptr<Graph> graph) : dim(dim), graph(graph) {}
  
  /**
     @param[in] edge The edge we want to compate the input number to
     \return true if the input number is the same, false otherwise
  */
  bool operator()(boost::graph_traits<Graph>::edge_descriptor edge) {
    return graph->operator[](edge)->inputDim==dim;
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

void WorkGraph::AddEdge(std::string const& nameFrom, unsigned int const outputDim, std::string const& nameTo, unsigned int const inputDim) {
  // get iterators to the upstream and downstream nodes (make sure they exist)
  auto itFrom = GetNodeIterator(nameFrom);
  assert(itFrom!=vertices(*graph).second);
  auto itTo = GetNodeIterator(nameTo);
  assert(itTo!=vertices(*graph).second);

  // either we don't know the number of outputs from "nameFrom" or the output dimension is less than the number of outputs
  if( graph->operator[](*itFrom)->piece->numOutputs>=0 && outputDim>=graph->operator[](*itFrom)->piece->numOutputs ) {
    std::cerr << std::endl << "ERROR: The number of outputs for node '" << nameFrom << "' is " << graph->operator[](*itFrom)->piece->numOutputs << " but the output required by 'WorkGraph::AddEdge' is " << outputDim << std::endl << std::endl;
    assert(graph->operator[](*itFrom)->piece->numOutputs<0 || outputDim<graph->operator[](*itFrom)->piece->numOutputs);
  }
  
  // either we don't know the number of inputs to "nameTo" or the input dimension is less than the number of inputs
  if( graph->operator[](*itTo)->piece->numInputs>=0 && inputDim>=graph->operator[](*itTo)->piece->numInputs ) {
    std::cerr << std::endl << "ERROR: The number of inputs for node '" << nameTo << "' is " << graph->operator[](*itTo)->piece->numInputs << " but the input required by 'WorkGraph::AddEdge' is " << inputDim << std::endl << std::endl;
    assert(graph->operator[](*itTo)->piece->numInputs<0 || inputDim<graph->operator[](*itTo)->piece->numInputs);
  }

  // either we don't know the input and/or output type or they match
  if(graph->operator[](*itTo)->piece->InputType(inputDim).compare("")!=0 && // we don't know the input type
     graph->operator[](*itFrom)->piece->OutputType(outputDim).compare("")!=0 && // we don't know the output type
     graph->operator[](*itTo)->piece->InputType(inputDim).compare(graph->operator[](*itFrom)->piece->OutputType(outputDim))!=0 ) { // the types must match
    std::cerr << std::endl << "ERROR: Types do not match in 'WorkGraph::AddEdge'.  The input type node '" << nameTo << "' is " << graph->operator[](*itTo)->piece->InputType(inputDim) << " but the output type for node '" << nameFrom << "' is " << graph->operator[](*itFrom)->piece->OutputType(outputDim) << std::endl << std::endl;
    assert(graph->operator[](*itTo)->piece->InputType(inputDim).compare(graph->operator[](*itFrom)->piece->OutputType(outputDim))==0); // the types must match
  }
    
  // remove any other edge going into dimension inputDim of the nameTo node
  boost::remove_in_edge_if(*itTo, SameInputDim(inputDim, graph), *graph);

  // try to add the new edge, if an edge already exists, notFound will be false and we need to delete the current edge first
  auto temp = boost::add_edge(*itFrom, *itTo, *graph);

  // set the edge to have the current dimension
  graph->operator[](temp.first) = std::make_shared<WorkGraphEdge>(outputDim, inputDim);
}

boost::graph_traits<Graph>::vertex_iterator WorkGraph::GetNodeIterator(std::string const& name) const {
  assert(graph);

  // get iterators to the begining and end of the graph
  boost::graph_traits<Graph>::vertex_iterator v, v_end;
  boost::tie(v, v_end) = vertices(*graph);

  // return the iterator with this name (it is end if that node does not exist)
  return std::find_if(v, v_end, NodeNameFinder(name, graph));
}

std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > WorkGraph::GraphOutputs() const {
  // create an empty vector to hold outputs
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > outputs;

    // loop through the vertices
  boost::graph_traits<Graph>::vertex_iterator v, v_end;
  for( std::tie(v, v_end)=vertices(*graph); v!=v_end; ++v ) {
    // a vector of the outputs that are set
    std::vector<int> isSet;

    // number of outputs (negative indcates we don't know)
    const int numOutputs = graph->operator[](*v)->piece->numOutputs;

    // if possible, reserve memory for the outputs that are set (the size reserved is the total number of outputs, however, if it is negative we don't know how many outputs there are so we reserve whatever the compiler did by default...)
    isSet.reserve(std::max((int)isSet.capacity(), numOutputs));

    // for each vertex, loop over the input nodes and figure out if the outputs are set
    boost::graph_traits<Graph>::out_edge_iterator e, e_end;
    for( tie(e, e_end)=out_edges(*v, *graph); e!=e_end; ++e ) {
      // we have an input, so store it in the vector
      isSet.push_back(graph->operator[](*e)->outputDim);
    }

    // if an input to this ModPiece is not set, it will be stored as an output to the graph
    for( int i=0; i<numOutputs; ++i ) {
      if( std::find(isSet.begin(), isSet.end(), i)==isSet.end() ) { // if the output is not set ..
	// ... store this vertex and the output number
	outputs.push_back(std::make_pair(*v, i));
      }
    }

    // if we don't know the number of outputs
    if( numOutputs<0 ) {
      outputs.push_back(std::make_pair(*v, -1));

      if( isSet.size()>0 ) { // if some outputs have been set ...
	// the maximum output number that is set
	const unsigned int maxOut = *std::max_element(isSet.begin(), isSet.end());

	// loop through all the outputs less than the max input
	for( unsigned int i=0; i<maxOut; ++i ) {
	  if( std::find(isSet.begin(), isSet.end(), i)==isSet.end() ) { // if an output is not set ...
	    // ... it must be a graph output
	    outputs.push_back(std::make_pair(*v, i));
	  }
	}
      }
    }
  }

  return outputs;
}

std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > WorkGraph::GraphInputs() const {
  // create an empty vector to hold inputs
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > inputs;

  // loop through the vertices
  boost::graph_traits<Graph>::vertex_iterator v, v_end;
  for( std::tie(v, v_end)=vertices(*graph); v!=v_end; ++v ) {
    // a vector of the inputs that are set
    std::vector<int> isSet;

    // number of inputs (negative indcates we don't know)
    const int numInputs = graph->operator[](*v)->piece->numInputs;

    // if possible, reserve memory for the inputs that are set (the size reserved is the total number of inputs, however, if it is negative we don't know how many inputs there are so we reserve whatever the compiler did by default...)
    isSet.reserve(std::max((int)isSet.capacity(), numInputs));
    
    // for each vertex, loop over the input nodes and figure out if the inputs are set
    boost::graph_traits<Graph>::in_edge_iterator e, e_end;
    for( tie(e, e_end)=in_edges(*v, *graph); e!=e_end; ++e ) {
      // we have an input, so store it in the vector
      isSet.push_back(graph->operator[](*e)->inputDim);
    }
    
    // if an input to this ModPiece is not set, it will be stored as an input to the graph
    for( int i=0; i<numInputs; ++i ) {
      if( std::find(std::begin(isSet), std::end(isSet), i)==isSet.end() ) { // if the input is not set ..
	// ... store this vertex and the input number
	inputs.push_back(std::make_pair(*v, i));
      }
    }

    // if we don't know the number of inputs
    if( numInputs<0 ) {
      inputs.push_back(std::make_pair(*v, -1));

      if( isSet.size()>0 ) { // if some inputs have been set ...
	// the maximum input number that is set
	const unsigned int maxIn = *std::max_element(isSet.begin(), isSet.end());
	
	// loop through all the inputs less than the max input
	for( unsigned int i=0; i<maxIn; ++i ) {
	  if( std::find(isSet.begin(), isSet.end(), i)==isSet.end() ) { // if an input is not set ...
	    // ... it must be a graph input
	    inputs.push_back(std::make_pair(*v, i));
	  }
	}
      }
    }
  }

  return inputs;
}

class MyVertexWriter {
public:
  
  MyVertexWriter(std::shared_ptr<const Graph> graph) : graph(graph) {}

  void operator()(std::ostream& out, const boost::graph_traits<Graph>::vertex_descriptor& v) const {
    int status;

    // the WorkPiece associated with this node
    auto workPtr = graph->operator[](v)->piece;

    // the name of this work piece
    const std::string nodeName = workPtr->Name();

    // style for the node visualization 
    const std::string style = "colorscheme=pastel16,color=2, style=filled";

    // label the node
    out << "[label=\"" << graph->operator[](v)->name << " : " << nodeName << "\", " << style << "]";
  }

private:
  // the graph we are visualizing
  std::shared_ptr<const Graph> graph;
};

class MyEdgeWriter {
public:

  MyEdgeWriter(std::shared_ptr<const Graph> graph) : graph(graph) {}

  void operator()(std::ostream& out, const boost::graph_traits<Graph>::edge_descriptor& e) const {
    const unsigned int inputDim = graph->operator[](e)->inputDim;
    
    const unsigned int outputDim = graph->operator[](e)->outputDim;

    // first, write the name as the label
    out << "[label=\" [out, in]: [" << outputDim << ", " << inputDim << "]\"]";
  }

private:
  // the graph we are visualizing
  std::shared_ptr<const Graph> graph;
};

class MyGraphWriter {
public:

  MyGraphWriter(std::shared_ptr<const Graph> graph) : graph(graph) {}

  void operator()(std::ostream& out) const {
    out << "splines = true;" << std::endl;
  }

private:
  std::shared_ptr<const Graph> graph;
};

void WorkGraph::Visualize(std::string const& filename) const {
  // split the graph (name and extension)
  std::vector<std::string> strs;
  boost::split(strs, filename, boost::is_any_of("."));

  // is the extension something we expect?
  const bool knownExtension = (strs.end()-1)->compare("png") || (strs.end()-1)->compare("jpg") || (strs.end()-1)->compare("tif") || (strs.end()-1)->compare("eps") || (strs.end()-1)->compare("pdf") || (strs.end()-1)->compare("svg");

  // open a file stream
  std::ofstream fout;

  const std::string tempname = *strs.begin() + "_temp.dot";
  if( knownExtension ) { // if we know the file extension ...
    // ... open a temporary file 
    fout.open(tempname);
  } else { // otherwise ... 
    // ... just open the file as the user asked
    fout.open(filename.c_str());
  }

  typedef std::map<boost::graph_traits<Graph>::vertex_descriptor, size_t> IndexMap;
  IndexMap mapIndex;
  boost::associative_property_map<IndexMap> propmapIndex(mapIndex);

  int vertexNum = 0;
  BGL_FORALL_VERTICES(v, *graph, Graph) {
    put(propmapIndex, v, vertexNum++);
  }

  boost::write_graphviz(fout, *graph, MyVertexWriter(graph), MyEdgeWriter(graph), MyGraphWriter(graph), propmapIndex);

  fout.seekp(-2, std::ios_base::cur); //back up so the rest is inside the brackets

  // loop over all the inputs and draw them
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > graphInputs = GraphInputs();

  int in = 0;
  for( auto aPair : graphInputs ) { // loop through the graph inputs
    vertexNum++;
    if( aPair.second<0 ) { // if we do not know the input number 
      fout << vertexNum << "[label=\"Unfixed input\", shape=invhouse,colorscheme=pastel13,color=1, style=filled];" << std::endl;
      fout << vertexNum << "->" << propmapIndex[aPair.first] << std::endl;
    } else { // if we know the input number
      fout << vertexNum << "[label=\"Input #" << in << "\", shape=invhouse,colorscheme=pastel13,color=1, style=filled];" << std::endl;
      fout << vertexNum << "->" << propmapIndex[aPair.first] << "[label=\" in: " << aPair.second << "\"];" << std::endl;
      ++in;
    }
  }

  // loop over all the outputs and draw them
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > graphOutputs = GraphOutputs();

  int out = 0;
  for( auto aPair : graphOutputs ) { // loop through the graph outputs
    vertexNum++;
    if( aPair.second<0 ) { // if we do not know the output number 
      fout << vertexNum << "[label=\"Unfixed output\", shape=box,colorscheme=pastel16,color=1, style=filled];" << std::endl;
      fout << propmapIndex[aPair.first]  << "->" << vertexNum << std::endl;
    } else { // if we know the input number
      fout << vertexNum << "[label=\"Output #" << out << "\", shape=box,colorscheme=pastel16,color=1, style=filled];" << std::endl;
      fout << propmapIndex[aPair.first]  << "->" << vertexNum << "[label=\" out: " << aPair.second << "\"];" << std::endl;
      ++out;
    }
  }

  fout << "}" << std::endl;

  // close the file stream
  fout.close();
  
  if( knownExtension ) {
    // move the *.dot file into the *.[whatever extension file]
    std::system(("dot -T" + *(strs.end() - 1) + " " + tempname + " -o " + filename).c_str());
    
    // remove the temporary *.dot file
    std::system(("rm "+tempname).c_str());
  }
}
