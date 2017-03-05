#ifndef WORKGRAPH_H_
#define WORKGRAPH_H_

#include <boost/graph/adjacency_list.hpp>

#include "MUQ/Modeling/Core/WorkGraphNode.h"
#include "MUQ/Modeling/Core/WorkGraphEdge.h"

namespace muq {
  namespace Modeling { 
    namespace Core {

      // Define a directed graph type
      typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, std::shared_ptr<WorkGraphNode>, std::shared_ptr<WorkGraphEdge> > Graph;

      class WorkGraph {
      public:
	
	WorkGraph();

	virtual ~WorkGraph();

	/// Get the number of nodes in the graph
	/**
	   \return The number of nodes
	 */
	unsigned int NumNodes() const;

	/// Get the number of edgess in the graph
	/**
	   \return The number of edges
	 */
	unsigned int NumEdges() const;

	/// Add a new node to the graph
	/**
	   @param[in] input A pointer to a muq::Modelling::Core::WorkPiece that will be called when this node is evaluated.
	   @param[in] name A string representing a <b>unique</b> name for this node.
	 */
	void AddNode(std::shared_ptr<WorkPiece> input, std::string const& name);

	/// Is the given node in the graph?
	/**
	   @param[in] name The name of the node 
	   \return true if the graph has a node named 'name', false otherwise
	 */
	bool HasNode(std::string const& name) const;

      private:
	
	/// The directed graph that represents this muq::Modeling::Core::WorkGraph
	std::shared_ptr<Graph> graph;

      };
      
    } // namespace Core
  } // namespace Modeling
} // namespace muq

#endif
