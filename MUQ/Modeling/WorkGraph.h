#ifndef WORKGRAPH_H_
#define WORKGRAPH_H_

#include <boost/graph/adjacency_list.hpp>

#include "MUQ/Modeling/WorkGraphNode.h"
#include "MUQ/Modeling/WorkGraphEdge.h"

namespace muq {
  namespace Modeling { 
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
      
      /// Add a new edge to the graph
      /**
	 @param[in] nameFrom The name of the upstream node.
	 @param[in] nameTo The name of the downstream node.
	 @param[in] inputDim The input dimension of "nameTo" that will be given the output of "nameFrom"
	 */
      void AddEdge(std::string const& nameFrom, std::string const& nameTo, unsigned int const inputDim);
      
      /// Is the given node in the graph?
      /**
	 @param[in] name The name of the node 
	 \return true if the graph has a node named 'name', false otherwise
      */
      bool HasNode(std::string const& name) const;
      
    private:

      /// Get a vertex_iterator to the node with name "name"
      /**
	 @param[in] name The name of a node in the graph
	 \return An iterator to that node (returns the "end" iterator if that node does not exist
      */
      boost::graph_traits<Graph>::vertex_iterator GetNodeIterator(std::string const& name) const;
      
      /// The directed graph that represents this muq::Modeling::Core::WorkGraph
      std::shared_ptr<Graph> graph;
    };
  } // namespace Modeling
} // namespace muq

#endif
