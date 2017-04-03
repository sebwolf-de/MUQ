#ifndef WORKGRAPH_H_
#define WORKGRAPH_H_

#include <boost/graph/adjacency_list.hpp>

#include "MUQ/Modeling/WorkGraphNode.h"
#include "MUQ/Modeling/WorkGraphEdge.h"
#include "MUQ/Modeling/WorkGraphPiece.h"

namespace muq {
  namespace Modeling { 
    // Define a directed graph type
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, std::shared_ptr<WorkGraphNode>, std::shared_ptr<WorkGraphEdge> > Graph;

    /// A graph of connected muq::Modeling::WorkPiece's 
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
	 @param[in] outputDim The output dimension of "nameFrom" that will be given to "nameTo"
	 @param[in] nameTo The name of the downstream node.
	 @param[in] inputDim The input dimension of "nameTo" that will be given the output of "nameFrom"
	 */
      void AddEdge(std::string const& nameFrom, unsigned int const outputDim, std::string const& nameTo, unsigned int const inputDim);
      
      /// Is the given node in the graph?
      /**
	 @param[in] name The name of the node 
	 \return true if the graph has a node named 'name', false otherwise
      */
      bool HasNode(std::string const& name) const;

      /// Visualize the graph
      /**
	 @param[in] filename The name of the file where the graph visualization is to be stored
       */
      void Visualize(std::string const& filename) const;

      /// Create a new graph cutting any of the nodes that do not affect the output node
      /**
	 @param[in] nameOut The output node we care about
	 \return A new graph where all of the nodes that do not affect the output node are cut
       */
      std::shared_ptr<WorkGraph> DependentCut(std::string const& nameOut) const;

      /// Create a muq::Modeling::WorkPiece whose output matches a given node
      /**
	 @param[in] node The name of the output node
	 \return A muq::Modeling::WorkPiece whose outputs are the same as the output node
       */
      std::shared_ptr<WorkPiece> CreateWorkPiece(std::string const& node) const;

      /// Check to see if a node is constant?
      /**
	 A constant node is nonrandom, has no random dependencies, and either has no inputs or all of its inputs are given by const ant upstream nodes
	 @param[in] node The name of the node
	 \return true: the node is constant, false: the node is not constant
       */
      bool Constant(std::string const& node) const;
      
    private:

      /// Check to see if a node is constant?
      /**
	 A constant node is nonrandom, has no random dependencies, and either has no inputs or all of its inputs are given by const ant upstream nodes
	 @param[in] node We want to know if this node is constant.
	 \return true: the node is constant, false: the node is not constant
       */
      bool Constant(boost::graph_traits<Graph>::vertex_descriptor node) const;

      /// Recursively go upstream from a node, copying nodes 
      /**
	 @param[in] vOld The old vertex, where we beging copying from 
	 @param[in] vNew The new vertex, a copy of the old one on the new graph
	 @param[in] newGraph All of the upstreams nodes from vOld will be copied are copied to this new graph
       */
      void RecursiveCut(const boost::graph_traits<Graph>::vertex_descriptor& vOld, const boost::graph_traits<Graph>::vertex_descriptor& vNew, std::shared_ptr<WorkGraph> newGraph) const;

      /// Is the given node in the graph?
      /**
	 @param[out] iter If the node exists, this is an iterator to it.  If the node does not exist, this is the end iterator
	 @param[in] name The name of the node 
	 \return true if the graph has a node named 'name', false otherwise
      */
      bool HasNode(boost::graph_traits<Graph>::vertex_iterator& iter, std::string const& name) const;

      /// Is there an edge between two vertices?
      /**
	 @param[in] vOut The node whose output is going to vIn
	 @param[in] vIn The node whose is recieving output from vOut
	 @param[in] inputDim The dimension of the vIn that the output is going to
	 \return true if there is an edge between vOut and vIn, false otherwise
      */
      bool HasEdge(boost::graph_traits<Graph>::vertex_descriptor const& vOut, boost::graph_traits<Graph>::vertex_descriptor const& vIn, int const inputDim) const;

      /// Get a vertex_iterator to the node with name "name"
      /**
	 @param[in] name The name of a node in the graph
	 \return An iterator to that node (returns the "end" iterator if that node does not exist
      */
      boost::graph_traits<Graph>::vertex_iterator GetNodeIterator(std::string const& name) const;

      /// Find the inputs to the graph
      /**
	 \return The inputs to the graph, each input is the vertex and the input dimension that has not be set
       */
      std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > GraphInputs() const;

      /// Find the outputs to the graph
      /**
	 \return The outputs to the graph, each input is the vertex and the output dimension that has not be set
       */
      std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > GraphOutputs() const;

      /// The directed graph that represents this muq::Modeling::Core::WorkGraph
      std::shared_ptr<Graph> graph;
    };
  } // namespace Modeling
} // namespace muq

#endif
