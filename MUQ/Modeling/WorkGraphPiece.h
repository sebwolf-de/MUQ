#ifndef WORKGRAPHPIECE_H_
#define WORKGRAPHPIECE_H_

#include <boost/graph/adjacency_list.hpp>

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/ConstantPiece.h"

#include "MUQ/Modeling/WorkGraphNode.h"
#include "MUQ/Modeling/WorkGraphEdge.h"

namespace muq {
  namespace Modeling {
    /// Define a directed graph type
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, std::shared_ptr<WorkGraphNode>, std::shared_ptr<WorkGraphEdge> > Graph;
    
    /// A muq::Modeling::WorkPiece created from a muq::Modeling::WorkGraph
    class WorkGraphPiece : public WorkPiece {
    public:

      /// Construct a muq::Modeling::WorkGraphPiece
      /**
	 Typically muq::Modeling::WorkGraphPiece's are constructed by calling muq::Modeling::WorkGraph::CreateWorkPiece
	 @param[in] graph From inputs to the output-of-interest
	 @param[in] constantPieces Pointers to the muq::Modeling::ConstantPiece's that hold the graph's inputs
	 @param[in] inTypes The input types (if known)
	 @param[in] outputNode The muq::Modeling::WorkPiece that we ultimately want to evaluate
       */
      WorkGraphPiece(std::shared_ptr<Graph> graph, std::vector<std::shared_ptr<ConstantPiece> > const& constantPieces, std::map<unsigned int, std::string> const& inTypes, std::shared_ptr<WorkPiece> outputNode);
      
    private:

      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

      /// Run order computed during construction (input->output order)
      std::deque<boost::graph_traits<Graph>::vertex_descriptor> runOrder;

      /// The WorkGraph associated with this WorkGraphPiece
      std::shared_ptr<Graph> graph;

      /// The ID of the WorkPiece corresponding to the output node
      unsigned int outputID;

      /// The muq::Modeling::ConstantPiece's that store the inputs
      std::vector<std::shared_ptr<ConstantPiece> > constantPieces;
      
    };
  } // namespace Modeling
} // namespace muq

#endif
