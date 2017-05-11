#include "MUQ/Modeling/WorkGraphPiece.h"

#include <boost/graph/topological_sort.hpp>

using namespace muq::Modeling;

WorkGraphPiece::WorkGraphPiece(std::shared_ptr<Graph> graph, std::vector<std::shared_ptr<ConstantPiece> > const& constantPieces,  std::map<unsigned int, std::string> const& inTypes, std::shared_ptr<WorkPiece> outputPiece) : WorkPiece(inTypes, constantPieces.size(), outputPiece->OutputTypes(), outputPiece->numOutputs), graph(graph), outputID(outputPiece->ID()), constantPieces(constantPieces) {
  // build the run order
  assert(graph);
  boost::topological_sort(*graph, std::front_inserter(runOrder));
}

void WorkGraphPiece::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // get the inputs and set them to the ConstantPiece nodes
  assert(inputs.size()==constantPieces.size());
  for( unsigned int i=0; i<inputs.size(); ++i ) {
    constantPieces[i]->SetOutputs(inputs[i]);
  }

  // a map from the WorkPiece ID to its outputs
  std::map<unsigned int, ref_vector<boost::any> > valMap;

  // loop over the run order
  for( auto it : runOrder ) {
    // how many inputs does this node require?
    const int numIns = graph->operator[](it)->piece->numInputs;

    // 
    std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > inMap;

    // loop though the input nodes
    boost::graph_traits<Graph>::in_edge_iterator e, e_end;
    for( tie(e, e_end)=boost::in_edges(it, *graph); e!=e_end; ++e ) {
      // get the WorkPiece id number, the output that it supplies, and the input that receives it
      const unsigned int id = graph->operator[](boost::source(*e, *graph))->piece->ID();
      const unsigned int inNum = graph->operator[](*e)->inputDim;
      const unsigned int outNum = graph->operator[](*e)->outputDim;

      // try to find the WorkPiece in the other upstream nodes
      auto it = inMap.find(id);
      
      if( it==inMap.end() ) { // if we have not yet needed this WorkPiece ...
	// ... add it to the list and store the input/output pair
	inMap[id] = std::vector<std::pair<unsigned int, unsigned int> >(1, std::pair<unsigned int, unsigned int>(inNum, outNum));
      } else { // we have needed this WorkPiece
	// ... add the input/output pair
	inMap[id].push_back(std::pair<unsigned int, unsigned int>(inNum, outNum));
      }
    }

    // the inputs to this WorkPiece
    boost::any empty(nullptr);
    ref_vector<boost::any> ins(numIns, std::cref(empty));

    // loop through the edges again, now we know which outputs supply which inputs
    for( auto edge : inMap ) {
      // loop over the input/output pairs supplied by this input
      for( auto in_out : edge.second ) {
	ins[in_out.first] = valMap[edge.first][in_out.second];
      }
    }

    // evaluate the current map and store the value
    graph->operator[](it)->piece->Evaluate(ins);
    valMap[graph->operator[](it)->piece->ID()] = ToRefVector(graph->operator[](it)->piece->outputs);
  }

  // store the result in the output vector
  outputs.resize(valMap[outputID].size());
  for(int i=0; i<outputs.size(); ++i)
      outputs.at(i) = valMap[outputID].at(i).get();
  
}

