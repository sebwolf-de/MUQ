#include "MUQ/Modeling/WorkGraphPiece.h"

#include <boost/range/adaptor/reversed.hpp>

#include <boost/graph/topological_sort.hpp>

#include <Eigen/Core>

using namespace muq::Modeling;

DependentPredicate::DependentPredicate() {}

DependentPredicate::DependentPredicate(boost::graph_traits<Graph>::vertex_descriptor const& baseNode, std::shared_ptr<const Graph> graph) {
  DownstreamNodes(baseNode, graph);
}

bool DependentPredicate::operator()(const boost::graph_traits<Graph>::vertex_descriptor& node) const {
  // is the node in the vector of depends?
  return std::find(doesDepend.begin(), doesDepend.end(), node)!=doesDepend.end();
}

void DependentPredicate::DownstreamNodes(const boost::graph_traits<Graph>::vertex_descriptor& baseNode, std::shared_ptr<const Graph> graph) {
  assert(graph);
	
  // add this node to the list of downstream nodes
  doesDepend.push_back(baseNode);
  
  // loop through its dependencies
  boost::graph_traits<Graph>::out_edge_iterator e, e_end;
  boost::tie(e, e_end) = boost::out_edges(baseNode, *graph);
  for( ; e!=e_end; ++e ) {
    // recursivley add its dependendcies to the downstream node list
    DownstreamNodes(boost::target(*e, *graph), graph);
  }
}

DependentEdgePredicate::DependentEdgePredicate() {}

DependentEdgePredicate::DependentEdgePredicate(std::shared_ptr<DependentPredicate> nodePred, std::shared_ptr<const Graph> graph) : nodePred(nodePred), graph(graph) {}

bool DependentEdgePredicate::operator()(const boost::graph_traits<Graph>::edge_descriptor& edge) const {
  assert(nodePred);
  assert(graph);
  
  // check to see if the source is a downstream node
  return nodePred->operator()(source(edge, *graph));
}

WorkGraphPiece::WorkGraphPiece(std::shared_ptr<const Graph> graph, std::vector<std::shared_ptr<ConstantPiece> > const& constantPieces,  std::vector<std::string> const& inputNames, std::map<unsigned int, std::string> const& inTypes, std::shared_ptr<WorkPiece> outputPiece, std::shared_ptr<const AnyAlgebra> algebra) : WorkPiece(inTypes, constantPieces.size(), outputPiece->OutputTypes(), outputPiece->numOutputs), graph(graph), outputID(outputPiece->ID()), constantPieces(constantPieces), algebra(algebra) {
  // build the run order
  assert(graph);
  boost::topological_sort(*graph, std::front_inserter(runOrder));

  // make sure we know the number of inputs
  assert(numInputs>=0);

  // each input only needs to loop over its downstream nodes when computing derivatives
  derivRunOrders.resize(numInputs);
  filtered_graphs.resize(numInputs);

  // compute a run order for each of the inputs so we only have to loop over their downstream nodes
  assert(numInputs==inputNames.size());
  for( unsigned int i=0; i<numInputs; ++i ) { // loop through the inputs
    // get iterators to the begining and end of the graph
    boost::graph_traits<Graph>::vertex_iterator v, v_end;
    boost::tie(v, v_end) = vertices(*graph);

    // determine the downstream nodes of this input
    auto nFilt = std::make_shared<DependentPredicate>(*std::find_if(v, v_end, NodeNameFinder(inputNames[i], graph)), graph);
    auto eFilt = std::make_shared<DependentEdgePredicate>(nFilt, graph);

    // filter the graph, we only care about downstream nodes of this input
    filtered_graphs[i] = std::make_shared<boost::filtered_graph<Graph, DependentEdgePredicate, DependentPredicate> >(*graph, *eFilt, *nFilt);
     
    // build specialized run order for each input dimension
    boost::topological_sort(*filtered_graphs[i], std::back_inserter(derivRunOrders[i]));
  }
}

WorkGraphPiece::~WorkGraphPiece() {}

void WorkGraphPiece::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // set the inputs
  SetInputs(inputs);
  
  // fill the map from the WorkPiece ID to its outputs
  OutputMap();

  // store the result in the output vector
  outputs.resize(valMap[outputID].size());
  for(int i=0; i<outputs.size(); ++i) {
    outputs.at(i) = valMap[outputID].at(i).get();
  }
}

void WorkGraphPiece::JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) {
  // set the inputs
  SetInputs(inputs);

  // fill the map from the WorkPiece ID to its outputs
  OutputMap();
  
  // a map from the WorkPiece ID a map from the output number to the jacobian of that out wrt the specified input
  std::map<unsigned int, std::map<unsigned int, boost::any> > jacMap;

  // loop through each downstream node
  for( auto node : boost::adaptors::reverse(derivRunOrders[wrtIn]) ) {
    // the ID of the current node
    const unsigned int nodeID = filtered_graphs[wrtIn]->operator[](node)->piece->ID();
    
    // get the outputs of this node that depend on the specified input
    std::vector<unsigned int> requiredOuts;

    if( nodeID==outputID ) { // if it is the output node ...
      // ... the user specifies the output derivative
      requiredOuts.push_back(wrtOut);
    } else {
      // loop though the output nodes
      boost::graph_traits<FilteredGraph>::out_edge_iterator eout, eout_end;
      for( tie(eout, eout_end)=boost::out_edges(node, *filtered_graphs[wrtIn]); eout!=eout_end; ++eout ) {
	// get the output number 
	const unsigned int outNum = filtered_graphs[wrtIn]->operator[](*eout)->outputDim;
	
	// if we have not already required this output, save it
	auto it = std::find(requiredOuts.begin(), requiredOuts.end(), outNum);
	if( it==requiredOuts.end() ) {
	  requiredOuts.push_back(outNum);
	}
      }
    }

    // how many inputs does this node require?
    const int numIns = filtered_graphs[wrtIn]->operator[](node)->piece->numInputs;

    // get the inputs for this node --- the input WorkPiece ID, the output number, and the input number
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > requiredIns;
    requiredIns.reserve(numIns);

    // loop though the output nodes
    boost::graph_traits<FilteredGraph>::in_edge_iterator ein, ein_end;
    for( tie(ein, ein_end)=boost::in_edges(node, *filtered_graphs[wrtIn]); ein!=ein_end; ++ein ) {
      // get the WorkPiece id number, the output that it supplies, and the input that receives it
      const unsigned int id = graph->operator[](boost::source(*ein, *graph))->piece->ID();
      const unsigned int outNum = filtered_graphs[wrtIn]->operator[](*ein)->outputDim;
      const unsigned int inNum = filtered_graphs[wrtIn]->operator[](*ein)->inputDim;

      // store the requried input
      requiredIns.push_back(std::tuple<unsigned int, unsigned int, unsigned int>(id, outNum, inNum));
    }

    // the inputs to this WorkPiece
    const ref_vector<boost::any>& ins = Inputs(node);

    // compute the jacobian of each required output wrt each input
    for( auto out : requiredOuts ) {
      if( requiredIns.size()==0 ) {
	// if there are no inputs, it is the input so set the Jacobian to the identity
	jacMap[nodeID][out] = algebra->IdentityBase(valMap[nodeID][out]);
      } else { 
	// initize the jacobian to nothing
	jacMap[nodeID][out] = boost::none;

	for( auto in : requiredIns ) {
	  // compute the Jacobian with respect to each required input
	  graph->operator[](node)->piece->Jacobian(std::get<2>(in), out, ins);
	  
	  // use chain rule to get the jacobian wrt to the required input
	  const boost::any& tempJac = algebra->MultiplyBase(*(graph->operator[](node)->piece->jacobian), jacMap[std::get<0>(in)][std::get<1>(in)]);
	  jacMap[nodeID][out] = algebra->AddBase(jacMap[nodeID][out], tempJac);
	}
      }
    }
  }

  // set the Jacobian for this WorkPiece
  jacobian = jacMap[outputID][wrtOut];
}

void WorkGraphPiece::SetInputs(ref_vector<boost::any> const& inputs) {
  // get the inputs and set them to the ConstantPiece nodes
  assert(inputs.size()==constantPieces.size());
  for( unsigned int i=0; i<inputs.size(); ++i ) {
    constantPieces[i]->SetOutputs(inputs[i]);
  }
}

std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > WorkGraphPiece::InputNodes(boost::graph_traits<Graph>::vertex_descriptor const& node) const {
  // the map of input nodes
  std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > inMap;
  
  // loop though the input nodes
  boost::graph_traits<Graph>::in_edge_iterator e, e_end;
  for( tie(e, e_end)=boost::in_edges(node, *graph); e!=e_end; ++e ) {
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

  return inMap;
}

void WorkGraphPiece::OutputMap() {
  // clear the map from the WorkPiece ID to its outputs
  valMap.clear();
  
  // loop over the run order
  for( auto it : runOrder ) {
    // the inputs to this WorkPiece
    const ref_vector<boost::any>& ins = Inputs(it);

    // evaluate the current map and store the value
    graph->operator[](it)->piece->Evaluate(ins);
    valMap[graph->operator[](it)->piece->ID()] = ToRefVector(graph->operator[](it)->piece->outputs);
  }
}

ref_vector<boost::any> WorkGraphPiece::Inputs(boost::graph_traits<Graph>::vertex_descriptor node) const {
  // how many inputs does this node require?
  const int numIns = graph->operator[](node)->piece->numInputs;
  
  // get the inputs for this node
  const std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > >& inMap = InputNodes(node);
  
  boost::any empty = boost::none;
  ref_vector<boost::any> ins(numIns, std::cref(empty));
  
  // loop through the edges again, now we know which outputs supply which inputs
  for( auto edge : inMap ) {
    // loop over the input/output pairs supplied by this input
    for( auto in_out : edge.second ) {
      // use at instead of operator[] because this is a const function
      ins[in_out.first] = valMap.at(edge.first)[in_out.second];
    }
  }

  return ins;
}
