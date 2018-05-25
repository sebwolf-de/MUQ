#include "gtest/gtest.h"

#include "MUQ/Modeling/WorkGraph.h"

#include "WorkPieceTestClasses.h"

#include "MUQ/Modeling/ConstantPiece.h"
#include "MUQ/Modeling/AnyAlgebra.h"

using namespace muq::Modeling;

TEST(WorkGraphTests, UnfixedInOut) {
  // create test WorkPieces
  auto test0 = std::make_shared<UnfixedMod>();
  auto test1 = std::make_shared<UnfixedMod>();
  auto test2 = std::make_shared<UnfixedMod>();

  // create and empty graph
  auto graph = std::make_shared<WorkGraph>();

  // make sure the graph is empty
  EXPECT_EQ(graph->NumNodes(), 0);
  EXPECT_EQ(graph->NumEdges(), 0);

  // add WorkPieces to the graph
  EXPECT_FALSE(graph->HasNode("test 0"));
  graph->AddNode(test0, "test 0");
  EXPECT_EQ(graph->NumNodes(), 1);
  EXPECT_TRUE(graph->HasNode("test 0"));

  EXPECT_FALSE(graph->HasNode("test 1"));
  graph->AddNode(test1, "test 1");
  EXPECT_EQ(graph->NumNodes(), 2);
  EXPECT_TRUE(graph->HasNode("test 1"));

  EXPECT_FALSE(graph->HasNode("test 2"));
  graph->AddNode(test2, "test 2");
  EXPECT_EQ(graph->NumNodes(), 3);
  EXPECT_TRUE(graph->HasNode("test 2"));

  // connect test0 to test1 and test2
  graph->AddEdge("test 0", 1, "test 1", 0);
  EXPECT_EQ(graph->NumEdges(), 1);
  graph->AddEdge("test 0", 0, "test 2", 2);
  EXPECT_EQ(graph->NumEdges(), 2);

  graph->Visualize("modules/Modeling/test/WorkGraphVisualizations/UnfixedInOut.pdf");
}

TEST(WorkGraphTests, FixedInOutNum) {
  // create test WorkPieces
  auto test0 = std::make_shared<FixedInOutMod>(2, 3);
  auto test1 = std::make_shared<FixedInOutMod>(1, 1);
  auto test2 = std::make_shared<FixedInOutMod>(2, 0);
  auto test3 = std::make_shared<FixedInOutMod>(1, 2);

  // create and empty graph
  auto graph = std::make_shared<WorkGraph>();

  // make sure the graph is empty
  EXPECT_EQ(graph->NumNodes(), 0);
  EXPECT_EQ(graph->NumEdges(), 0);

  // add WorkPieces to the graph
  EXPECT_FALSE(graph->HasNode("test 0"));
  graph->AddNode(test0, "test 0");
  EXPECT_EQ(graph->NumNodes(), 1);
  EXPECT_TRUE(graph->HasNode("test 0"));

  EXPECT_FALSE(graph->HasNode("test 1"));
  graph->AddNode(test1, "test 1");
  EXPECT_EQ(graph->NumNodes(), 2);
  EXPECT_TRUE(graph->HasNode("test 1"));

  EXPECT_FALSE(graph->HasNode("test 2"));
  graph->AddNode(test2, "test 2");
  EXPECT_EQ(graph->NumNodes(), 3);
  EXPECT_TRUE(graph->HasNode("test 2"));

  EXPECT_FALSE(graph->HasNode("test 3"));
  graph->AddNode(test3, "test 3");
  EXPECT_EQ(graph->NumNodes(), 4);
  EXPECT_TRUE(graph->HasNode("test 3"));

  // connect test0 to test1 and test2 and test 3
  graph->AddEdge("test 0", 1, "test 1", 0);
  EXPECT_EQ(graph->NumEdges(), 1);
  graph->AddEdge("test 0", 0, "test 2", 1);
  EXPECT_EQ(graph->NumEdges(), 2);
  graph->AddEdge("test 0", 0, "test 3", 0);
  EXPECT_EQ(graph->NumEdges(), 3);

  graph->Visualize("modules/Modeling/test/WorkGraphVisualizations/FixedInOutNum.pdf");
}

TEST(WorkGraphTests, FixedInOutType) {
  // the input types
  std::vector<std::string> inTypes({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});
  // the output types
  std::vector<std::string> outTypes({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test0 = std::make_shared<FixedInOutMod>(inTypes, outTypes);
  auto test1 = std::make_shared<FixedInOutMod>(inTypes, outTypes);

  // create and empty graph
  auto graph = std::make_shared<WorkGraph>();

  // make sure the graph is empty
  EXPECT_EQ(graph->NumNodes(), 0);
  EXPECT_EQ(graph->NumEdges(), 0);

  // add WorkPieces to the graph
  EXPECT_FALSE(graph->HasNode("test 0"));
  graph->AddNode(test0, "test 0");
  EXPECT_EQ(graph->NumNodes(), 1);
  EXPECT_TRUE(graph->HasNode("test 0"));

  EXPECT_FALSE(graph->HasNode("test 1"));
  graph->AddNode(test1, "test 1");
  EXPECT_EQ(graph->NumNodes(), 2);
  EXPECT_TRUE(graph->HasNode("test 1"));

  // connect test0 to test1
  graph->AddEdge("test 0", 1, "test 1", 1);
  EXPECT_EQ(graph->NumEdges(), 1);
  graph->AddEdge("test 0", 0, "test 1", 0);
  EXPECT_EQ(graph->NumEdges(), 2);

  graph->Print();
  graph->Visualize("modules/Modeling/test/WorkGraphVisualizations/FixedInOutType.pdf");
}

TEST(WorkGraphTests, DependentCut) {
  // create test WorkPieces
  auto test0 = std::make_shared<FixedInOutMod>(2, 3);
  auto test1 = std::make_shared<FixedInOutMod>(1, 1);
  auto test2 = std::make_shared<FixedInOutMod>(2, 0);
  auto test3 = std::make_shared<FixedInOutMod>(1, 2);

  // create and empty graph
  auto graph = std::make_shared<WorkGraph>();

  // make sure the graph is empty
  EXPECT_EQ(graph->NumNodes(), 0);
  EXPECT_EQ(graph->NumEdges(), 0);

  // add WorkPieces to the graph
  EXPECT_FALSE(graph->HasNode("test 0"));
  graph->AddNode(test0, "test 0");
  EXPECT_EQ(graph->NumNodes(), 1);
  EXPECT_TRUE(graph->HasNode("test 0"));

  EXPECT_FALSE(graph->HasNode("test 1"));
  graph->AddNode(test1, "test 1");
  EXPECT_EQ(graph->NumNodes(), 2);
  EXPECT_TRUE(graph->HasNode("test 1"));

  EXPECT_FALSE(graph->HasNode("test 2"));
  graph->AddNode(test2, "test 2");
  EXPECT_EQ(graph->NumNodes(), 3);
  EXPECT_TRUE(graph->HasNode("test 2"));

  EXPECT_FALSE(graph->HasNode("test 3"));
  graph->AddNode(test3, "test 3");
  EXPECT_EQ(graph->NumNodes(), 4);
  EXPECT_TRUE(graph->HasNode("test 3"));

  // connect test0, test1, test2, and test 3
  graph->AddEdge("test 0", 1, "test 1", 0);
  EXPECT_EQ(graph->NumEdges(), 1);
  graph->AddEdge("test 0", 0, "test 2", 1);
  EXPECT_EQ(graph->NumEdges(), 2);
  graph->AddEdge("test 1", 0, "test 3", 0);
  EXPECT_EQ(graph->NumEdges(), 3);

  // cut the node that does not affect test 3 (test 2)
  auto newGraph = graph->DependentCut("test 3");

  // check the number of nodes in each graph
  EXPECT_EQ(graph->NumNodes(), 4);
  EXPECT_EQ(newGraph->NumNodes(), 3);

  // check the number of edges in each graph
  EXPECT_EQ(graph->NumEdges(), 3);
  EXPECT_EQ(newGraph->NumEdges(), 2);

  newGraph->Visualize("modules/Modeling/test/WorkGraphVisualizations/DependentCut.pdf");
}

TEST(WorkGraphTests, IsConstant) {
  // the input types
  std::vector<std::string> inTypes({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});
  // the output types
  std::vector<std::string> outTypes({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test0 = std::make_shared<FixedInOutMod>(inTypes, outTypes);
  auto test1 = std::make_shared<FixedInOutMod>(2, 2);
  auto test2 = std::make_shared<FixedInOutMod>(3, 2);
  auto test4 = std::make_shared<FixedInOutMod>(1, 1);

    // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(2.0);
  obj->flag = true;

  // make a constant parameter
  auto test5 = std::make_shared<ConstantPiece>(obj, 1);
  auto test3 = std::make_shared<ConstantPiece>((std::string)"string", 3.0);

  // create and empty graph
  auto graph = std::make_shared<WorkGraph>();

  // make sure the graph is empty
  EXPECT_EQ(graph->NumNodes(), 0);
  EXPECT_EQ(graph->NumEdges(), 0);

  // add WorkPieces to the graph
  EXPECT_FALSE(graph->HasNode("test 0"));
  graph->AddNode(test0, "test 0");
  EXPECT_EQ(graph->NumNodes(), 1);
  EXPECT_TRUE(graph->HasNode("test 0"));

  EXPECT_FALSE(graph->HasNode("test 1"));
  graph->AddNode(test1, "test 1");
  EXPECT_EQ(graph->NumNodes(), 2);
  EXPECT_TRUE(graph->HasNode("test 1"));

  EXPECT_FALSE(graph->HasNode("test 2"));
  graph->AddNode(test2, "test 2");
  EXPECT_EQ(graph->NumNodes(), 3);
  EXPECT_TRUE(graph->HasNode("test 2"));

  EXPECT_FALSE(graph->HasNode("test 3"));
  graph->AddNode(test3, "test 3");
  EXPECT_EQ(graph->NumNodes(), 4);
  EXPECT_TRUE(graph->HasNode("test 3"));

  EXPECT_FALSE(graph->HasNode("test 4"));
  graph->AddNode(test4, "test 4");
  EXPECT_EQ(graph->NumNodes(), 5);
  EXPECT_TRUE(graph->HasNode("test 4"));

  EXPECT_FALSE(graph->HasNode("test 5"));
  graph->AddNode(test5, "test 5");
  EXPECT_EQ(graph->NumNodes(), 6);
  EXPECT_TRUE(graph->HasNode("test 5"));

  // connect test0 to test1
  graph->AddEdge("test 4", 0, "test 1", 0);
  EXPECT_EQ(graph->NumEdges(), 1);
  graph->AddEdge("test 0", 1, "test 4", 0);
  EXPECT_EQ(graph->NumEdges(), 2);
  graph->AddEdge("test 2", 0, "test 1", 1);
  EXPECT_EQ(graph->NumEdges(), 3);
  graph->AddEdge("test 3", 0, "test 2", 0);
  EXPECT_EQ(graph->NumEdges(), 4);
  graph->AddEdge("test 3", 1, "test 2", 1);
  EXPECT_EQ(graph->NumEdges(), 5);
  graph->AddEdge("test 5", 0, "test 2", 2);
  EXPECT_EQ(graph->NumEdges(), 6);

  // make sure the nodes are constant (or not)
  EXPECT_FALSE(graph->Constant("test 0"));
  EXPECT_FALSE(graph->Constant("test 4"));
  EXPECT_FALSE(graph->Constant("test 1"));
  EXPECT_TRUE(graph->Constant("test 2"));
  EXPECT_TRUE(graph->Constant("test 3"));
  EXPECT_TRUE(graph->Constant("test 5"));

  // get the outputs for the ConstantPiece node
  const std::vector<boost::any>& outputs = graph->GetConstantOutputs("test 5");
  EXPECT_EQ(outputs.size(), 2);
  EXPECT_DOUBLE_EQ(boost::any_cast<std::shared_ptr<AnObject> >(outputs[0])->value, 2.0);
  EXPECT_TRUE(boost::any_cast<std::shared_ptr<AnObject> >(outputs[0])->flag);
  EXPECT_EQ(boost::any_cast<int>(outputs[1]), 1);

  const std::vector<boost::any>& outputs2 = graph->GetConstantOutputs("test 3");
  EXPECT_EQ(outputs2.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs2[0]).compare((std::string)"string")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs2[1]), 3.0);

  // get the constant parameters of a down stream node
  const std::vector<boost::any>& outputs3 = graph->GetConstantOutputs("test 2");

  EXPECT_EQ(outputs3.size(), 2);
  EXPECT_TRUE(boost::any_cast<std::string>(outputs3[0]).compare((std::string)"string")==0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(outputs3[1]), 3.0);
}

TEST(WorkGraphTests, ConstantDependentCut) {
  // the input types
  std::vector<std::string> inTypes({typeid(std::string).name(), typeid(double).name(), typeid(std::shared_ptr<AnObject>).name()});
  // the output types
  std::vector<std::string> outTypes({typeid(std::string).name(), typeid(double).name()});

    // create the test WorkPiece
  auto test0 = std::make_shared<FixedInOutMod>(inTypes, outTypes);
  auto test1 = std::make_shared<FixedInOutMod>(inTypes, outTypes);

    // create a version of "AnObject"
  auto obj = std::make_shared<AnObject>(2.0);
  obj->flag = true;

  // make a constant parameter
  auto test2 = std::make_shared<ConstantPiece>(obj, 1);
  auto test3 = std::make_shared<ConstantPiece>((std::string)"string", 3.0);

  // create and empty graph
  auto graph = std::make_shared<WorkGraph>();

  // add WorkPieces to the graph
  graph->AddNode(test0, "test 0");
  graph->AddNode(test1, "test 1");
  graph->AddNode(test2, "test 2");
  graph->AddNode(test3, "test 3");
  EXPECT_EQ(graph->NumNodes(), 4);

  // connect test0 to test1
  graph->AddEdge("test 1", 0, "test 0", 0);
  graph->AddEdge("test 1", 1, "test 0", 1);
  graph->AddEdge("test 2", 0, "test 1", 2);
  graph->AddEdge("test 3", 0, "test 1", 0);
  graph->AddEdge("test 3", 1, "test 1", 1);
  EXPECT_EQ(graph->NumEdges(), 5);

  auto newGraph0 = graph->DependentCut("test 1");

  EXPECT_EQ(newGraph0->NumNodes(), 1);
  EXPECT_EQ(newGraph0->NumEdges(), 0);

  auto newGraph1 = graph->DependentCut("test 0");

  EXPECT_EQ(newGraph1->NumNodes(), 2);
  EXPECT_EQ(newGraph1->NumEdges(), 2);

  graph->Visualize("modules/Modeling/test/WorkGraphVisualizations/ConstantDependentCut_BeforeCut.pdf");
  newGraph1->Visualize("modules/Modeling/test/WorkGraphVisualizations/ConstantDependentCut_AfterCut.pdf");
}
