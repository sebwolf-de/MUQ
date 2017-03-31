#include "gtest/gtest.h"

#include "MUQ/Modeling/WorkGraph.h"

#include "WorkPieceTestClasses.h"

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
  std::vector<std::string> inTypes({typeid(std::string).name(), typeid(std::shared_ptr<AnObject>).name(), typeid(double).name()});
  // the output types
  std::vector<std::string> outTypes({typeid(std::string).name(), typeid(double).name()});

  // create the test WorkPiece
  auto test0 = std::make_shared<FixedTypesMod>(inTypes, outTypes);
  auto test1 = std::make_shared<FixedTypesMod>(inTypes, outTypes);
  
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
  graph->AddEdge("test 0", 1, "test 1", 2);
  EXPECT_EQ(graph->NumEdges(), 1);
  graph->AddEdge("test 0", 0, "test 1", 0);
  EXPECT_EQ(graph->NumEdges(), 2);

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
