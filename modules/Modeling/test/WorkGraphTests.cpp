#include "gtest/gtest.h"

#include "MUQ/Modeling/WorkGraph.h"

#include "WorkPieceTestClasses.h"

using namespace muq::Modeling;

TEST(WorkGraphTests, ConstructGraph) {
  // create a test WorkPiece
  auto test0 = std::make_shared<UnfixedMod>();
  auto test1 = std::make_shared<UnfixedMod>();

  // create and empty graph
  auto graph = std::make_shared<WorkGraph>();

  // make sure the graph is empty
  EXPECT_EQ(graph->NumNodes(), 0);
  EXPECT_EQ(graph->NumEdges(), 0);

  // add WorkPieces to the graph
  EXPECT_FALSE(graph->HasNode("test0"));
  graph->AddNode(test0, "test0");
  EXPECT_EQ(graph->NumNodes(), 1);
  EXPECT_TRUE(graph->HasNode("test0"));

  EXPECT_FALSE(graph->HasNode("test1"));
  graph->AddNode(test1, "test1");
  EXPECT_EQ(graph->NumNodes(), 2);
  EXPECT_TRUE(graph->HasNode("test1"));

  // connect test0 to test1
  graph->AddEdge("test0", 0, "test1", 0);
  EXPECT_EQ(graph->NumEdges(), 1);

  graph->Visualize("graph.pdf");
}
