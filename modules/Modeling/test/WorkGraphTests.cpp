#include "gtest/gtest.h"

#include "MUQ/Modeling/Core/WorkGraph.h"

#include "WorkPieceTestClasses.h"

using namespace muq::Modeling::Core;

TEST(WorkGraphTests, ConstructGraph) {
  // create a test WorkPiece
  auto test0 = std::make_shared<UnfixedMod>();

  // create and empty graph
  auto graph = std::make_shared<WorkGraph>();

  // make sure the graph is empty
  EXPECT_EQ(graph->NumNodes(), 0);
  EXPECT_EQ(graph->NumEdges(), 0);

  // add a WorkPiece to the graph
  EXPECT_FALSE(graph->HasNode("test0"));
  graph->AddNode(test0, "test0");
  EXPECT_EQ(graph->NumNodes(), 1);
  EXPECT_TRUE(graph->HasNode("test0"));
}
