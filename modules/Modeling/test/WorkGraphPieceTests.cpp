#include "gtest/gtest.h"

#include "MUQ/Modeling/WorkGraph.h"

#include "WorkPieceTestClasses.h"

using namespace muq::Modeling;

TEST(WorkGraphPiece, FixedInOutNum) {
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

  graph->Visualize("modules/Modeling/test/WorkGraphVisualizations/FixedInOutNum_WorkPiece.pdf");
}
