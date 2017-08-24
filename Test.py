import pymuqModeling as mm
import numpy as np

id = mm.IdentityPiece(1)
c = mm.ConstantPiece([1.0,2.0])

graph = mm.WorkGraph()
graph.AddNode(c,'x')
graph.AddNode(id, 'y')
graph.AddEdge('x', 0, 'y', 0)

print('Graph size = (%d,%d)'%(graph.NumNodes(), graph.NumEdges()))

fullPiece = graph.CreateWorkPiece('y')
