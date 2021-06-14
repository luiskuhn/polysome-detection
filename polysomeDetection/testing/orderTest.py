import numpy as np
import networkx as nx
g = nx.DiGraph()

g.add_nodes_from(range(4))

g.add_edge(0,1, weight = 1)
g.add_edge(1,2, weight = 1)
g.add_edge(2,3, weight = 1)

g.add_edge(1,4, weight = .1)
g.add_edge(4,2, weight = .1)
g.add_edge(2,1, weight = .01)

#g.remove_edge(1,4)

