import networkx as nx
import numpy as np
from scipy.spatial import KDTree

def getInputNeighbors(q, tree, knn, rnn):

    d, i = tree.query(q, k=knn, eps=0.0, p=2.0, distance_upper_bound=rnn)
    d = list(d)
    i = list(i)

    nn = 0
    if d.count(float('inf')) > 0:
        infIndex = d.index(float('inf'))
        nn = i[:infIndex]
    else:
        nn = i

    for index in range(len(nn)):
        nn[index] = int(nn[index])

    return nn

def buildGraph(pl, knn, rnn, prior, useShifts=False):

    pointList = []

    for p in pl:

        x = p.getPickPosition().getX()
        y = p.getPickPosition().getY()
        z = p.getPickPosition().getZ()
        if useShifts:
            x = x + p.getShift().getX()
            y = y + p.getShift().getY()
            z = z + p.getShift().getZ()

        pointList.append((x, y, z))
    
    pointTree = KDTree(list(pointList))
    
    g = nx.DiGraph()

    g.add_nodes_from(range(len(pl)))
    
    for i in range(len(pl)):

        x = pl[i].getPickPosition().getX()
        y = pl[i].getPickPosition().getY()
        z = pl[i].getPickPosition().getZ()
        if useShifts:
            x = x + pl[i].getShift().getX()
            y = y + pl[i].getShift().getY()
            z = z + pl[i].getShift().getZ()
            
        q = (x, y, z)

        neighbors = getInputNeighbors(q, pointTree, knn, rnn)

        for j in neighbors:

            #still check for self-entry
            if i == j:
                continue

            w = prior.weightEdge(i, j, pl, useShifts=useShifts)

            g.add_edge(i, j, weight = w)
            g[i][j]['neg_log'] = -np.log(w)
    
    return g

