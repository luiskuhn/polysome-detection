
import networkx as nx
from scipy.spatial import KDTree
import numpy as np

from polysomeDetection.postProcFunctions import *

def weightFnc(i, j, pl, poslen):

    w = 0.1

    a = np.array(pl[i])
    b = np.array(pl[j])

    dist = np.linalg.norm(a-b)

    if dist <= 2 and i < j and i < poslen and j < poslen:
        w = 0.99

    return w

def getInputNeighbors(q, tree):

    knn = 200
    rnn = 3.0

    d, i = tree.query(q, k=knn, eps=0.0, p=2.0, distance_upper_bound=rnn)
    d = list(d)
    i = list(i)

    nn = 0
    if d.count(float('inf')) > 0:
        infIndex = d.index(float('inf'))
        nn = i[:infIndex]
    else:
        nn = i

    return nn

def buildGraph(pl):

    import random

    poslen = len(pl)
    randPoints = 20
    for i in range(randPoints):
        x = random.uniform(0, 15)
        y = random.uniform(0, 15)
        pl.append((x,y))

    pointTree = KDTree(list(pl))
    g = nx.DiGraph()

    g.add_nodes_from(range(len(pl)))

    for i in range(len(pl)):

        neighbors = getInputNeighbors(pl[i], pointTree)

        for j in neighbors:

            #still check for self-entry
            if i == j:
                continue

            w = weightFnc(i, j, pl, poslen)

            g.add_edge(i, j, weight = w)
            g[i][j]['neg_log'] = -np.log(w)


    return g

def paintGraph(pl, g):

    import matplotlib.pyplot as plt

    pointArray = np.array(pl)

    plt.plot(pointArray[:,0], pointArray[:,1], 'ro')
    for e in g.edges():
        plt.plot([pointList[e[0]][0], pointList[e[1]][0]], [pointList[e[0]][1], pointList[e[1]][1]], 'r', alpha=g[e[0]][e[1]]['weight'])

    #plt.savefig("edge_colormap.png") # save as png
    plt.ylim([-3,20])
    plt.xlim([-3,20])

    plt.show() # display

###################################################

from polysomeDetection.polysomeDetectionStructures import PolysomeDetectionMRF

###################################################

pointList = [(0, 0),
             (1, 1),
             (10, 6),
             (11, 7),
             (12, 8),
             (3, 8),
             (4, 9),
             (5, 10),
             (6, 11)]



it = 1
print 'it', it
g = buildGraph(pointList)

paintGraph(pointList, g)

gm = PolysomeDetectionMRF(g)

mapVec = gm.getMAP(it)

# for i in range(len(pointList)):
#     print i, gm.getBeliefs(i)

print mapVec

# pm = makePolysomeMap(mapVec)
# print pm
# transitivityCorrect(pm)
# print 'trans corr:'
# print pm



