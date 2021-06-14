import networkx as nx
import numpy as np
from scipy.spatial import KDTree
import random

def weightFnc(i, j, pl, poslen, w_pos):
    w_neg = 0.05
    #w_pos = 0.99
    max_dist = 2.0

    w = w_neg

    a = np.array(pl[i])
    b = np.array(pl[j])

    dist = np.linalg.norm(a-b)

    if dist <= max_dist and i == (j -1) and i < poslen and j < poslen:
        w = w_pos

    return w

def getInputNeighbors(q, tree):
    knn = 500
    rnn = 4.0

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

def addRandPoints(pl, randPoints):
    x_lim = (0, 50)
    y_lim = (0, 50)

    poslen = len(pl)
    
    for i in range(randPoints):
        x = random.uniform(x_lim[0], x_lim[1])
        y = random.uniform(y_lim[0], y_lim[1])
        pl.append((x,y))

def buildGraphDirect(pl, w_pos, poslen):
    pointTree = KDTree(list(pl))
    g = nx.DiGraph()

    g.add_nodes_from(range(len(pl)))

    for i in range(len(pl)):

        neighbors = getInputNeighbors(pl[i], pointTree)

        for j in neighbors:

            #still check for self-entry
            if i == j:
                continue

            w = weightFnc(i, j, pl, poslen, w_pos)

            g.add_edge(i, j, weight = w)
            g[i][j]['neg_log'] = -np.log(w)

    return g

def buildGraph(pl, randPoints, w_pos):
    #randPoints = 200
    x_lim = (0, 20)
    y_lim = (0, 20)

    poslen = len(pl)
    
    for i in range(randPoints):
        x = random.uniform(x_lim[0], x_lim[1])
        y = random.uniform(y_lim[0], y_lim[1])
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

            w = weightFnc(i, j, pl, poslen, w_pos)

            g.add_edge(i, j, weight = w)
            g[i][j]['neg_log'] = -np.log(w)

    return g

def buildGraphRand(pl, randPoints, sd):
    #randPoints = 200
    x_lim = (0, 20)
    y_lim = (0, 20)

    poslen = len(pl)
    
    for i in range(randPoints):
        x = random.uniform(x_lim[0], x_lim[1])
        y = random.uniform(y_lim[0], y_lim[1])
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

            w = weightFnc(i, j, pl, poslen, max(0.05, 1.0 -abs(random.gauss(0, sd))) )

            g.add_edge(i, j, weight = w)
            g[i][j]['neg_log'] = -np.log(w)

    return g


def paintGraph(pl, g):

    import matplotlib.pyplot as plt

    pointArray = np.array(pl)

    #plt.plot(pointArray[:,0], pointArray[:,1], 'ko', markersize=8)
    plt.plot(pointArray[:,0], pointArray[:,1], 'bo')

    for e in g.edges():
        plt.plot([pointArray[e[0]][0], pointArray[e[1]][0]], [pointArray[e[0]][1], pointArray[e[1]][1]], 'b', alpha=g[e[0]][e[1]]['weight'])

    #plt.savefig("edge_colormap.png") # save as png
    plt.ylim([-3, 50])
    plt.xlim([-3, 50])

    plt.show() # display


def generatePositiveList():
    seq_num = 3
    init_size = 2
    rng = 100

    sec = rng/seq_num

    pointList = []

    for i in range(seq_num):
        x = random.randint(i*sec, (i +1)*sec)
        y = random.randint(rng -((i+1)*sec), rng -((i)*sec) )
        
        for j in range(init_size):

            pointList.append((x, y))

            x = x +1
            y = y +1

        init_size = init_size +1

    return pointList

