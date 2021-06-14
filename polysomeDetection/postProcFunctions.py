
def makePolysomeMap(mapVector):

    n = mapVector.shape[0]
    fm = {}
    for i in range(n):        
        x_i = int(mapVector[i])

        if fm.has_key(x_i):
            fm[x_i].append(i)
        else:
            fm[x_i] = [i]

    return fm

def makeLenMap(polysomeMap):

    lenMap = {}
    for k in polysomeMap.keys():

        polyLen = len(polysomeMap[k])

        if lenMap.has_key(polyLen):
            lenMap[polyLen].append(k)
        else:
            lenMap[polyLen] = [k]

    return lenMap

def transitivityCorrect(polysomeMap, mapVector):

    for k in polysomeMap.keys():
        if not k in polysomeMap[k]:
            
            transKey = 0
            for i in polysomeMap.keys():
                if k in polysomeMap[i]:
                    transKey = i
                    break
            
            polysomeMap[transKey].extend(polysomeMap[k])
            
            for j in polysomeMap[k]:
                mapVector[j] = transKey

            del polysomeMap[k]

def greedyPolysomeWalk(G, n):
    import numpy as np
        
    outEdges = G.out_edges(n, True)

    if len(outEdges) == 0:
        return [[n]]
    else:
        newWalks = []

        weight_vec = np.zeros(len(outEdges))
        for i in range(len(outEdges)):
            weight_vec[i] = outEdges[i][2]['neg_log']
        sortList = weight_vec.argsort().tolist()

        for i in range(len(outEdges)):
            edge = outEdges[sortList[i]]
            walks = greedyPolysomeWalk(G, edge[1])

            for w in range(len(walks)):
                newWalk = [n]
                newWalk.extend(walks[w])
                walks[w] = newWalk

            newWalks.extend(walks)

        return newWalks

def polysomeWalk(G, n):
        
    outEdges = G.out_edges(n)

    if len(outEdges) == 0:
        return [[n]]
    else:
        newWalks = []
        for i in range(len(outEdges)):
            edge = outEdges[i]
            walks = polysomeWalk(G, edge[1])

            for w in range(len(walks)):
                newWalk = [n]
                newWalk.extend(walks[w])
                walks[w] = newWalk

            newWalks.extend(walks)

        return newWalks

def inclusiveOrder(end, polysomeSeq, G):
    import networkx as nx
    import numpy as np

    poly_end = end    
    subG = G.subgraph(polysomeSeq)

    #print 'reducing sub-graph...'
    reduceGraph(subG)

    paths = nx.single_source_dijkstra_path(subG, poly_end, weight='neg_log')
    posList = []
    for vertex in polysomeSeq:
        if paths.has_key(vertex):
            posList.append( len(paths[vertex]) -1 )
        else:
            posList.append(1)


    return posList

def reduceGraph(g):
    import networkx as nx
    import numpy as np

    for n in g.nodes():
        in_e = g.in_edges(n, data=True)
        if len(in_e) > 0:
            neg_log_list = []
            for e in in_e:
                neg_log = e[2]['neg_log']
                neg_log_list.append(neg_log)
            neg_log_list = np.array(neg_log_list)
            min_index = neg_log_list.argmin()
            for e_i in range(len(in_e)):
                e = in_e[e_i]
                if e_i == min_index:
                    e[2]['neg_log'] = -np.log(0.999)
                else:
                    #g.remove_edge(e[0], e[1])
                    e[2]['neg_log'] = -np.log(1e-30)







def getMaxBranchLen(end, polysomeSeq, G):
    import networkx as nx
    import numpy as np

    poly_end = end
    subG = G.subgraph(polysomeSeq)

    ######

    endVec = np.zeros(len(polysomeSeq))
    for i in range(len(polysomeSeq)):

        n = polysomeSeq[i]            
        inEdges = subG.in_edges(n, True)
        
        if len(inEdges) == 0:
            continue

        logSum = 0.0
        for edge in inEdges:
            w = edge[2]['weight']
            logSum = logSum + (-np.log(1.0 - w))

        endVec[i] = logSum

    poly_end = polysomeSeq[endVec.argmin()]

    if endVec[polysomeSeq.index(end)] <= endVec[polysomeSeq.index(poly_end)]:
        poly_end = end
    
    ######

    paths = nx.single_source_dijkstra_path(subG, poly_end, weight='neg_log')
    maxLen = 0
    for vertex in polysomeSeq:
        pathLen = len(paths[vertex])
        if pathLen > maxLen:
            maxLen = pathLen

    return maxLen

def makeMaxBranchLenMap(polysomeMap, G):

    lenMap = {}
    for k in polysomeMap.keys():
        maxBranchLen = getMaxBranchLen(k, polysomeMap[k], G)

        if lenMap.has_key(maxBranchLen):
            lenMap[maxBranchLen].append(k)
        else:
            lenMap[maxBranchLen] = [k]
    return lenMap



# def orderLongestBranch(end, polysomeSeq, G):
#     import numpy as np
#     import networkx as nx

#     poly_end = 0
#     subG = G.subgraph(polysomeSeq)

#     if subG.in_degree(end) == 0:
#         poly_end = end
#     else:
#         endProbVec = np.ones(len(polysomeSeq))
#         for i in range(len(polysomeSeq)):

#             n = polysomeSeq[i]            
#             inEdges = subG.in_edges(n, True)
            
#             if len(inEdges) == 0:
#                 continue

#             endProb = 1.0
#             for edge in inEdges:
#                 w = edge[2]['weight']
#                 endProb = endProb * (1.0 - w)

#             endProbVec[i] = endProb

#         poly_end = polysomeSeq[endProbVec.argmax()]

#     ######

#     w_cutoff = 0.0
#     w_cutoff_step = 0.00000001

#     while nx.is_directed_acyclic_graph(subG) == False and len(subG.edges()) > 0:

#         edgeList = subG.edges(data=True)

#         for edge in edgeList:

#             w = edge[2]['weight']

#             if w <= w_cutoff:
#                 subG.remove_edge(edge[0], edge[1])

#                 if nx.is_weakly_connected(subG) == False:
#                     subG.add_edge(edge[0], edge[1], weight = w)

#         w_cutoff = w_cutoff + w_cutoff_step

#     ######    

#     lists = polysomeWalk(subG, poly_end)

#     #print lists

#     maxList = []
#     for l in lists:
#         if len(l) > len(maxList):
#             maxList = l
    
#     #print 'max', maxList

#     return maxList




