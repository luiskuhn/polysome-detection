from mrf.mrfStructuresPara import *

###########################
##parallel implementation
###########################

class ParaPolysomeDetectionMRF(ParaMRF):
    
    def __init__(self, G):

        ParaMRF.__init__(self, G)    
        
        self._pathProb = np.ndarray((self._totalNodes, self._totalNodes))
        self.buildPaths()


    def buildPaths(self):

        self._paths = nx.all_pairs_dijkstra_path(self._G, weight='neg_log')
        for s in self._G.nodes():
            for t in self._G.nodes():

                prob = 1e-50

                if s == t:
                    inEdges = self._G.in_edges(s, True)
                    
                    logProb = 0.0
                    for edge in inEdges:
                        w = edge[2]['weight']
                        logProb = logProb + np.log(1.0 - w)
                    endProb = np.exp(logProb)
                    prob = np.max((endProb, prob))     
                                
                elif self._paths[s].has_key(t):
                    path = self._paths[s][t]
                    
                    pathCost = 0.0             
                    for n in range(len(path) -1):                       
                        negLog = self._G[ path[n] ][ path[n +1] ]['neg_log']
                        pathCost = pathCost + negLog
                    pathProb = np.exp(-pathCost)
                    prob = np.max((pathProb, prob))

                self._pathProb[s, t] = prob


    def psiPotential(self, source, target, fs, ft):

        potential = 0       
        jumpProb = self._G.get_edge_data(source, target)['weight']

        ###############

        ftPath = []
        if self._paths[ft].has_key(target):
            ftPath = self._paths[ft][target]
                
        ###############

        inPath = False
        if source in ftPath:
            sourceIndex = ftPath.index(source)
            inPath = len(ftPath) >= (sourceIndex +1) +1 and ftPath[sourceIndex +1] == target 

        ###############

        if fs == ft and inPath:
            potential = -np.log(jumpProb)
        else:
            potential = -np.log(1e-20)
        
        return potential
        
    def phiPotential(self, source, fs):

        pathScore = self._pathProb[fs, source]
        fatherScore = self._pathProb[fs, fs]

        potential = 0

        if fs == source:
            potential = -np.log(fatherScore)

        elif pathScore > 0:
            potential =  -(np.log(fatherScore) + np.log(pathScore))

        else:
            potential = -np.log(1e-20)


        return potential


