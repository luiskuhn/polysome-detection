from mrf.mrfStructures import *

# class PolysomeDetectionMRF(MRF):
    
#     def __init__(self, G):

#         MRF.__init__(self, G)    
        
#         ####should go elsewhere
#         self._pathProb = np.ndarray((self._totalNodes, self._totalNodes))
#         self.buildPaths()
#         self._rho = 0.0001
#         #############

#     ######should go elsewhere

#     def buildPaths(self):

#         self._paths = nx.all_pairs_dijkstra_path(self._G, weight='neg_log')
#         for s in self._G.nodes():
#             for t in self._G.nodes():

#                 Pp = 0

#                 if s == t:
#                     inEdges = self._G.in_edges(s, True)

#                     # endProb = 0.0
#                     # for edge in inEdges:
#                     #     w = edge[2]['weight']
#                     #     endProb = endProb + (1.0 - w)
#                     # Pp = endProb

#                     endProb = 1.0
#                     for edge in inEdges:
#                         w = edge[2]['weight']
#                         endProb = endProb * (1.0 - w)
#                     Pp = endProb     
                                
#                 elif self._paths[s].has_key(t):
#                     path = self._paths[s][t]
#                     #print s, t, path                   
                    
#                     # min = 1             
#                     # for n in range(len(path) -1):                       
#                     #     w = self._G[ path[n] ][ path[n +1] ]['weight']
#                     #     if w < min:
#                     #         min = w

#                     # Pp = min

#                     pProb = 0.0             
#                     for n in range(len(path) -1):                       
#                         w = self._G[ path[n] ][ path[n +1] ]['neg_log']
#                         pProb = pProb + w                       

#                     Pp = np.e**(-pProb)
                    
#                     #print 'pathScore', Pp

#                 self._pathProb[s, t] = Pp

#     ########################

    
#     #########should go elsewhere
#     def psiPotential(self, source, target, fs, ft):

#         potential = 0       
#         jumpProb = self._G.get_edge_data(source, target)['weight']

#         ###############

#         ftPath = []
#         if self._paths[ft].has_key(target):
#             ftPath = self._paths[ft][target]
                
#         ###############

#         inPath = False
#         if source in ftPath:
#             sourceIndex = ftPath.index(source)
#             inPath = len(ftPath) >= (sourceIndex +1) +1 and ftPath[sourceIndex +1] == target 

#         ###############

#         if fs == ft and inPath:
#             potential = jumpProb
#         else:
#             potential = self._rho**2##0
        
#         return 1.0 - potential
        
#     def phiPotential(self, source, fs):
        
       
#         w0 = 0.5 #0.5, 0.65
#         w1 = 1.0 -w0

#         pathScore = self._pathProb[fs, source]
#         fatherScore = self._pathProb[fs, fs]

#         potential = 0

#         if fs == source:
#             potential = (w0*fatherScore) + w1
#             #potential = fatherScore
#         elif pathScore > 0:
#             comp0 = w0*fatherScore
#             comp1 = w1*pathScore
#             potential =  comp0 + comp1
#             #potential = pathScore*fatherScore
#         else:
#             potential = 0#self._rho**2##0
               

#         return 1.0 - potential            

#     #######################


class PolysomeDetectionMRF(MRF):
    
    def __init__(self, G):

        MRF.__init__(self, G)    
        
        self._pathProb = np.ndarray((self._totalNodes, self._totalNodes))
        self.buildPaths()


    def buildPaths(self):

        self._paths = nx.all_pairs_dijkstra_path(self._G, weight='neg_log')
        for s in self._G.nodes():
            for t in self._G.nodes():

                Pp = 1e-20

                if s == t:
                    inEdges = self._G.in_edges(s, True)

                    endProb = 1.0
                    for edge in inEdges:
                        w = edge[2]['weight']
                        endProb = endProb * (1.0 - w)
                    Pp = np.max((endProb, 1e-20))     
                                
                elif self._paths[s].has_key(t):
                    path = self._paths[s][t]
                    
                    pProb = 0.0             
                    for n in range(len(path) -1):                       
                        w = self._G[ path[n] ][ path[n +1] ]['neg_log']
                        pProb = pProb + w                       

                    Pp = np.e**(-pProb)

                self._pathProb[s, t] = Pp


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




