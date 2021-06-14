import networkx as nx
import numpy as np

from pytom.tools.ProgressBar import FixedProgBar

##sequential implementation
###########################

class MRF:

    def __init__(self, G):

        self._G = G

        self._totalNodes = self._G.order()

        self._beliefs = np.ndarray((self._totalNodes, self._totalNodes))

        self._msgMap = {}
        for n in self._G.nodes():
            self._msgMap[n] = np.ones((len(self._G.in_edges(n)), self._totalNodes))#/totalNodes

    def msgMAP(self, source, target):

        msgRow = self._G.out_edges(source).index((source, target))
        target_msgRow = -1
        if (target, source) in self._G.in_edges(source):
            target_msgRow = self._G.in_edges(source).index((target, source))

        phi_vec = np.zeros(self._totalNodes)
        sum_n_vec = np.zeros(self._totalNodes)
        for x_j in range(self._totalNodes):
            phi_vec[x_j] = self.phiPotential(source, x_j)
            red_val = 0.0
            if target_msgRow >= 0:
                red_val = self._msgInMap[source][target_msgRow, x_j] #N(j)-{i}
            sum_n_vec[x_j] = np.sum(self._msgInMap[source][:, x_j]) - red_val

        msg = np.zeros(self._totalNodes)
        for x_i in range(self._totalNodes):
            x_j_vec = np.zeros(self._totalNodes)
            for x_j in range(self._totalNodes):
                psi = self.psiPotential(source, target, x_j, x_i)
                x_j_vec[x_j] = phi_vec[x_j] + sum_n_vec[x_j] + psi
            msg[x_i] = x_j_vec.min()

        # #norm
        # msg = msg/msg.sum()

        self._msgMap[target][msgRow, :] = msg



    def lbpMAP(self, maxIterations):

        bar = FixedProgBar(0, int(maxIterations*self._totalNodes), 'LBP MAP')
        barCnt = 1
    
        it = 0
        
        while it < maxIterations:
            
            for s in self._G.nodes():
                
                out_edges = self._G.out_edges(s)
                
                for t_index in range(len(out_edges)):
                    
                    t = out_edges[t_index][1]
                    
                    self.msgMAP(s, t)

                bar.update(barCnt)
                barCnt = barCnt +1
                                
            it = it +1

    def computeBeliefsMAP(self):

        for i in self._G.nodes():
            for x_i in self._G.nodes():
                sum_x_i = np.sum(self._msgMap[i][:, x_i])
                self._beliefs[i, x_i] = self.phiPotential(i, x_i) + sum_x_i

            # #norm
            # self._beliefs[i, :] = self._beliefs[i, :]/self._beliefs[i, :].sum()

    def getMAP(self, maxIterations):

        self.lbpMAP(maxIterations)
        self.computeBeliefsMAP()

        mapVector = np.zeros(self._totalNodes)
        for i in self._G.nodes():            
            mapVector[i] = self._beliefs[i, :].argmin()

        return mapVector

    def getBeliefs(self, n):
        
        return self._beliefs[n, :]
