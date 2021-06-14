import networkx as nx
import numpy as np

from pytom.tools.ProgressBar import FixedProgBar

###########################
##parallel implementation
###########################

import time

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


class ParaMRF:    

    def __init__(self, G):

        self._G = G
        self._totalNodes = self._G.order()
        
        self._vertex2NodeMap = {}
        self.buildVertex2Node()

        if rank > 0:            
            self._vertexList = []
            self.buildVertexList()

            self._beliefs = np.ndarray((len(self._vertexList), self._totalNodes))

            self._msgInMap = {}
            self._msgOutMap = {}
            for n in self._vertexList:
                self._msgInMap[n] = np.ones((len(self._G.in_edges(n)), self._totalNodes))#/totalNodes
                self._msgOutMap[n] = np.ones((len(self._G.out_edges(n)), self._totalNodes))#/totalNodes


    def buildVertex2Node(self):
        
        for vertex in range(self._totalNodes):
            nodeId = ( vertex%( comm.Get_size() -1) ) +1

            self._vertex2NodeMap[vertex] = nodeId

    def buildVertexList(self):

        for vertex in range(self._totalNodes):
            nodeId = self._vertex2NodeMap[vertex]
            if nodeId == rank:
                self._vertexList.append(vertex)

    def compMsgMAP(self, source, target):

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

        self._msgOutMap[source][msgRow, :] = msg
    
    def sendMsgMAP(self, source, target):

        msgRow = self._G.out_edges(source).index((source, target))

        comm.Send(self._msgOutMap[source][msgRow, :], dest = self._vertex2NodeMap[target], tag = source +1)

    def recvMsgMAP(self, source, target):

        msgRow = self._G.in_edges(target).index((source, target))

        msg = np.empty(self._totalNodes)
        comm.Recv(msg, source = self._vertex2NodeMap[source], tag = source +1)

        self._msgInMap[target][msgRow, :] = msg


    def lbpMAP(self, maxIterations):

        barComp = 0
        compCnt = 1
        if rank > 0:
            barComp = FixedProgBar(0, int(len(self._vertexList)*maxIterations), 'RANK ' + str(rank) + ': LBP MAP COMP')


        it = 0        
        while it < maxIterations:

            if rank > 0:

                #compute msgs
                for s in self._vertexList:

                    out_edges = self._G.out_edges(s)

                    for t_index in range(len(out_edges)):

                        t = out_edges[t_index][1]
                        
                        self.compMsgMAP(s, t)

                    barComp.update(compCnt)
                    compCnt = compCnt +1
                
            comm.Barrier()

            if rank == 0:

                barDist = FixedProgBar(0, self._totalNodes, 'LBP MAP DIST')
                distCnt = 0

                #dist msgs
                for s in range(self._totalNodes):
                
                    out_edges = self._G.out_edges(s)
                    
                    for t_index in range(len(out_edges)):
                        
                        t = out_edges[t_index][1]
                        
                        comm.send((0, s, t), dest = self._vertex2NodeMap[s], tag = 1)
                        comm.send((1, s, t), dest = self._vertex2NodeMap[t], tag = 1)
                        
                        done = comm.recv(source = self._vertex2NodeMap[t], tag = 1)

                    barDist.update(distCnt)
                    distCnt = distCnt +1

                for n in range(1, comm.Get_size()):
                    comm.send((-1, 0, 0), dest = n, tag = 1)

            else:

                while True:
                    
                    cmd = comm.recv(source = 0, tag = 1)
                    
                    op = cmd[0]
                    source = cmd[1]
                    target = cmd[2]

                    if op == 0:

                        self.sendMsgMAP(source, target)
                        
                    elif op == 1:

                        self.recvMsgMAP(source, target)                      
                        comm.send(0, dest = 0, tag = 1)

                    elif op == -1:
                        break

            comm.Barrier()

            it = it +1

    def computeBeliefsMAP(self):

        for i in range(len(self._vertexList)):
            for x_i in self._G.nodes():
                sum_x_i = np.sum(self._msgInMap[self._vertexList[i]][:, x_i])
                self._beliefs[i, x_i] = self.phiPotential(self._vertexList[i], x_i) + sum_x_i

            # #norm
            # self._beliefs[i, :] = self._beliefs[i, :]/self._beliefs[i, :].sum()

    def sendMinMAP(self):

        for i in range(len(self._vertexList)):
            minLabel = self._beliefs[i, :].argmin()
            comm.send(minLabel, dest = 0, tag = self._vertexList[i] +1)
        
    def getMAP(self, maxIterations):

        returnValue = 0

        comm.Barrier()
        self.lbpMAP(maxIterations)
        comm.Barrier()
        if rank > 0:
            self.computeBeliefsMAP()
        comm.Barrier()
        if rank > 0:
            self.sendMinMAP()
        comm.Barrier()

        if rank == 0:
            mapVector = np.zeros(self._totalNodes)
            for i in self._G.nodes():
                mapVector[i] = comm.recv(source = self._vertex2NodeMap[i], tag= i +1)

            returnValue = mapVector

        return returnValue
 

