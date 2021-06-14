from polysomeDetection.polysomeDetectionStructuresPara import *
from polysomeDetection.graphFunctions import *
from polysomeDetection.structuralPriors.erEukarioticPrior import ErEukarioticPrior
from polysomeDetection.basicStructures import PolysomeParticleList

import time
import sys

##################

from pytom.basic.structures import ParticleList

##################

pl_in = sys.argv[1]
pl_out = sys.argv[2]

ps = float(sys.argv[3]) #nm
knn= int(sys.argv[4])
rnn = float(sys.argv[5]) #nm

pl = ParticleList()
pl.fromXMLFile(pl_in)

prior = ErEukarioticPrior(ps)

g = buildGraph(pl, knn, rnn/ps, prior, useShifts=True)

comm.Barrier()

it = 1 #3
gm = ParaPolysomeDetectionMRF(g)

comm.Barrier()

if rank == 0:
	#print pl_path
	print 'init', time.ctime()

mapVec = gm.getMAP(it)

if rank == 0:

	print 'end', time.ctime()

	#poly_pl = PolysomeParticleList(mapVec, g, pl = pl, correctLoops = True, order = 'longestBranchDFS')
	poly_pl = PolysomeParticleList(mapVec, g, pl = pl, correctLoops = True)
	#poly_pl = PolysomeParticleList(mapVec, g, pl = pl)
				
	for k in poly_pl._polysomeMap.keys():
		print k, poly_pl._polysomeMap[k]	

	poly_pl.toXMLFile(pl_out)

	

	


