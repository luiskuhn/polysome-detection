from polysomeDetection.polysomeDetectionStructuresPara import *
from polysomeDetection.graphFunctions import *
from polysomeDetection.postProcFunctions import *
from polysomeDetection.structuralPriors.cytoBacterialPrior import CytoBacterialPrior
from polysomeDetection.basicStructures import PolysomeParticleList

import time

##################

from pytom.visualization.particleSpace import *
from pytom.basic.structures import ParticleList

from pytom_volume import read
from pytom.visualization.particleSpace import *

##################


#pl_path = '../../../../polysomeDetection/dataset/cytosolic_bacterial/testing/8/pl_labeled_ccFiltered_oriCorr.xml'
pl_path = '../../../../polysomeDetection/dataset/cytosolic_bacterial/testing/sim/4/pl_labeled_filtered.xml'

pl = ParticleList('.')
pl.fromXMLFile(pl_path)

prior = CytoBacterialPrior(2.24)

g = buildGraph(pl, 500, 60/2.24, prior)

comm.Barrier()

it = 1
gm = ParaPolysomeDetectionMRF(g)

if rank == 0:
	print pl_path
	print 'init', time.ctime()

mapVec = gm.getMAP(it)

if rank == 0:

	print 'end', time.ctime()

	poly_pl = PolysomeParticleList(mapVec, g, pl=pl)

	print 'pl', len(pl)
	
	for k in poly_pl._polysomeMap.keys():
		if len(poly_pl._polysomeMap[k]) >= 6:
			print poly_pl._polysomeMap[k]	

	poly_pl.toXMLFile('polysomeList.xml')

	print '###########################'

	file_poly_pl = PolysomeParticleList()
	file_poly_pl.fromXMLFile('polysomeList.xml')

	for k in file_poly_pl._polysomeMap.keys():
		if len(file_poly_pl._polysomeMap[k]) >= 6:
			print file_poly_pl._polysomeMap[k]

	print '###########################'

	print file_poly_pl.getPolysomesBySize(6)

	

	