from polysomeDetection.polysomeDetectionStructuresPara import *
from polysomeDetection.graphFunctions import *
from polysomeDetection.postProcFunctions import *
from polysomeDetection.structuralPriors.erEukarioticPrior import ErEukarioticPrior

import time

##################

from pytom.visualization.particleSpace import *
from pytom.basic.structures import ParticleList

from pytom_volume import read
from pytom.visualization.particleSpace import *

##################


pl_path = '../../../../polysomeDetection/tool/priorDataset/MOTL_15.0_mirr_scaled.xml'

pl = ParticleList('.')
pl.fromXMLFile(pl_path)

prior = ErEukarioticPrior(1.0)

g = buildGraph(pl, 500, 40, prior)

comm.Barrier()

it = 1
gm = ParaPolysomeDetectionMRF(g)

if rank == 0:
	print pl_path
	print 'init', time.ctime()

mapVec = gm.getMAP(it)

if rank == 0:

	print 'end', time.ctime()

	templatePath = 'template_emd_1093_18.8A_4mumCTF.em'
	template = read(templatePath)#*-1

	pm = makePolysomeMap(mapVec)

	for k in pm.keys():
		print k, pm[k]

		if len(pm[k]) < 3:
			continue

		poly_pl = []
		for i in pm[k]:
			poly_pl.append(pl[i])

			vol = particleVolume(poly_pl, template, 550, 550, 200)
			vol.write('poly_' + str(k) + '.em')


