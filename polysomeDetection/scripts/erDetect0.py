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

templatePath = 'un_mirr.em'

pl = ParticleList()
#pl.fromMOTL('/fs/pool/pool-titomo1/InVitroTranslocRetic/MEF/tomo02/template_matching/motlfiles/motl_80S_bin2_mirr_complete_r9.em', '', '')
pl.fromXMLFile('motl_80S_bin2_mirr_complete_r9_mirr.xml')

ps = 2.1
knn= 500
rnn = 100.0 #nm

prior = ErEukarioticPrior(ps)

g = buildGraph(pl, knn, rnn/ps, prior)

comm.Barrier()

it = 1
gm = ParaPolysomeDetectionMRF(g)

if rank == 0:
	#print pl_path
	print 'init', time.ctime()

mapVec = gm.getMAP(it)

if rank == 0:

	print 'end', time.ctime()

	template = read(templatePath)

	pm = makePolysomeMap(mapVec)

	for k in pm.keys():
		print k, pm[k]

		if len(pm[k]) < 3:
			continue

		poly_pl = []
		for i in pm[k]:
			poly_pl.append(pl[i])

			vol = particleVolume(poly_pl, template, 463, 463, 128)
			vol.write('poly_' + str(k) + '.em')


