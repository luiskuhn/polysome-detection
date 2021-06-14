from polysomeDetection.polysomeDetectionStructuresPara import *
from polysomeDetection.graphFunctions import *
from polysomeDetection.postProcFunctions import *
from polysomeDetection.structuralPriors.cytoBacterialPrior import CytoBacterialPrior
from polysomeDetection.eval.evalFunctions import *

import time

##################

from pytom.visualization.particleSpace import *
from pytom.basic.structures import ParticleList

from pytom_volume import read
from pytom.visualization.particleSpace import *

##################


pl_path = '../../../../polysomeDetection/dataset/cytosolic_bacterial/testing/8/pl_labeled_ccFiltered_oriCorr.xml'

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

	# for i in range(len(pl)):
	# 	pl[i].setClass(mapVec[i])

	# pl.toXMLFile('polysomeDetected.xml')

	# true_pl = ParticleList('.')
	# true_pl.fromXMLFile(pl_path)

	# detected_pl = ParticleList('.')
	# detected_pl.fromXMLFile('polysomeDetected.xml')

	# trueSets, detectedSets, noiseSet = getSetsFromParticleLists(true_pl, detected_pl, 5)
	# truePositive, trueNegative, falsePositive, falseNegative = confusionValues(trueSets, detectedSets, noiseSet)

	# print 'A', getMeasuresA(truePositive, trueNegative, falsePositive, falseNegative)
	# print 'B', getMeasuresB(truePositive, trueNegative, falsePositive, falseNegative)

	#################


	templatePath = '70S_bin3.em'
	template = read(templatePath)*-1

	pm = makePolysomeMap(mapVec)

	for k in pm.keys():
		print k, pm[k]

		if len(pm[k]) < 7:
			continue

		poly_pl = []
		for i in pm[k]:
			poly_pl.append(pl[i])

			vol = particleVolume(poly_pl, template, 512, 512, 150)
			vol.write('poly_' + str(k) + '.em')


