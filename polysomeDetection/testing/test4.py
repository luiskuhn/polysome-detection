from polysomeDetection.polysomeDetectionStructures import PolysomeDetectionMRF
from polysomeDetection.graphFunctions import *
from polysomeDetection.postProcFunctions import *
from polysomeDetection.structuralPriors.cytoBacterialPrior import CytoBacterialPrior

##################

from pytom.visualization.particleSpace import *
from pytom.basic.structures import ParticleList

from pytom_volume import read
from pytom.visualization.particleSpace import *

# ##################

templatePath = '70S_bin3.em'
template = read(templatePath)*-1

pl_path = '../../../../polysomeDetection/dataset/cytosolic_bacterial/testing/8/pl_labeled_ccFiltered_oriCorr.xml'

print pl_path

pl = ParticleList('.')
pl.fromXMLFile(pl_path)

print 'particles', len(pl)

#####

prior = CytoBacterialPrior(2.24)

#####

g = buildGraph(pl, 500, 30/2.24, prior)

paint2dGraph(pl, g)

# it = 1
# print 'it', it

# gm = PolysomeDetectionMRF(g)

# mapVec = gm.getMAP(it)

# #print mapVec

# pm = makePolysomeMap(mapVec)

# for k in pm.keys():
# 	print k, pm[k]

# 	poly_pl = []
# 	for i in pm[k]:
# 		poly_pl.append(pl[i])

# 		vol = particleVolume(poly_pl, template, 512, 512, 150)
# 		vol.write('poly_' + str(k) + '.em')

