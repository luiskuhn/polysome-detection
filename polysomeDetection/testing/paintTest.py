from polysomeDetection.graphFunctions import *
from polysomeDetection.structuralPriors.erEukarioticPrior import ErEukarioticPrior
from polysomeDetection.structuralPriors.cytoEukarioticPrior import CytoEukarioticPrior
from pytom.basic.structures import ParticleList

from polysomeDetection.basicStructures import PolysomeParticleList
from polysomeDetection.visualizationFunctions import *

##################

#pl_path = '../../../datasets/inVitro/1_aligned_full_bin1_mirr.xml'
#pl_path = '../../../datasets/inVitro/1_motl_80S_bin2_complete_r9_mirr.xml'
#pl_path = '../../../datasets/inVitro/MOTL_1.0_mirr_binned.xml'

pl_path = '../../../datasets/inVitro/poly_cyto_1_aligned_full_bin1_mirr.xml'

print pl_path

# pl = ParticleList()
# pl.fromXMLFile(pl_path)

pl = PolysomeParticleList()
pl.fromXMLFile(pl_path)

print 'particles', len(pl)

ps = 0.288*2
print 'ps', ps

rnn = 70.0/ps
print 'rnn', rnn

knn = 100

#####

prior = ErEukarioticPrior(ps)
prior = CytoEukarioticPrior(ps)
g = buildGraph(pl, knn, rnn, prior, useShifts=True)

#####

polysomeId = 358
binFactor = 1.0/4


graphToMarkerFile('graph.cmm', pl, g, binFactor=binFactor, xOffset=256, yOffset=256, zOffset=48, useShifts=True, name='graph')
peptideExitsToBildFile('exits.bild', pl, polysomeId, binFactor=binFactor, xOffset=256, yOffset=256, zOffset=48, useShifts=True, riboType='80S')
mRNA2MarkerFile('mRNA.cmm', pl, polysomeId, binFactor=binFactor, xOffset=256, yOffset=256, zOffset=48, useShifts=True)

