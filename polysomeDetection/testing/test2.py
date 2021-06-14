from polysomeDetection.polysomeDetectionStructures import PolysomeDetectionMRF
from polysomeDetection.sim.simFunctions import generatePositiveList, buildGraphRand, paintGraph, buildGraph
from polysomeDetection.postProcFunctions import *
import random

pointList = generatePositiveList()
print pointList

it = 1
print 'it', it

#g = buildGraphRand(pointList, 9, 1.0)
g = buildGraph(pointList, 10, 0.9)

#paintGraph(pointList, g)

gm = PolysomeDetectionMRF(g)

mapVec = gm.getMAP(it)

print mapVec

pm = makePolysomeMap(mapVec)

for k in pm.keys():
	random.shuffle(pm[k])
	print k, pm[k], 'ordered:', orderPolysomeSeq(k, pm[k], g)




