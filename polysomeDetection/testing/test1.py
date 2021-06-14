from polysomeDetection.polysomeDetectionStructures import PolysomeDetectionMRF
from polysomeDetection.sim.simFunctions import generatePositiveList, buildGraph, paintGraph

import time

pointList = generatePositiveList()
print pointList

it = 1
# print 'it', it

g = buildGraph(pointList, 100, 0.9)

#paintGraph(pointList, g)

gm = PolysomeDetectionMRF(g)

print 'init', time.ctime()

mapVec = gm.getMAP(it)

print 'end', time.ctime()

print mapVec

