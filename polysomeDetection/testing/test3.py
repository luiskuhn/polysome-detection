from polysomeDetection.sim.simFunctions import generatePositiveList, addRandPoints, buildGraphDirect, paintGraph

pointList = generatePositiveList()
poslen = len(pointList)
addRandPoints(pointList, 2000)

g = buildGraphDirect(pointList, 0.9, poslen)

paintGraph(pointList, g)

