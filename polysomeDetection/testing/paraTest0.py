from polysomeDetection.polysomeDetectionStructuresPara import *
from polysomeDetection.sim.simFunctions import generatePositiveList, addRandPoints, buildGraphDirect

import time

#print 'size:', comm.Get_size()

# if rank == 0:
# 	pointList = generatePositiveList()
# 	poslen = len(pointList)
# 	addRandPoints(pointList, 100)

# 	comm.send(pointList, dest=1, tag=1)
# 	comm.send(poslen, dest=1, tag=2)

# 	g = buildGraphDirect(pointList, 0.9, poslen)

# 	print 'rank', rank, g.adjacency_list()

# elif rank == 1:
# 	pointList = comm.recv(source=0, tag=1)
# 	poslen = comm.recv(source=0, tag=2)
	
# 	g = buildGraphDirect(pointList, 0.9, poslen)
# 	print 'rank', rank, g.adjacency_list()


if rank == 0:
	pointList = generatePositiveList()
	poslen = len(pointList)
	addRandPoints(pointList, 500)

else:
	pointList = None
	poslen = None
	
pointList = comm.bcast(pointList, root=0)
poslen = comm.bcast(poslen, root=0)

#print 'rank', rank, pointList
g = buildGraphDirect(pointList, 0.9, poslen)


it = 1
gm = ParaPolysomeDetectionMRF(g)

if rank == 0:
	print 'init', time.ctime()

mapVec = gm.getMAP(it)

if rank == 0:
	print 'end', time.ctime()
	for i in range(len(pointList)):
		print i, ':', mapVec[i]




