from polysomeDetection.polysomeDetectionStructures import PolysomeDetectionMRF
from polysomeDetection.sim.simFunctions import generatePositiveList, buildGraphRand, paintGraph
import numpy as np

posNum = 9
mapIt = 1

# neighborIt = 25
# for i in range(1, neighborIt +1):

# 	sr_sum = 0.0

# 	for j in range(10):

# 		pointList = generatePositiveList()
# 		g = buildGraphRand(pointList, i*posNum, 0.0)
# 		#paintGraph(pointList, g)

# 		gm = PolysomeDetectionMRF(g)
# 		mapVec = gm.getMAP(mapIt)

# 		ground = [0, 0, 2, 2, 2, 5, 5, 5, 5]
# 		ground.extend(range(9, 9 + i*posNum))
# 		n = len(ground)
# 		ground = np.array(ground)

# 		res = (mapVec == ground)
		
# 		sr = float(res.sum())/n

# 		sr_sum = sr_sum + sr

# 	sr_mean = sr_sum/10

# 	print 'nn', i*posNum, sr_mean

# 	f = open('stats_neighbors.txt', 'a')
# 	f.write(str( float(posNum)/(i*posNum) ) + ' ' + str(sr_mean) + '\n')
# 	f.close()

# print 'done...'
# ##########

sdIt = 500
sdInit = 0.05
for i in range(1, sdIt +1):

	sr_sum = 0.0

	for j in range(10):		

		pointList = generatePositiveList()
		g = buildGraphRand(pointList, 50, i*sdInit)
		#paintGraph(pointList, g)

		gm = PolysomeDetectionMRF(g)
		mapVec = gm.getMAP(mapIt)

		ground = [0, 0, 2, 2, 2, 5, 5, 5, 5]
		ground.extend(range(9, 9 + 50))
		n = len(ground)
		ground = np.array(ground)

		res = (mapVec == ground)

		sr = float(res.sum())/n

		sr_sum = sr_sum + sr

	sr_mean = sr_sum/10

	print 'sd', i*sdInit, sr_mean

	f = open('stats_sd_1.txt', 'a')
	f.write(str( i*sdInit ) + ' ' + str(sr_mean) + '\n')
	f.close()


print 'done...'
##########