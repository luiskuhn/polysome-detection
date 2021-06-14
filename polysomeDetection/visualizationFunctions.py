
def paint2dGraph(pl, g):

    import matplotlib.pyplot as plt

    pointList = []
    min_x = pl[0].getPickPosition().getX()
    min_y = pl[0].getPickPosition().getY()
    max_x = pl[0].getPickPosition().getX()
    max_y = pl[0].getPickPosition().getY()

    for p in pl:

        x = p.getPickPosition().getX()
        y = p.getPickPosition().getY()
        pointList.append((x, y))

        if x < min_x:
            min_x = x
        if x > max_x:
            max_x = x
        if y < min_y:
            min_y = y
        if y > max_y:
            max_y = y


    pointArray = np.array(pointList)

    #plt.plot(pointArray[:,0], pointArray[:,1], 'ko', markersize=8)
    plt.plot(pointArray[:,0], pointArray[:,1], 'bo')

    max_w = 0.0
    for e in g.edges():
        if g[e[0]][e[1]]['weight'] > max_w:
            max_w = g[e[0]][e[1]]['weight']
    print 'max_w', max_w

    for e in g.edges():
        plt.plot([pointArray[e[0]][0], pointArray[e[1]][0]], [pointArray[e[0]][1], pointArray[e[1]][1]], 'b', alpha=g[e[0]][e[1]]['weight'])
        #(10000000*g[e[0]][e[1]]['weight']))
        #(0.5 + g[e[0]][e[1]]['weight'])) #(g[e[0]][e[1]]['weight']/max_w))

    #plt.savefig("edge_colormap.png") # save as png
    plt.ylim([min_y, max_y])
    plt.xlim([min_x, max_x])
 
    plt.show() # display

def paint3dGraph(pl, g):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import time

    pointList = []
    
    for p in pl:

        x = p.getPickPosition().getX()
        y = p.getPickPosition().getY()
        z = p.getPickPosition().getZ()
        pointList.append((x, y, z))
      
    pointArray = np.array(pointList)

    plt.ion()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot3D(pointArray[:,0], pointArray[:,1], pointArray[:,2], 'ko')
    for e in g.edges():
        w = g[e[0]][e[1]]['weight']
        ax.plot3D([pointArray[e[0]][0], pointArray[e[1]][0]], [pointArray[e[0]][1], pointArray[e[1]][1]], [pointArray[e[0]][2], pointArray[e[1]][2]],
            'b', color=(w, 0, 1.0 -w), alpha=0.2)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    for angle in range(0, 90):
        #time.sleep(0.0001)
        ax.view_init(30, angle)        
        plt.draw()

    plt.show(block=True)

def graphToMarkerFile(filePath, pl, g, binFactor=1, xOffset=0, yOffset=0, zOffset=0, useShifts=False, name='graph'):
    import math

    mf = open(filePath, 'w')
    mf.write('<marker_set name=\"' + name + '\">')
    mf.write('\n')

    pointList = []
    
    for i in range(len(pl)):
        p = pl[i]
        x = (p.getPickPosition().getX()*binFactor) + xOffset
        y = (p.getPickPosition().getY()*binFactor) + yOffset
        z = (p.getPickPosition().getZ()*binFactor) + zOffset
        if useShifts:
            x = x + (p.getShift().getX()*binFactor)
            y = y + (p.getShift().getY()*binFactor)
            z = z + (p.getShift().getZ()*binFactor)

        mf.write('<marker id=\"' + str(i) + '\" x=\"' + str(x) + '\" y=\"' + str(y) + '\" z=\"' + str(z) + '\"' + 
            ' r=\"0\" g=\"0\" b=\"1\" radius=\"1.0\"/>')
        mf.write('\n')

    ##
    min_w = 1.0
    max_w = 0.0
    for e in g.edges():
        s = e[0]
        t = e[1]
        w = g[s][t]['weight']
        if w < min_w:
            min_w = w
        if w > max_w:
            max_w = w

    links = []
    for e in g.edges():
        s = e[0]
        t = e[1]
        w = g[s][t]['weight']

        if (t, s) in links:
            continue

        if g[t].has_key(s):
            op_w = g[t][s]['weight']
            if op_w > w:
                w = op_w
        
        if w == min_w:
            w = 0.0

        w = w/max_w

        mf.write('<link id1=\"' + str(s) + '\" id2=\"' + str(t) + '\" r=\"' + str(w) + '\" g=\"0\" b=\"' + str(1.0 -w) + '\" radius="0.5"/>')
        mf.write('\n')

        links.append((s, t))

    mf.write('</marker_set>')
    mf.close()

def peptideExitsToBildFile(filePath, pl, polysomeId, pixelSizeFactor=1.0, binFactor=1, xOffset=0, yOffset=0, zOffset=0, useShifts=False, riboType='80S', coneRadius=1.0):
	import numpy as np
	from pytom.angles.angleFnc import pointRotateZXZ

	templateCenter = 0
	point3 = 0
	point4 = 0
	
	if riboType == '80S':
		templateCenter = np.array([10.0, 10.0, 10.0])
		point3 = (np.array([10.0, 15.0, 8.0]) - templateCenter)*pixelSizeFactor #peptide exit
		point4 = (np.array([10.0, 17.0, 6.0]) - templateCenter)*pixelSizeFactor # peptide exit top
	elif riboType == '70S':
		templateCenter = np.array([7.5, 7.5, 7.5])
		point3 = (np.array([7.0, 11.0, 6.0]) - templateCenter)*pixelSizeFactor #peptide exit
		point4 = (np.array([7.0, 14.0, 6.0]) - templateCenter)*pixelSizeFactor # peptide exit top

	polysomeList = []
	for i in pl._polysomeMap[polysomeId]:
		polysomeList.append(pl[i])

	f = open(filePath, 'w')
	f.write('.color 0 1 0')	
	f.write('\n')

	for p in polysomeList:
		x = (p.getPickPosition().getX()*binFactor) + xOffset
		y = (p.getPickPosition().getY()*binFactor) + yOffset
		z = (p.getPickPosition().getZ()*binFactor) + zOffset
		if useShifts:
			x = x + (p.getShift().getX()*binFactor)
			y = y + (p.getShift().getY()*binFactor)
			z = z + (p.getShift().getZ()*binFactor)

		t = np.array([x, y, z])
		z1 = p.getRotation().getZ1()
		z2 = p.getRotation().getZ2()
		x1 = p.getRotation().getX()

		pepExit = np.array(pointRotateZXZ(point3.tolist(), z1, z2, x1)) + t
		pepExitTop = np.array(pointRotateZXZ(point4.tolist(), z1, z2, x1)) + t

		f.write('.cone ' + str(pepExitTop[0]) + ' ' + str(pepExitTop[1]) + ' ' + str(pepExitTop[2]) + ' ' +
			str(pepExit[0]) + ' ' + str(pepExit[1]) + ' ' + str(pepExit[2]) + ' ' + str(coneRadius))
		f.write('\n')

	f.close()

def mRNA2MarkerFile(filePath, pl, polysomeId, pixelSizeFactor=1.0, binFactor=1, xOffset=0, yOffset=0, zOffset=0, useShifts=False, riboType='80S', markerRadius=0.5, linkRadius=0.2, knn=3, rnn=50.0):
	import numpy as np
	from pytom.angles.angleFnc import pointRotateZXZ
	from scipy.spatial import KDTree

	if knn < 3:
		knn = 3

	templateCenter = 0
	point1 = 0
	point2 = 0

	if riboType == '80S':
		templateCenter = np.array([10.0, 10.0, 10.0])
		point1 = (np.array([10.0, 5.0, 10.0]) - templateCenter)*pixelSizeFactor ##entry
		point2 = (np.array([8.0, 6.0, 9.0]) - templateCenter)*pixelSizeFactor ##exit
	elif riboType == '70S':
		templateCenter = np.array([7.5, 7.5, 7.5])
		point1 = (np.array([8.0, 3.0, 8.0]) - templateCenter)*pixelSizeFactor ##entry
		point2 = (np.array([4.0, 3.0, 8.0]) - templateCenter)*pixelSizeFactor ##exit

	polysomeList = []
	for i in pl._polysomeMap[polysomeId]:
		polysomeList.append(pl[i])

	f = open(filePath, 'w') #.cmm
	f.write('<marker_set name=\"mRNA ' + str(polysomeId) + '\">')
	f.write('\n')

	marker_id = 0
	treeData = []
	pointData = []

	for p in polysomeList:
		x = (p.getPickPosition().getX()*binFactor) + xOffset
		y = (p.getPickPosition().getY()*binFactor) + yOffset
		z = (p.getPickPosition().getZ()*binFactor) + zOffset
		if useShifts:
			x = x + (p.getShift().getX()*binFactor)
			y = y + (p.getShift().getY()*binFactor)
			z = z + (p.getShift().getZ()*binFactor)

		t = np.array([x, y, z])
		z1 = p.getRotation().getZ1()
		z2 = p.getRotation().getZ2()
		x1 = p.getRotation().getX()
		
		entryPoint = np.array(pointRotateZXZ(point1.tolist(), z1, z2, x1)) + t
		exitPoint = np.array(pointRotateZXZ(point2.tolist(), z1, z2, x1)) + t

		f.write('<marker id=\"' + str(marker_id) + '\" x=\"' + str(entryPoint[0]) + '\" y=\"' + str(entryPoint[1]) + '\" z=\"' + str(entryPoint[2]) + '\"' + 
			' r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(markerRadius) + '\"/>')
		f.write('\n')
		marker_id = marker_id + 1
		f.write('<marker id=\"' + str(marker_id) + '\" x=\"' + str(exitPoint[0]) + '\" y=\"' + str(exitPoint[1]) + '\" z=\"' + str(exitPoint[2]) + '\"' + 
			' r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(markerRadius) + '\"/>')
		f.write('\n')
		marker_id = marker_id + 1

		treeData.append((entryPoint[0], entryPoint[1], entryPoint[2]))
		pointData.append([entryPoint, exitPoint])

	tree = KDTree(list(treeData))

	for i in range(len(polysomeList)):

		entryPoint = pointData[i][0]		
		exitPoint = pointData[i][1]
		
		f.write('<link id1=\"' + str(i*2) + '\" id2=\"' + str((i*2) + 1) + '\" r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(linkRadius) + '\"/>')
		f.write('\n')

		query_tuple = tuple(exitPoint)

		#####

		nnd, nni = tree.query(query_tuple, k=knn, eps=0.0, p=2.0, distance_upper_bound=rnn)				
		nnd = list(nnd)
		nni = list(nni)

		nn = 0
		if nnd.count(float('inf')) > 0:
			infIndex = nnd.index(float('inf'))
			nn = nni[:infIndex]
		else:
			nn = nni

		#####

		nCount = 0	
		for j in nn:
			j = int(j)

			if j == i:
				continue

			f.write('<link id1=\"' + str((i*2) + 1) + '\" id2=\"' + str((j*2)) + '\" r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(linkRadius) + '\"/>')
			f.write('\n')

	f.write('</marker_set>')
	f.close()

def mRNA2MarkerFile2(filePath, pl, polysomeId, pixelSizeFactor=1.0, binFactor=1, xOffset=0, yOffset=0, zOffset=0, useShifts=False, riboType='80S', markerRadius=0.5, linkRadius=0.2, knn=3, rnn=50.0):
	import numpy as np
	from pytom.angles.angleFnc import pointRotateZXZ
	import networkx as nx
	
	templateCenter = 0
	point1 = 0
	point2 = 0

	if riboType == '80S':
		templateCenter = np.array([10.0, 10.0, 10.0])
		point1 = (np.array([10.0, 5.0, 10.0]) - templateCenter)*pixelSizeFactor ##entry
		point2 = (np.array([8.0, 6.0, 9.0]) - templateCenter)*pixelSizeFactor ##exit
	elif riboType == '70S':
		templateCenter = np.array([7.5, 7.5, 7.5])
		point1 = (np.array([8.0, 3.0, 8.0]) - templateCenter)*pixelSizeFactor ##entry
		point2 = (np.array([4.0, 3.0, 8.0]) - templateCenter)*pixelSizeFactor ##exit

	polysomeList = []
	polysomeSeq = pl._polysomeMap[polysomeId]
	for i in polysomeSeq:
		polysomeList.append(pl[i])

	subG = pl._G.subgraph(polysomeSeq)
	paths = nx.single_source_dijkstra_path(subG, polysomeId, weight='neg_log')

	f = open(filePath, 'w') #.cmm
	f.write('<marker_set name=\"mRNA ' + str(polysomeId) + '\">')
	f.write('\n')

	marker_id = 0
	
	for p in polysomeList:
		x = (p.getPickPosition().getX()*binFactor) + xOffset
		y = (p.getPickPosition().getY()*binFactor) + yOffset
		z = (p.getPickPosition().getZ()*binFactor) + zOffset
		if useShifts:
			x = x + (p.getShift().getX()*binFactor)
			y = y + (p.getShift().getY()*binFactor)
			z = z + (p.getShift().getZ()*binFactor)

		t = np.array([x, y, z])
		z1 = p.getRotation().getZ1()
		z2 = p.getRotation().getZ2()
		x1 = p.getRotation().getX()
		
		entryPoint = np.array(pointRotateZXZ(point1.tolist(), z1, z2, x1)) + t
		exitPoint = np.array(pointRotateZXZ(point2.tolist(), z1, z2, x1)) + t

		f.write('<marker id=\"' + str(marker_id) + '\" x=\"' + str(entryPoint[0]) + '\" y=\"' + str(entryPoint[1]) + '\" z=\"' + str(entryPoint[2]) + '\"' + 
			' r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(markerRadius) + '\"/>')
		f.write('\n')
		marker_id = marker_id + 1
		f.write('<marker id=\"' + str(marker_id) + '\" x=\"' + str(exitPoint[0]) + '\" y=\"' + str(exitPoint[1]) + '\" z=\"' + str(exitPoint[2]) + '\"' + 
			' r=\"1\" g=\"0.5\" b=\"0\" radius=\"' + str(markerRadius) + '\"/>')
		f.write('\n')
		marker_id = marker_id + 1
	
	edgeSet = set()
	for i in range(len(polysomeList)):
		
		f.write('<link id1=\"' + str(i*2) + '\" id2=\"' + str((i*2) + 1) + '\" r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(linkRadius) + '\"/>')
		f.write('\n')

		vertex_id = polysomeSeq[i]
		path = paths[vertex_id]

		# print vertex_id, path

		for j in range(len(path) -1):

			edgeTuple = (path[j], path[j +1])
			if not edgeTuple in edgeSet:

				s = polysomeSeq.index(edgeTuple[0])				
				t = polysomeSeq.index(edgeTuple[1])

				f.write('<link id1=\"' + str((s*2) + 1) + '\" id2=\"' + str((t*2)) + '\" r=\"1\" g=\"0\" b=\"0\" radius="0.2"/>')
				f.write('\n')

				edgeSet.add(edgeTuple)

	f.write('</marker_set>')
	f.close()

def mRNA2MarkerFile3(filePath, pl, polysomeId, entryPoint, exitPoint, pixelSizeFactor=1.0, binFactor=1, xOffset=0, yOffset=0, zOffset=0, useShifts=False, markerRadius=0.5, linkRadius=0.2, knn=3, rnn=50.0):
	import numpy as np
	from pytom.angles.angleFnc import pointRotateZXZ
	import networkx as nx
		
	point1 = entryPoint
	point2 = exitPoint
	
	polysomeList = []
	polysomeSeq = pl._polysomeMap[polysomeId]
	for i in polysomeSeq:
		polysomeList.append(pl[i])

	subG = pl._G.subgraph(polysomeSeq)
	paths = nx.single_source_dijkstra_path(subG, polysomeId, weight='neg_log')

	f = open(filePath, 'w') #.cmm
	f.write('<marker_set name=\"mRNA ' + str(polysomeId) + '\">')
	f.write('\n')

	marker_id = 0
	
	for p in polysomeList:
		x = (p.getPickPosition().getX()*binFactor) + xOffset
		y = (p.getPickPosition().getY()*binFactor) + yOffset
		z = (p.getPickPosition().getZ()*binFactor) + zOffset
		if useShifts:
			x = x + (p.getShift().getX()*binFactor)
			y = y + (p.getShift().getY()*binFactor)
			z = z + (p.getShift().getZ()*binFactor)

		t = np.array([x, y, z])
		z1 = p.getRotation().getZ1()
		z2 = p.getRotation().getZ2()
		x1 = p.getRotation().getX()
		
		entryPoint = np.array(pointRotateZXZ(point1.tolist(), z1, z2, x1)) + t
		exitPoint = np.array(pointRotateZXZ(point2.tolist(), z1, z2, x1)) + t

		f.write('<marker id=\"' + str(marker_id) + '\" x=\"' + str(entryPoint[0]) + '\" y=\"' + str(entryPoint[1]) + '\" z=\"' + str(entryPoint[2]) + '\"' + 
			' r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(markerRadius) + '\"/>')
		f.write('\n')
		marker_id = marker_id + 1
		f.write('<marker id=\"' + str(marker_id) + '\" x=\"' + str(exitPoint[0]) + '\" y=\"' + str(exitPoint[1]) + '\" z=\"' + str(exitPoint[2]) + '\"' + 
			' r=\"1\" g=\"0.5\" b=\"0\" radius=\"' + str(markerRadius) + '\"/>')
		f.write('\n')
		marker_id = marker_id + 1
	
	edgeSet = set()
	for i in range(len(polysomeList)):
		
		f.write('<link id1=\"' + str(i*2) + '\" id2=\"' + str((i*2) + 1) + '\" r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(linkRadius) + '\"/>')
		f.write('\n')

		vertex_id = polysomeSeq[i]
		path = paths[vertex_id]

		# print vertex_id, path

		for j in range(len(path) -1):

			edgeTuple = (path[j], path[j +1])
			if not edgeTuple in edgeSet:

				s = polysomeSeq.index(edgeTuple[0])				
				t = polysomeSeq.index(edgeTuple[1])

				f.write('<link id1=\"' + str((s*2) + 1) + '\" id2=\"' + str((t*2)) + '\" r=\"1\" g=\"0\" b=\"0\" radius="0.2"/>')
				f.write('\n')

				edgeSet.add(edgeTuple)

	f.write('</marker_set>')
	f.close()

def mRNA2MarkerFile4(filePath, pl, polysomeId, entry_Point, exit_Point, binFactor=1, xOffset=0, yOffset=0, zOffset=0, useShifts=False, markerRadius=0.5, linkRadius=0.2, knn=3, rnn=50.0):
	import numpy as np
	from pytom.angles.angleFnc import pointRotateZXZ
	from scipy.spatial import KDTree

	if knn < 3:
		knn = 3
	
	point1 = entry_Point
	point2 = exit_Point
	
	polysomeList = []
	for i in pl._polysomeMap[polysomeId]:
		polysomeList.append(pl[i])

	f = open(filePath, 'w') #.cmm
	f.write('<marker_set name=\"mRNA ' + str(polysomeId) + '\">')
	f.write('\n')

	marker_id = 0
	treeData = []
	pointData = []

	for p in polysomeList:
		x = (p.getPickPosition().getX()*binFactor) + xOffset
		y = (p.getPickPosition().getY()*binFactor) + yOffset
		z = (p.getPickPosition().getZ()*binFactor) + zOffset
		if useShifts:
			x = x + (p.getShift().getX()*binFactor)
			y = y + (p.getShift().getY()*binFactor)
			z = z + (p.getShift().getZ()*binFactor)

		t = np.array([x, y, z])
		z1 = p.getRotation().getZ1()
		z2 = p.getRotation().getZ2()
		x1 = p.getRotation().getX()
		
		entryPoint = np.array(pointRotateZXZ(point1.tolist(), z1, z2, x1)) + t
		exitPoint = np.array(pointRotateZXZ(point2.tolist(), z1, z2, x1)) + t

		f.write('<marker id=\"' + str(marker_id) + '\" x=\"' + str(entryPoint[0]) + '\" y=\"' + str(entryPoint[1]) + '\" z=\"' + str(entryPoint[2]) + '\"' + 
			' r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(markerRadius) + '\"/>')
		f.write('\n')
		marker_id = marker_id + 1
		f.write('<marker id=\"' + str(marker_id) + '\" x=\"' + str(exitPoint[0]) + '\" y=\"' + str(exitPoint[1]) + '\" z=\"' + str(exitPoint[2]) + '\"' + 
			' r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(markerRadius) + '\"/>')
		f.write('\n')
		marker_id = marker_id + 1

		treeData.append((entryPoint[0], entryPoint[1], entryPoint[2]))
		pointData.append([entryPoint, exitPoint])

	tree = KDTree(list(treeData))

	for i in range(len(polysomeList)):

		entryPoint = pointData[i][0]		
		exitPoint = pointData[i][1]
		
		f.write('<link id1=\"' + str(i*2) + '\" id2=\"' + str((i*2) + 1) + '\" r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(linkRadius) + '\"/>')
		f.write('\n')

		query_tuple = tuple(exitPoint)

		#####

		nnd, nni = tree.query(query_tuple, k=knn, eps=0.0, p=2.0, distance_upper_bound=rnn)				
		nnd = list(nnd)
		nni = list(nni)

		nn = 0
		if nnd.count(float('inf')) > 0:
			infIndex = nnd.index(float('inf'))
			nn = nni[:infIndex]
		else:
			nn = nni

		#####

		nCount = 0	
		for j in nn:
			j = int(j)

			if j == i:
				continue

			f.write('<link id1=\"' + str((i*2) + 1) + '\" id2=\"' + str((j*2)) + '\" r=\"1\" g=\"0\" b=\"0\" radius=\"' + str(linkRadius) + '\"/>')
			f.write('\n')

	f.write('</marker_set>')
	f.close()











