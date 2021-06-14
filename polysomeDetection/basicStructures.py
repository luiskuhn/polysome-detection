
from pytom.basic.structures import ParticleList
from polysomeDetection.postProcFunctions import makePolysomeMap, makeLenMap, transitivityCorrect
from polysomeDetection.postProcFunctions import inclusiveOrder, makeMaxBranchLenMap

from polysomeDetection.graphFunctions import buildGraph
from polysomeDetection.structuralPriors.cytoBacterialPrior import CytoBacterialPrior
from polysomeDetection.structuralPriors.cytoEukarioticPrior import CytoEukarioticPrior
from polysomeDetection.structuralPriors.erEukarioticPrior import ErEukarioticPrior


class PolysomeParticleList(ParticleList):

    def __init__(self, mapVector=None, G=None, pl=None, correctLoops=False, order='None', directory=None):
        
        ParticleList.__init__(self, directory, pl)

        import numpy as np

        if mapVector != None and G != None:            

            self._mapVector = np.array(mapVector, dtype = 'int')
            self._polysomeMap = makePolysomeMap(mapVector)

            if correctLoops:
            	self.correctLoops()

            self._lenMap = makeLenMap(self._polysomeMap)

            self._order = order            
            self.orderPolysomeSequences(G, order)


    def toXML(self):

        from lxml import etree

        directory_element = etree.Element('ParticleList')
        
        if len(self._particleList) > 0:
            for i in range(len(self._particleList)):
                particle_element = self._particleList[i].toXML()

                polysomeId = self._mapVector[i]
                polysomePos = -1

                if i in self._polysomeMap[self._mapVector[i]]:
                    polysomePos = self._polysomeMap[polysomeId].index(i)

                particle_element.append(etree.Element('Polysome', id = str(polysomeId), order = self._order, pos =str(polysomePos) ))
                directory_element.append(particle_element)
        
        return directory_element

    def fromXML(self,xmlObj):

    	from pytom.basic.structures import Particle
        from lxml.etree import _Element
        import numpy as np
        
        if xmlObj.__class__ != _Element :
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML object.')
    
        if xmlObj.tag == 'ParticleList':
            directory_element = xmlObj
        else:
            directory_element = xmlObj.xpath('ParticleList')
            directory_element = directory_element[0]
    
        self._particleList = []
        particles = directory_element.xpath('Particle')
    
        if len(particles) > 0:
            
            self._mapVector = np.zeros(len(particles), dtype = 'int')
            posVector = np.zeros(len(particles))            

            for i in range(len(particles)):
                p = particles[i]
                
                pp = Particle('')
                pp.fromXML(p)
                self._particleList.append(pp)

                particle_element = 0
                if p.tag == 'Particle':
                    particle_element = p
                else:
                    particle_element = p.xpath('Particle')
                    particle_element = particle_element[0]

                poly_element = particle_element.xpath('Polysome')
                poly_element = poly_element[0]
                polysomeId = int(str(poly_element.get('id')))
                polysomePos = int(str(poly_element.get('pos')))

                self._mapVector[i] = polysomeId
                posVector[i] = polysomePos

                if i == 0:
                    self._order = str(poly_element.get('order'))

            self._polysomeMap = makePolysomeMap(self._mapVector)
            self._lenMap = makeLenMap(self._polysomeMap)

            for k in self._polysomeMap.keys():

                polysomeList = self._polysomeMap[k]

                orderSize = 0
                for i in range(len(polysomeList)):
                    pos = posVector[polysomeList[i]]
                    if pos != -1:
                        orderSize = orderSize + 1

                polyVector = np.zeros(orderSize, dtype = 'int')

                for i in range(len(polysomeList)):
                    pos = posVector[polysomeList[i]]
                    if pos != -1:
                        polyVector[pos] = polysomeList[i]

                self._polysomeMap[k] = polyVector.tolist()

    def getPolysomesBySize(self, size):

        polysomeList = []

        if self._lenMap.has_key(size):
            idList = self._lenMap[size]
            for poly_id in idList:
                polysomeList.append(self._polysomeMap[poly_id])

        return polysomeList

    def correctLoops(self):

        transitivityCorrect(self._polysomeMap, self._mapVector)

    def orderPolysomeSequences(self, G, order):

        print 'hi'


class PolysomeParticleList2(ParticleList):

    def __init__(self, mapVector=None, G=None, pl=None, correctLoops=False, order='None', directory=None, polysomeType='None'):
        
        ParticleList.__init__(self, directory, pl)

        import numpy as np
        import copy

        self._mapVector = np.zeros(1)
        self._polysomeMap = {}
        self._lenMap = {}
        self._polysomeTypes = []
        self._particlePositions = []

        if pl != None and not type(pl) is PolysomeParticleList:
            self._mapVector = np.zeros(len(self._particleList))
            for i in range(len(self._particleList)):
                self._polysomeTypes.append(polysomeType)
                self._particlePositions.append(-1)
        elif type(pl) is PolysomeParticleList:
            self._mapVector = copy.deepcopy(pl._mapVector)
            self._polysomeMap = copy.deepcopy(pl._polysomeMap)
            self._lenMap = copy.deepcopy(pl._lenMap)
            self._polysomeTypes = copy.deepcopy(pl._polysomeTypes)
            self._particlePositions = copy.deepcopy(pl._particlePositions)
        elif mapVector != None and G != None:            

            self._mapVector = np.array(mapVector, dtype = 'int')
            self._polysomeMap = makePolysomeMap(mapVector)

            if correctLoops:
                self.correctLoops()

            self._lenMap = makeLenMap(self._polysomeMap)
            
            self.orderPolysomeSequences(G, order)

    def toXML(self):

        from lxml import etree

        directory_element = etree.Element('ParticleList')

        if len(self._particleList) > 0:
            for i in range(len(self._particleList)):
                particle_element = self._particleList[i].toXML()

                polysomeId = self._mapVector[i]
                polysomePos = self._particlePositions[i]

                particle_element.append(etree.Element('Polysome',
                    type = str(self._polysomeTypes[i]),
                    id = str(polysomeId), pos = str(polysomePos) ))

                directory_element.append(particle_element)

        return directory_element

    def fromXML(self, xmlObj):
        from pytom.basic.structures import Particle
        from lxml.etree import _Element
        import numpy as np

        self._mapVector = np.zeros(1)
        self._polysomeMap = {}
        self._lenMap = {}
        self._polysomeTypes = []
        self._particlePositions = []

        if xmlObj.__class__ != _Element :
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML object.')
    
        if xmlObj.tag == 'ParticleList':
            directory_element = xmlObj
        else:
            directory_element = xmlObj.xpath('ParticleList')
            directory_element = directory_element[0]
    
        self._particleList = []
        particles = directory_element.xpath('Particle')
    
        if len(particles) > 0:
            
            self._mapVector = np.zeros(len(particles), dtype = 'int')
            
            for i in range(len(particles)):
                p = particles[i]
                
                pp = Particle('')
                pp.fromXML(p)
                self._particleList.append(pp)

                particle_element = 0
                if p.tag == 'Particle':
                    particle_element = p
                else:
                    particle_element = p.xpath('Particle')
                    particle_element = particle_element[0]

                poly_element = particle_element.xpath('Polysome')
                poly_element = poly_element[0]


                polysomeId = int(str(poly_element.get('id')))
                polysomePos = int(str(poly_element.get('pos')))
                polysomeType = str(poly_element.get('type'))

                self._mapVector[i] = polysomeId
                self._particlePositions.append(polysomePos)
                self._polysomeTypes.append(polysomeType)

                self._particlePositions[i] = polysomePos

            self._polysomeMap = makePolysomeMap(self._mapVector)
            self._lenMap = makeLenMap(self._polysomeMap)

            #order polysomeLists?


    def getPolysomesBySize(self, size):

        polysomeList = []

        if self._lenMap.has_key(size):
            idList = self._lenMap[size]
            for poly_id in idList:
                polysomeList.append(self._polysomeMap[poly_id])

        return polysomeList

    def correctLoops(self):

        transitivityCorrect(self._polysomeMap, self._mapVector)

    def orderPolysomeSequences(self, G, order):
        for i in range(len(self._particleList)):
            polysomeId = self._mapVector[i]
            polysomePos = -1

            if i in self._polysomeMap[self._mapVector[i]]:
                polysomePos = self._polysomeMap[polysomeId].index(i)

            self._particlePositions[i] = polysomePos


class PolysomeParticleList3(ParticleList):

    def __init__(self, pl=None, mapVector=None, knn=None, rnn=None, ps=None, prior=None, priorType=None, useShifts=False, correctLoops=False, order=False, directory=None):
        
        ParticleList.__init__(self, directory, pl)

        import numpy as np
        import copy

        if pl == None:
            self._mapVector = 0
            self._polysomeMap = 0
            self._lenMap = 0
            self._particlePositions = 0
            self._knn = 0
            self._rnn = 0
            self._priorType = 0
            self._G = 0
            self._ps = 0
            self._shifts = useShifts
            self._maxBranchLenMap = 0

        else:
            self._mapVector = np.array(mapVector, dtype = 'int')
            self._polysomeMap = makePolysomeMap(mapVector)
            self._lenMap = {}
            self._particlePositions = []
            self._knn = knn
            self._rnn = rnn
            self._priorType = priorType
            self._ps = ps
            self._shifts = useShifts
            self._maxBranchLenMap = {}

            for i in range(len(self._particleList)):
                self._particlePositions.append(-1)

            if correctLoops:
                self.correctLoops()

            self._lenMap = makeLenMap(self._polysomeMap)

            self._G = buildGraph(pl, knn, rnn/ps, prior, useShifts=useShifts)

            if order:
                self.orderPolysomeSequences(self._G)

            self._maxBranchLenMap = makeMaxBranchLenMap(self._polysomeMap, self._G)

    def toXML(self):

        from lxml import etree

        directory_element = etree.Element('ParticleList',
                                            type = str(self._priorType),
                                            knn = str(self._knn),
                                            rnn = str(self._rnn),
                                            ps = str(self._ps),
                                            shifts = str(int(self._shifts)))

        if len(self._particleList) > 0:
            for i in range(len(self._particleList)):
                particle_element = self._particleList[i].toXML()

                polysomeId = self._mapVector[i]
                polysomePos = self._particlePositions[i]

                particle_element.append(etree.Element('Polysome', id = str(polysomeId), pos = str(polysomePos) ))

                directory_element.append(particle_element)

        return directory_element

    def fromXML(self, xmlObj):
        from pytom.basic.structures import Particle
        from lxml.etree import _Element
        import numpy as np

        self._mapVector = np.zeros(1)
        self._polysomeMap = {}
        self._lenMap = {}
        self._priorType = ''        
        self._particlePositions = []

        if xmlObj.__class__ != _Element :
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML object.')
    
        if xmlObj.tag == 'ParticleList':
            directory_element = xmlObj
        else:
            directory_element = xmlObj.xpath('ParticleList')
            directory_element = directory_element[0]

        self._priorType = str(directory_element.get('type'))
        self._knn = int(str(directory_element.get('knn')))
        self._rnn = float(str(directory_element.get('rnn')))
        self._ps = float(str(directory_element.get('ps')))
        self._shifts = bool(int(str(directory_element.get('shifts'))))
    
        self._particleList = []
        particles = directory_element.xpath('Particle')
    
        if len(particles) > 0:
            
            self._mapVector = np.zeros(len(particles), dtype = 'int')
            self._particlePositions = np.zeros(len(particles), dtype = 'int')
            
            for i in range(len(particles)):
                p = particles[i]
                
                pp = Particle('')
                pp.fromXML(p)
                self._particleList.append(pp)

                particle_element = 0
                if p.tag == 'Particle':
                    particle_element = p
                else:
                    particle_element = p.xpath('Particle')
                    particle_element = particle_element[0]

                poly_element = particle_element.xpath('Polysome')
                poly_element = poly_element[0]

                polysomeId = int(str(poly_element.get('id')))
                polysomePos = int(str(poly_element.get('pos')))
                
                self._mapVector[i] = polysomeId               
                self._particlePositions[i] = polysomePos

            self._polysomeMap = makePolysomeMap(self._mapVector)
            self._lenMap = makeLenMap(self._polysomeMap)

        prior = 0        
        if self._priorType == 'cytoBacterialPrior':
            prior = CytoBacterialPrior(self._ps)
        elif self._priorType == 'cytoEukarioticPrior':
            prior = CytoEukarioticPrior(self._ps)
        elif self._priorType == 'erEukarioticPrior':
            prior = ErEukarioticPrior(self._ps)

        # self._G = buildGraph(self._particleList, self._knn, self._rnn/self._ps, prior, useShifts=self._shifts)
        # self._maxBranchLenMap = makeMaxBranchLenMap(self._polysomeMap, self._G)

    def getPolysomesBySize(self, size):

        polysomeList = []

        if self._lenMap.has_key(size):
            idList = self._lenMap[size]
            for poly_id in idList:
                polysomeList.append(self._polysomeMap[poly_id])

        return polysomeList

    def getPolysomesByMaxBranchSize(self, size):

        polysomeList = []

        if self._maxBranchLenMap.has_key(size):
            idList = self._maxBranchLenMap[size]
            for poly_id in idList:
                polysomeList.append(self._polysomeMap[poly_id])

        return polysomeList

    def correctLoops(self):

        transitivityCorrect(self._polysomeMap, self._mapVector)

    def orderPolysomeSequences(self, G):

        for k in self._polysomeMap.keys():
            posList = inclusiveOrder(k, self._polysomeMap[k], G)

            for i in range(len(self._polysomeMap[k])):
                particle = self._polysomeMap[k][i]
                pos = posList[i]
                self._particlePositions[particle] = pos

    def buildGraph(self, prior):
        self._G = buildGraph(self._particleList, self._knn, self._rnn/self._ps, prior, useShifts=self._shifts)
        self._maxBranchLenMap = makeMaxBranchLenMap(self._polysomeMap, self._G)




class PolysomeParticleList4(ParticleList):

    def __init__(self, mapVector=None, G=None, pl=None, correctLoops=False, order='None', directory=None, polysomeType='None'):
        
        ParticleList.__init__(self, directory, pl)

        import numpy as np
        import copy

        self._mapVector = np.zeros(1)
        self._polysomeMap = {}
        self._lenMap = {}
        self._polysomeTypes = []
        self._particlePositions = []

        if pl != None and not type(pl) is PolysomeParticleList:
            self._mapVector = np.zeros(len(self._particleList))
            for i in range(len(self._particleList)):
                self._polysomeTypes.append(polysomeType)
                self._particlePositions.append(-1)
        elif type(pl) is PolysomeParticleList:
            self._mapVector = copy.deepcopy(pl._mapVector)
            self._polysomeMap = copy.deepcopy(pl._polysomeMap)
            self._lenMap = copy.deepcopy(pl._lenMap)
            self._polysomeTypes = copy.deepcopy(pl._polysomeTypes)
            self._particlePositions = copy.deepcopy(pl._particlePositions)
        elif mapVector != None and G != None:            

            self._mapVector = np.array(mapVector, dtype = 'int')
            self._polysomeMap = makePolysomeMap(mapVector)

            if correctLoops:
                self.correctLoops()

            self._lenMap = makeLenMap(self._polysomeMap)
            
            self.orderPolysomeSequences(G, order)

    def toXML(self):

        from lxml import etree

        directory_element = etree.Element('ParticleList',
                                            type = str(self._priorType),
                                            knn = str(self._knn),
                                            rnn = str(self._rnn),
                                            ps = str(self._ps),
                                            shifts = str(int(self._shifts)))

        if len(self._particleList) > 0:
            for i in range(len(self._particleList)):
                particle_element = self._particleList[i].toXML()

                polysomeId = self._mapVector[i]
                polysomePos = self._particlePositions[i]

                particle_element.append(etree.Element('Polysome', id = str(polysomeId), pos = str(polysomePos) ))

                directory_element.append(particle_element)

        return directory_element

    def fromXML(self, xmlObj):
        from pytom.basic.structures import Particle
        from lxml.etree import _Element
        import numpy as np

        self._mapVector = np.zeros(1)
        self._polysomeMap = {}
        self._lenMap = {}
        self._polysomeTypes = []
        self._particlePositions = []

        if xmlObj.__class__ != _Element :
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML object.')
    
        if xmlObj.tag == 'ParticleList':
            directory_element = xmlObj
        else:
            directory_element = xmlObj.xpath('ParticleList')
            directory_element = directory_element[0]
    
        self._particleList = []
        particles = directory_element.xpath('Particle')
    
        if len(particles) > 0:
            
            self._mapVector = np.zeros(len(particles), dtype = 'int')
            
            for i in range(len(particles)):
                p = particles[i]
                
                pp = Particle('')
                pp.fromXML(p)
                self._particleList.append(pp)

                particle_element = 0
                if p.tag == 'Particle':
                    particle_element = p
                else:
                    particle_element = p.xpath('Particle')
                    particle_element = particle_element[0]

                poly_element = particle_element.xpath('Polysome')
                poly_element = poly_element[0]


                polysomeId = int(str(poly_element.get('id')))
                polysomePos = int(str(poly_element.get('pos')))
                polysomeType = str(poly_element.get('type'))

                self._mapVector[i] = polysomeId
                self._particlePositions.append(polysomePos)
                self._polysomeTypes.append(polysomeType)

                self._particlePositions[i] = polysomePos

            self._polysomeMap = makePolysomeMap(self._mapVector)
            self._lenMap = makeLenMap(self._polysomeMap)

            #order polysomeLists?

    def getPolysomesBySize(self, size):

        polysomeList = []

        if self._lenMap.has_key(size):
            idList = self._lenMap[size]
            for poly_id in idList:
                polysomeList.append(self._polysomeMap[poly_id])

        return polysomeList

    def correctLoops(self):

        transitivityCorrect(self._polysomeMap, self._mapVector)

    def orderPolysomeSequences(self, G):

        for k in self._polysomeMap.keys():
            posList = inclusiveOrder(k, self._polysomeMap[k], G)

            for i in range(len(self._polysomeMap[k])):
                particle = self._polysomeMap[k][i]
                pos = posList[i]
                self._particlePositions[particle] = pos

    def setup(self, priorType, knn, rnn, ps, shifts):

        self._knn = knn
        self._rnn = rnn
        self._priorType = priorType
        self._ps = ps
        self._shifts = shifts

        prior = 0        
        if self._priorType == 'cytoBacterialPrior':
            prior = CytoBacterialPrior(ps)
        elif self._priorType == 'cytoEukarioticPrior':
            prior = CytoEukarioticPrior(ps)
        elif self._priorType == 'erEukarioticPrior':
            prior = ErEukarioticPrior(ps)

        self._G = buildGraph(self._particleList, knn, rnn/ps, prior, useShifts=shifts)

        self.orderPolysomeSequences(self._G)











