
def confusionValues(trueSets, detectedSets, noiseSet):

    truePositive = 0.0
    trueNegative = 0.0
    falsePositive = 0.0
    falseNegative = 0.0

    detectedSet = set()

    for ps in trueSets:
        maxIntersection = 0
        setIndex = 0
        for j in range(len(detectedSets)):
            pd = detectedSets[j]
            intersectionSize = len(ps & pd)
            if intersectionSize > maxIntersection:
                maxIntersection = intersectionSize
                setIndex = j

        if len(detectedSets) > 0:
            pd = detectedSets[setIndex]
        else:
            pd = set()

        truePositive = truePositive + (len(ps & pd))
        falsePositive = falsePositive + (len(pd) - len(ps & pd))
        falseNegative = falseNegative + (len(ps) - len(ps & pd))

        detectedSet = detectedSet | pd ##

    trueNegative = float(len(noiseSet) - len(noiseSet & detectedSet))

    return truePositive, trueNegative, falsePositive, falseNegative

def getMeasuresA(truePositive, trueNegative, falsePositive, falseNegative):

    sensitivity = float(truePositive) / (truePositive + falseNegative)
    specificity = float(trueNegative) / (trueNegative + falsePositive)

    return sensitivity, specificity

def getMeasuresB(truePositive, trueNegative, falsePositive, falseNegative):

    accuracy = (truePositive + trueNegative)/(truePositive + falsePositive + trueNegative + falseNegative)
    precision = truePositive/(truePositive + falsePositive)
    recall = truePositive/(truePositive + falseNegative)

    return accuracy, precision, recall


def getSetsFromParticleLists(true_pl, detected_pl, minSize):

    classMap = {}
    polysomeSet = set()
    for i in range(len(true_pl)):

        p = true_pl[i]

        classId = int(float(p.getClass()))

        if classId < 0:
            continue

        if classMap.has_key(classId):
            classMap[classId].append(i)
        else:
            classMap[classId] = [i]

        polysomeSet.add(i)

    noiseSet = set(range(len(true_pl))) - polysomeSet

    trueSets = []
    for k in classMap.keys():
        polysomeSet = set(classMap[k])
        trueSets.append(polysomeSet)

    #################

    detectedSets = []

    polysomeMap = {}
    for i in range(len(detected_pl)):

        p = detected_pl[i]
        classId = int(float(p.getClass()))

        if polysomeMap.has_key(classId):
            polysomeMap[classId].append(i)
        else:
            polysomeMap[classId] = [i]


    for k in polysomeMap.keys():
        if len(polysomeMap[k]) >= minSize:

            detectedSets.append(set(polysomeMap[k]))

    return trueSets, detectedSets, noiseSet

