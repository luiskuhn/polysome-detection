
def rotationSum(a_z1, a_z2, a_x, b_z1, b_z2, b_x):
    from pytom.angles.angleFnc import zxzToMat, matToZXZ
        
    m1 = zxzToMat(a_z1, a_z2, a_x)
    m2 = zxzToMat(b_z1, b_z2, b_x)
    
    ms = m2*m1
    
    if ms[2,2] > 1:
        ms[2,2] = 1.0
    elif ms[2,2] < -1:
        ms[2,2] = -1.0
    
    zxz = matToZXZ(ms)
    
    z1 = zxz.getZ1()
    z2 = zxz.getZ2()
    x = zxz.getX()
    
    z1 = z1%360
    z2 = z2%360
    x = x%180
    
    return [z1, z2, x]

def binghamPDF(x, A, V, F):
    import numpy as np
    import math

    x = np.array(x)

    v1 = np.array(V[:,0])
    v2 = np.array(V[:,1])
    v3 = np.array(V[:,2])
    # d1 = np.dot(v1, x)
    # d2 = np.dot(v2, x)
    # d3 = np.dot(v3, x)
    d1 = v1[0]*x[0] + v1[1]*x[1] + v1[2]*x[2] + v1[3]*x[3]
    d2 = v2[0]*x[0] + v2[1]*x[1] + v2[2]*x[2] + v2[3]*x[3]
    d3 = v3[0]*x[0] + v3[1]*x[1] + v3[2]*x[2] + v3[3]*x[3]
    a1 = A[0]
    a2 = A[1]
    a3 = A[2]

    value = math.exp(a1*d1*d1 + a2*d2*d2 + a3*d3*d3) #/ F
        
    return value        

def euler2quat(z1, z2, x):
    from pytom.angles.quaternions import Quaternion
    
    q = Quaternion(z1, z2, x)
    
    return (q.getX(), q.getY(), q.getZ(), q.getW())

def multivariateNormalPDF(x, mean, cov):
    import numpy as np

    d = x.shape[0]
    p1 = np.exp(-0.5 * d * np.log(2 * np.pi))
    p2 = np.power(np.linalg.det(cov), -0.5)
    
    dev = x - mean
    
    p3 = np.exp(-0.5 * np.dot(np.dot(dev.transpose(), np.linalg.inv(cov)), dev))
    
    return p1*p2*p3

######################

class CytoEukarioticPrior:
    
    def __init__(self, pixelSize):
        import numpy as np

        originalPixelSize = 2.3 #nm        
        pixelSizeFactor = originalPixelSize/float(pixelSize)
                
        templateCenter = np.array([10.0, 10.0, 10.0])
        pointInput = (np.array([11.0, 5.0, 10.0]) - templateCenter) ##entry
        pointOutput = (np.array([7.0, 6.0, 9.0]) - templateCenter) ##exit
        self._pointInput = pointInput*pixelSizeFactor

        self._minW = 1e-30
        self._maxW = 1.0
                       
        ########

        self._vectorMean1 = np.array([-6.69137235, -4.95773035, -1.27676148])
        self._vectorCov1 = np.zeros((3, 3))
        self._vectorCov1[0, :] = [ 5.5287588,  -1.27370093, -0.70345909]
        self._vectorCov1[1, :] = [-1.27370093,  7.93021641, -0.43614245]
        self._vectorCov1[2, :] = [-0.70345909, -0.43614245,  9.50881794]        
        self._vectorMean1 = pixelSizeFactor*self._vectorMean1
        self._vectorCov1 = pixelSizeFactor*self._vectorCov1

        self._vectorCov1 = 0.1*self._vectorCov1 #1.0

        #########
        #cluster_0
        self._bingham1_mode_quat = np.array([ 0.59662938, -0.12447144, -0.48218198, -0.62931772])
        self._bingham1_A = 1000.0*np.array([-0.0576, -0.04,   -0.0256]) #1000.0
        self._bingham1_V = np.array([[-0.67692976,  0.23955312, -0.35834871],
                                     [-0.57557818, -0.12109165,  0.79909538],
                                     [-0.4030871,  -0.62690673, -0.46044465],
                                     [-0.21908139,  0.73139528, -0.14499481]])
        self._bingham1_F = 1.918781
        self._bingham1_mode = binghamPDF(self._bingham1_mode_quat, self._bingham1_A, self._bingham1_V, self._bingham1_F)

        #cluster_1       
        self._bingham2_mode_quat = np.array([-0.51409972,  0.3429061,  -0.54498523, -0.56666392])
        self._bingham2_A = 1000.0*np.array([-0.0676, -0.0484, -0.0361]) #1000.0
        self._bingham2_V = np.array([[-0.71836433, -0.1869535,   0.42977035],
                                     [-0.63955107, -0.06093326, -0.68532983],
                                     [ 0.27372425, -0.63699667, -0.47148847],
                                     [ 0.00146342,  0.7453662,  -0.35116818]])
        self._bingham2_F = 1.299603
        self._bingham2_mode = binghamPDF(self._bingham2_mode_quat, self._bingham2_A, self._bingham2_V, self._bingham2_F)

        #########

                

    def vectorPrior1(self, output, input):
        import numpy as np
        from pytom.angles.angleFnc import pointRotateZXZ

        z1 = output[0]
        z2 = output[1]
        x1 = output[2]
        t = np.array([output[3], output[4], output[5]])
                
        entryPoint = np.array(input[6])
        t_entryPoint = entryPoint - t
        templateEntryPoint = np.array(pointRotateZXZ(t_entryPoint.tolist(), -z2, -z1, -x1))

        ##########################
        
        value = multivariateNormalPDF(templateEntryPoint,
                                      self._vectorMean1,
                                      self._vectorCov1)/multivariateNormalPDF(self._vectorMean1, self._vectorMean1, self._vectorCov1)

        return value

    def rotPrior1(self, output, input):
        import numpy as np

        out_z1 = output[0]
        out_z2 = output[1]
        out_x1 = output[2]

        in_z1 = input[0]
        in_z2 = input[1]
        in_x1 = input[2]

        #rot = rotationSum(-out_z2, -out_z1, -out_x1, in_z1, in_z2, in_x1)
        rot = rotationSum(in_z1, in_z2, in_x1, -out_z2, -out_z1, -out_x1)
        
        q = np.array(euler2quat(rot[0], rot[1], rot[2]))

        value = binghamPDF(q, self._bingham1_A, self._bingham1_V, self._bingham1_F)/self._bingham1_mode

        return value

    def rotPrior2(self, output, input):
        import numpy as np

        out_z1 = output[0]
        out_z2 = output[1]
        out_x1 = output[2]

        in_z1 = input[0]
        in_z2 = input[1]
        in_x1 = input[2]

        #rot = rotationSum(-out_z2, -out_z1, -out_x1, in_z1, in_z2, in_x1)
        rot = rotationSum(in_z1, in_z2, in_x1, -out_z2, -out_z1, -out_x1)
        
        q = np.array(euler2quat(rot[0], rot[1], rot[2]))

        value = binghamPDF(q, self._bingham2_A, self._bingham2_V, self._bingham2_F)/self._bingham2_mode

        return value
     
    def weightEdge(self, outIndex, inIndex, pl, useShifts=False):
        from pytom.angles.angleFnc import pointRotateZXZ
        import numpy as np

        p = pl[outIndex]
        x = p.getPickPosition().getX()
        y = p.getPickPosition().getY()
        z = p.getPickPosition().getZ()
        if useShifts:
            x = x + p.getShift().getX()
            y = y + p.getShift().getY()
            z = z + p.getShift().getZ()
        z1 = p.getRotation().getZ1()
        z2 = p.getRotation().getZ2()
        x1 = p.getRotation().getX()
        outputData = [z1, z2, x1, x, y, z]
        
        p = pl[inIndex]
        x = p.getPickPosition().getX()
        y = p.getPickPosition().getY()
        z = p.getPickPosition().getZ()
        if useShifts:
            x = x + p.getShift().getX()
            y = y + p.getShift().getY()
            z = z + p.getShift().getZ()
        t = np.array([x, y, z])
        z1 = p.getRotation().getZ1()
        z2 = p.getRotation().getZ2()
        x1 = p.getRotation().getX()
        entryPoint = np.array(pointRotateZXZ(self._pointInput.tolist(), z1, z2, x1)) + t
        inputData = [z1, z2, x1, x, y, z, entryPoint]

        rotW1 = 0.5
        rotW2 = 1.0 -rotW1

        vectorComp1 = self.vectorPrior1(outputData, inputData)
        rotComp1 = self.rotPrior1(outputData, inputData)
        rotComp2 = self.rotPrior2(outputData, inputData)

        #weight = vectorComp1*(rotW1*rotComp1 + rotW2*rotComp2)
        # weight = vectorComp1*(rotComp1 + rotComp2)

        # weight = vectorComp1#*(rotComp1 + rotComp2)
        # weight = (rotComp1 + rotComp2)

        weight = vectorComp1
        
        weight = np.max((weight, self._minW))
        weight = np.min((weight, self._maxW))
                                
        return weight

