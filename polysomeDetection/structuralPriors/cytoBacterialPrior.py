
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

class CytoBacterialPrior:
    
    def __init__(self, pixelSize):
        import numpy as np

        originalPixelSize = 2.24 #nm        
        pixelSizeFactor = originalPixelSize/float(pixelSize)
                
        templateCenter = np.array([30.0, 30.0, 30.0])
        pointInput = (np.array([32.0, 13.0, 29.0]) - templateCenter) ##entry
        pointOutput = (np.array([22.0, 17.0, 30.0]) - templateCenter) ##exit
        self._pointInput = 0.25*pointInput*pixelSizeFactor
        
        self._minW = 1e-20
        self._maxW = 1.0
                       
        ########

        self._vectorMean1 = np.array([-20.05566142, -12.15328273,   1.14399976])
        self._vectorCov1 = np.zeros((3, 3))
        self._vectorCov1[0, :] = [ 36.47814307, -25.50411012,   3.44513998]
        self._vectorCov1[1, :] = [-25.50411012,  86.29420872,  -9.8524489 ]
        self._vectorCov1[2, :] = [  3.44513998,  -9.8524489,   56.10919657]        
        self._vectorMean1 = pixelSizeFactor*0.25*self._vectorMean1
        self._vectorCov1 = pixelSizeFactor*0.25*self._vectorCov1

        self._vectorCov1 = 1.5*self._vectorCov1 #7.0

        #########
        #side t-t
        self._bingham1_mode_quat = np.array([ 0.24119865, -0.19825854,  0.64864188,  0.69410408])
        self._bingham1_A = 0.02*np.array([-64.0, -25.0, -16.0]) #0.05
        self._bingham1_V = np.array([[-0.16586544,  0.39872731, -0.86909631],
                                     [-0.53816316,  0.72545473,  0.38051202],
                                     [-0.64950901, -0.37399833,  0.13238958],
                                     [ 0.51088854,  0.41815932,  0.28697597]])
        self._bingham1_F = 1.0
        self._bingham1_mode = binghamPDF(self._bingham1_mode_quat, self._bingham1_A, self._bingham1_V, self._bingham1_F)

        #mid t-b       
        self._bingham2_mode_quat = np.array([ 0.90676364, -0.34780511, -0.23025406, -0.06159851])
        self._bingham2_A = 0.1*np.array([-11.56,  -6.76,  -2.56]) #0.2
        self._bingham2_V = np.array([[-0.31369709,  0.2546413,  -0.12054728],
                                     [-0.48079395,  0.76487835,  0.25065894],
                                     [-0.64429479, -0.29047608, -0.66894764],
                                     [ 0.50528759,  0.51550225, -0.68930944]])
        self._bingham2_F = 1.0
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
        weight = vectorComp1*(rotComp1 + rotComp2)

        #weight = vectorComp1
        
        weight = np.max((weight, self._minW))
        weight = np.min((weight, self._maxW))
                                
        return weight

