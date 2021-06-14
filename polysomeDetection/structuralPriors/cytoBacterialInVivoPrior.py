
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

class CytoBacterialInVivoPrior:
    
    def __init__(self, pixelSize):
        import numpy as np

        originalPixelSize = 2.74 #nm        
        pixelSizeFactor = originalPixelSize/float(pixelSize)
        
        self._minW = 1e-20
        self._maxW = 1.0
                       
        ########

        self._vectorMean1 = np.array([-9.26800775, -2.0953042,  -1.08114208])
        self._vectorCov1 = np.zeros((3, 3))
        self._vectorCov1[0, :] = [  4.00581011e+00,  -1.13841342e+00,  -2.71465418e-03]
        self._vectorCov1[1, :] = [ -1.13841342e+00,   7.39390925e+00,  -5.22888570e-01]
        self._vectorCov1[2, :] = [ -2.71465418e-03,  -5.22888570e-01,   7.36591519e+00]        
        self._vectorMean1 = pixelSizeFactor*self._vectorMean1
        self._vectorCov1 = pixelSizeFactor*self._vectorCov1

        self._vectorCov1 = 0.05*self._vectorCov1 #

        #########
        self._bingham1_mode_quat = np.array([-0.44699263, -0.19779147,  0.51982296,  0.70061417])
        self._bingham1_A = 0.2*np.array([-7.84, -4.84, -3.24]) #
        self._bingham1_V = np.array([[ 0.69672425,  0.4084237,  -0.38465957],
                                    [ 0.45358399, -0.19204187,  0.84750222],
                                    [ 0.52462612, -0.60543262, -0.29665277],
                                    [ 0.1833148,   0.65556185,  0.2139489 ]])
        self._bingham1_F = 1.0
        self._bingham1_mode = binghamPDF(self._bingham1_mode_quat, self._bingham1_A, self._bingham1_V, self._bingham1_F)

        #########
      
        self._bingham2_mode_quat = np.array([-0.33077208,  0.34683736, -0.57158108, -0.66602458])
        self._bingham2_A = 0.5*np.array([-9.,   -5.76, -4.84]) #0.2
        self._bingham2_V = np.array([[ 0.91616835, -0.21659938, -0.06565125],
                                    [ 0.14635185, -0.18572617,  0.9076292 ],
                                    [-0.36772464, -0.72192544,  0.12998966],
                                    [-0.06320814,  0.63040793,  0.39370272]])
        self._bingham2_F = 1.0
        self._bingham2_mode = binghamPDF(self._bingham2_mode_quat, self._bingham2_A, self._bingham2_V, self._bingham2_F)

        #########

        print 'using bac t, R'


                

    def vectorPrior1(self, output, input):
        import numpy as np
        from pytom.angles.angleFnc import pointRotateZXZ

        z1 = output[0]
        z2 = output[1]
        x1 = output[2]
        t = np.array([output[3], output[4], output[5]])

        t_in = np.array([input[3], input[4], input[5]])
        t_in = t_in - t
        t_in = np.array(pointRotateZXZ(t_in.tolist(), -z2, -z1, -x1))

        ##########################
        
        value = multivariateNormalPDF(t_in,
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
        
        inputData = [z1, z2, x1, x, y, z]
        
        vectorComp1 = self.vectorPrior1(outputData, inputData)
        rotComp1 = self.rotPrior1(outputData, inputData)
        rotComp2 = self.rotPrior2(outputData, inputData)
        
        weight = vectorComp1*(rotComp1 + rotComp2)

        #weight = vectorComp1
        #weight = rotComp1 + rotComp2
        
        weight = np.max((weight, self._minW))
        weight = np.min((weight, self._maxW))
                                
        return weight

