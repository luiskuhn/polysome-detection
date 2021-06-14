
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

class CytoEukarioticMRNAPrior:
    
    def __init__(self, pixelSize):
        import numpy as np

        originalPixelSize = 2.1 #nm        
        pixelSizeFactor = originalPixelSize/float(pixelSize)
                
        templateCenter = np.array([11.0, 11.0, 11.0])
        point_mRNA = (np.array([10.0, 6.0, 11.0]) - templateCenter) ##MRNA
        self._point_mRNA = point_mRNA*pixelSizeFactor

        self._minW = 1e-30
        self._maxW = 1.0
                       
        ########

        self._vectorMean1 = np.array([-6.14650148, -4.28081272, -2.91793258])
        self._vectorCov1 = np.zeros((3, 3))
        self._vectorCov1[0, :] = [ 2.23378593, -0.73153362,  0.1821682 ]
        self._vectorCov1[1, :] = [-0.73153362,  3.09091353, -0.68683336]
        self._vectorCov1[2, :] = [ 0.1821682,  -0.68683336,  2.47298855]       
        self._vectorMean1 = pixelSizeFactor*self._vectorMean1
        self._vectorCov1 = pixelSizeFactor*self._vectorCov1

        self._vectorCov1 = 1.0*self._vectorCov1 #1.0



        #########
        #cluster_0
        self._bingham1_mode_quat = np.array([ 0.37985225, -0.42835177,  0.61285126,  0.54464701])
        self._bingham1_A = 1.0*np.array([-25.,   -16.,    -6.76]) #0.1
        self._bingham1_V = np.array([[-0.196459,    0.38589442, -0.81743601],
                                    [-0.42914233,  0.75183927,  0.25901608],
                                    [-0.67847745, -0.1654781,   0.36972785],
                                    [ 0.56294678,  0.50837015,  0.35778535]])
        self._bingham1_F = 1.0
        self._bingham1_mode = binghamPDF(self._bingham1_mode_quat, self._bingham1_A, self._bingham1_V, self._bingham1_F)

        #########
        #cluster_1       
        self._bingham2_mode_quat = np.array([-0.19245504, -0.3688796,   0.56778855,  0.7102852 ])
        self._bingham2_A = 1.0*np.array([-36., -16.,  -9.]) #0.1
        self._bingham2_V = np.array([[-0.75263739,  0.59748914, -0.19875802],
                                    [ 0.3386434,   0.53430042,  0.68100775],
                                    [ 0.42705248,  0.59675108, -0.373002  ],
                                    [-0.36943713, -0.03765545,  0.59799098]])
        self._bingham2_F = 1.0
        self._bingham2_mode = binghamPDF(self._bingham2_mode_quat, self._bingham2_A, self._bingham2_V, self._bingham2_F)

        #########
        #cluster_2       
        self._bingham3_mode_quat = np.array([ 0.72847174, -0.34938671, -0.13177843, -0.57436251])
        self._bingham3_A = 1.0*np.array([-30.25, -12.96,  -6.76]) #0.1
        self._bingham3_V = np.array([[-0.26903031,  0.13930369, -0.6144478 ],
                                    [-0.8520159,   0.3870764,   0.04657991],
                                    [-0.36552127, -0.89943638, -0.20010711],
                                    [ 0.26093256,  0.14758234, -0.76173576]])
        self._bingham3_F = 1.0
        self._bingham3_mode = binghamPDF(self._bingham2_mode_quat, self._bingham2_A, self._bingham2_V, self._bingham2_F)

        #########

        print 'using mRNA'

   
    def vectorPrior1(self, output, input):
        import numpy as np
        from pytom.angles.angleFnc import pointRotateZXZ

        z1 = output[0]
        z2 = output[1]
        x1 = output[2]
        t = np.array([output[3], output[4], output[5]])
                
        in_t = np.array(input[6])
        in_t = in_t - t
        in_t = np.array(pointRotateZXZ(in_t.tolist(), -z2, -z1, -x1))

        ##########################
        
        value = multivariateNormalPDF(in_t,
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

    def rotPrior3(self, output, input):
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

        value = binghamPDF(q, self._bingham3_A, self._bingham3_V, self._bingham3_F)/self._bingham3_mode

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
        entryPoint = np.array(pointRotateZXZ(self._point_mRNA.tolist(), z1, z2, x1)) + t
        
        inputData = [z1, z2, x1, x, y, z, entryPoint]

        vectorComp1 = self.vectorPrior1(outputData, inputData)
        rotComp1 = self.rotPrior1(outputData, inputData)
        rotComp2 = self.rotPrior2(outputData, inputData)
        rotComp3 = self.rotPrior3(outputData, inputData)


        weight = vectorComp1*(rotComp1 + rotComp2 + rotComp3)

        #weight = (rotComp1 + rotComp2 + rotComp3)
        #weight = vectorComp1
        #weight = rotComp2
        
        weight = np.max((weight, self._minW))
        weight = np.min((weight, self._maxW))
                                
        return weight
