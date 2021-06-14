
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
    x = x%360
    
    return [z1, z2, x]

def multivariateNormalPDF(x, mean, cov):
    import numpy as np

    d = x.shape[0]
    p1 = np.exp(-0.5 * d * np.log(2 * np.pi))
    p2 = np.power(np.linalg.det(cov), -0.5)
    
    dev = x - mean
    
    p3 = np.exp(-0.5 * np.dot(np.dot(dev.transpose(), np.linalg.inv(cov)), dev))
    
    return p1*p2*p3

def von_mises_pdf_scaled(theta, theta_mean, m):
    import math

    rad_theta = math.radians(theta)
    rad_theta_mean = math.radians(theta_mean)

    value = math.exp(m*math.cos(rad_theta - rad_theta_mean))/math.exp(m)

    return value

#################################################################################

class ErEukarioticPrior:
    
    def __init__(self, pixelSize):
        import numpy as np
        import math

        originalPixelSize = 2.3 #nm        
        pixelSizeFactor = originalPixelSize/float(pixelSize)
        #pixelSizeFactor = 1.0

        templateCenter = np.array([10.0, 10.0, 10.0])
        pointInput = (np.array([11.0, 5.0, 10.0]) - templateCenter) ##entry
        pointOutput = (np.array([7.0, 6.0, 9.0]) - templateCenter) ##exit
        self._pointInput = pointInput*pixelSizeFactor

        self._minW = 1e-30
        self._maxW = 1.0
        
        #########

        # self._vectorMean = np.array([-12.12809919, -3.80563887,  0.23996123])
        # self._vectorCov = np.zeros((3, 3))
        # self._vectorCov[0, :] = [ 3.74897312,-1.50098483, 0.02833292]
        # self._vectorCov[1, :] = [-1.50098483, 5.76472403, 4.26525181]
        # self._vectorCov[2, :] = [ 0.02833292, 4.26525181, 9.30667892]        
        # self._vectorMean = pixelSizeFactor*self._vectorMean
        # self._vectorCov = pixelSizeFactor*self._vectorCov

        # self._vectorCov = 0.5*self._vectorCov #1.5

        self._vectorMean = np.array([-8.91149144, -3.72384487, 0.6975887])
        self._vectorCov = np.zeros((3, 3))
        self._vectorCov[0, :] = [ 1.85014052, -0.86921626, -0.55161837]
        self._vectorCov[1, :] = [-0.86921626,  5.24707746,  5.80015939]
        self._vectorCov[2, :] = [-0.55161837,  5.80015939, 10.82075923]        
        self._vectorMean = pixelSizeFactor*self._vectorMean
        self._vectorCov = pixelSizeFactor*self._vectorCov

        self._vectorCov = 1.5*self._vectorCov #1.5

        #########

        self._thetaMean = 330.0
        self._thetaSd = 25.0 #20.0
        self._thetaM = 1.0/(math.radians(self._thetaSd)**2)

        ########

        self._pointDown = (np.array([10.0, 15.0, 8.0]) - templateCenter) ##down
        self._pointUp = (np.array([10.0, 3.0, 16.0]) - templateCenter) ##up

        self._normalCutoff = 30.0 #40.0


    def vectorPrior(self, output, input):
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
                                      self._vectorMean,
                                      self._vectorCov)/multivariateNormalPDF(self._vectorMean, self._vectorMean, self._vectorCov)

        return value

    def rotPrior(self, output, input):
        from scipy.stats import norm

        plane_rot_s = [0, 0, 123]
        
        out_z1 = output[0]
        out_z2 = output[1]
        out_x1 = output[2]

        p_r = rotationSum(-out_z2, -out_z1, -out_x1, plane_rot_s[0], plane_rot_s[1], plane_rot_s[2])
        p_r = [-p_r[1], -p_r[0], -p_r[2]]

        in_z1 = input[0]
        in_z2 = input[1]
        in_x1 = input[2]

        n_r = rotationSum(-in_z2, -in_z1, -in_x1, plane_rot_s[0], plane_rot_s[1], plane_rot_s[2])
        n_r = [-n_r[1], -n_r[0], -n_r[2]]

        #rot = rotationSum(in_z1, in_z2, in_x1, -out_z2, -out_z1, -out_x1)
        rot = rotationSum(n_r[0], n_r[1], n_r[2], -p_r[1], -p_r[0], -p_r[2])
        sz1 = rot[0]%360
        sz2 = rot[1]%360
        
        ###############

        theta = (sz1 + sz2)%360
        value = von_mises_pdf_scaled(theta, self._thetaMean, self._thetaM)

        return value

    def weightEdge(self, outIndex, inIndex, pl, useShifts=False):
        from pytom.angles.angleFnc import pointRotateZXZ
        from numpy import linalg as LA
        import numpy as np
        import math

        weight = self._minW
        
        #################################################
        
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

        t = np.array([x, y, z])
        in_downPoint = np.array(pointRotateZXZ(self._pointDown.tolist(), z1, z2, x1)) + t
        in_upPoint = np.array(pointRotateZXZ(self._pointUp.tolist(), z1, z2, x1)) + t

        #################################################

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

        t = np.array([x, y, z])
        out_downPoint = np.array(pointRotateZXZ(self._pointDown.tolist(), z1, z2, x1)) + t
        out_upPoint = np.array(pointRotateZXZ(self._pointUp.tolist(), z1, z2, x1)) + t

        #################################################

        #filter by membrane normal
        #################################################
        normalFlag = True 

        downPoint = out_downPoint
        upPoint = out_upPoint

        n_downPoint = in_downPoint
        n_upPoint = in_upPoint
        
        t_downPoint = downPoint - t
        t_upPoint = upPoint - t
        t_n_downPoint = n_downPoint - t
        t_n_upPoint = n_upPoint - t
        r_downPoint = np.array(pointRotateZXZ(t_downPoint.tolist(), -z2, -z1, -x1))
        r_upPoint = np.array(pointRotateZXZ(t_upPoint.tolist(), -z2, -z1, -x1))
        r_n_downPoint = np.array(pointRotateZXZ(t_n_downPoint.tolist(), -z2, -z1, -x1))
        r_n_upPoint = np.array(pointRotateZXZ(t_n_upPoint.tolist(), -z2, -z1, -x1))
        
        vectorA = r_upPoint - r_downPoint
        vectorB = r_n_upPoint - r_n_downPoint
        vectorA = vectorA/LA.norm(vectorA)
        vectorB = vectorB/LA.norm(vectorB)

        x = np.inner(vectorA, vectorB)
        if x < -1:
            x = -1
        elif x > 1:
            x = 1
        x = math.degrees(math.acos(x))
        if x > self._normalCutoff:
            normalFlag = False

        #################################################

        if normalFlag:
            vectorComp = self.vectorPrior(outputData, inputData)
            rotComp = self.rotPrior(outputData, inputData)

            weight = vectorComp*rotComp

            weight = np.max((weight, self._minW))
            weight = np.min((weight, self._maxW))

                                
        return weight

