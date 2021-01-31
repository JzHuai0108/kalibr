from random import gauss
import sys
import numpy as np
import matplotlib.pyplot as plt

import aslam_backend as aopt
import aslam_cv as acv
import aslam_cameras_april as acv_april
import sm
import kalibr_common as kc

import BSplineIO

def removeUncoveredTimes(timeArray, timeShift):
    length = timeArray.shape[0]      
    if timeShift > 0:
        finishIndex = length-1
        finishTime = timeArray[-1] - timeShift
        for jack in range(length-1, -1, -1):
            if timeArray[jack] < finishTime:
                finishIndex = jack
                break
        if finishIndex == length-1 or finishIndex == 0:
            raise ValueError('finishIndex is at the %d-th element of the timeArray of length %d though timeShift is %.9f' % (finishIndex, length, timeShift))
        return 0, finishIndex+1
    elif timeShift < 0:
        startIndex = 0
        startTime = timeArray[0] - timeShift
        for jack in range(0, length, 1):
            if timeArray[jack] > startTime:
                startIndex = jack
                break
        if startIndex == length-1 or startIndex == 0:
            raise ValueError('startIndex is at the %d-th element of the timeArray of length %d though timeShift is %.9f' % (startIndex, length, timeShift))
        return startIndex, length
    else:
        return 0, length


def getCameraPoseAt(timeScalar, poseSplineDv, T_b_c=sm.Transformation()):
    timeOffsetPadding = 0.0
    dv = aopt.Scalar(timeScalar)
    timeExpression = dv.toExpression()
        
    if timeScalar <= poseSplineDv.spline().t_min() or timeScalar >= poseSplineDv.spline().t_max():
        print "Warn: %.9f time out of range [%.9f, %.9f]" % (timeScalar, poseSplineDv.spline().t_min(), poseSplineDv.spline().t_max())
        return sm.Transformation(), False
       
    T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
    sm_T_w_c = sm.Transformation(T_w_b.toTransformationMatrix())*T_b_c  
    return sm_T_w_c, True


class RsVisualInertialMeasViaBsplineSimulator(object):
    '''simulate visual(rolling shutter) inertial measurements with a predefined BSpline model representing a realistic motion'''
    def __init__(self, modelCoeffFile, targetYaml= "april_6x6.yaml", \
      chainYaml='camchain.yaml', time_offset=0.0, time_readout=30e-3):
        self.modelFilename = modelCoeffFile
        self.simulPoseSplineDv = BSplineIO.loadBsplineModel(modelCoeffFile)      
        refStateFile = modelCoeffFile.replace("knot_coeff", "states", 1)
        fullRefStates = BSplineIO.loadArrayWithHeader(refStateFile)

        if fullRefStates.shape[1] != 20:
            raise ValueError('Each row of states is expected to have 20 columns rather than %d' % fullRefStates.shape[1])

        fullRefStateTimes = fullRefStates[:, 2]
        if time_offset > 0:
            timeShift = (time_offset + time_readout)*1.2 
        else:
            timeShift = (time_offset - time_readout)*1.2
        startIndex, finishIndex = removeUncoveredTimes(fullRefStateTimes, timeShift)
       
        self.refStateTimes = fullRefStateTimes[startIndex:finishIndex]
        print >> sys.stdout, "  Camera frame sampling start time %.9f and finish time %.9f" % (self.refStateTimes[0], self.refStateTimes[-1])
        pathStub = ''
        if self.modelFilename.find('/') != -1:
            pathStub = self.modelFilename.rsplit('/', 1)[0]   
        
        trimmedStateFile = pathStub + "/initial_states_td" + str(time_offset).replace('.','',1) + '.txt'       
        trimmedStateStream = open(trimmedStateFile, 'w')
        for iota in range(startIndex, finishIndex, 1):
            print >> trimmedStateStream, '%d %.9f' % (iota-startIndex, fullRefStates[iota, 2]), \
              ' '.join(map(str, fullRefStates[iota, 3:19])), '%d' % fullRefStates[iota, 19]
        trimmedStateStream.close()
        print "  Written simulated states to", trimmedStateFile
        print "  Reading camera chain:", chainYaml
        chain = kc.CameraChainParameters(chainYaml)        
        camNr=0
        self.T_c0_imu = chain.getExtrinsicsImuToCam(camNr)  
        self.T_imu_c0 = self.T_c0_imu.inverse()
        camConfig = chain.getCameraParameters(camNr)
        camConfig.printDetails()
        resolution = camConfig.getResolution()
        print '  Camera resolution', resolution
        camera = kc.AslamCamera.fromParameters( camConfig)    
        self.camGeometry = camera.geometry

        targetConfig = kc.CalibrationTargetParameters(targetYaml)
        print("Target used in the simulation:")
        targetConfig.printDetails()
        self.setupCalibrationTarget(targetConfig, showExtraction=False, showReproj=False, imageStepping=False)
        self.imageHeight = resolution[1]  
        self.timeOffset = time_offset 
        self.lineDelay = time_readout/self.imageHeight
        print("Line delay {}, camera time offset {}".format(self.lineDelay, self.timeOffset))
        self.landmark2ObservationList = dict()
        for iota in range(self.simulatedObs.getTotalTargetPoint()):
            self.landmark2ObservationList[iota]=list()

    def setupCalibrationTarget(self, targetConfig, showExtraction=False, showReproj=False, imageStepping=False):
        '''copied from IccCamera class'''
        #load the calibration target configuration
        targetParams = targetConfig.getTargetParams()
        targetType = targetConfig.getTargetType()
    
        if targetType == 'checkerboard':
            options = acv.CheckerboardOptions() 
            options.filterQuads = True
            options.normalizeImage = True
            options.useAdaptiveThreshold = True        
            options.performFastCheck = False
            options.windowWidth = 5
            options.showExtractionVideo = showExtraction
            grid = acv.GridCalibrationTargetCheckerboard(targetParams['targetRows'], 
                                                            targetParams['targetCols'], 
                                                            targetParams['rowSpacingMeters'], 
                                                            targetParams['colSpacingMeters'],
                                                            options)
        elif targetType == 'circlegrid':
            options = acv.CirclegridOptions()
            options.showExtractionVideo = showExtraction
            options.useAsymmetricCirclegrid = targetParams['asymmetricGrid']
            grid = acv.GridCalibrationTargetCirclegrid(targetParams['targetRows'],
                                                          targetParams['targetCols'], 
                                                          targetParams['spacingMeters'], 
                                                          options)
        elif targetType == 'aprilgrid':
            options = acv_april.AprilgridOptions() 
            options.showExtractionVideo = showExtraction
            options.minTagsForValidObs = int( np.max( [targetParams['tagRows'], targetParams['tagCols']] ) + 1 )
            
            grid = acv_april.GridCalibrationTargetAprilgrid(targetParams['tagRows'],
                                                            targetParams['tagCols'], 
                                                            targetParams['tagSize'], 
                                                            targetParams['tagSpacing'], 
                                                            options)
        else:
            raise RuntimeError( "Unknown calibration target." )
                          
        options = acv.GridDetectorOptions() 
        options.imageStepping = imageStepping
        options.plotCornerReprojection = showReproj
        options.filterCornerOutliers = True
        
        self.simulatedObs = acv.GridCalibrationTargetObservation(grid) #this line should work although 
        # in the constructor GridCalibrationTargetObservation(GridCalibrationTarget::Ptr) is
        # inconsistent with GridCalibrationTargetBase::Ptr
        self.allTargetCorners = np.array(self.simulatedObs.getAllCornersTargetFrame()) # nx3
        assert self.allTargetCorners.shape[0] == self.simulatedObs.getTotalTargetPoint()
   
    def checkNaiveVsNewtonRsProjection(self):
        for frameId in range(1):
            state_time = self.refStateTimes[0]
            line_delay = self.lineDelay             
            imageCornersNaive = self.naiveMethodToRsProjection(state_time, line_delay, False)
            imageCornersNewton, unusedKeypoints, unusedImageOffset = \
                    self.newtonMethodToRsProjection(state_time, line_delay, 1.0, 0.0, False)

        assert np.allclose(imageCornersNaive[:,:,0], imageCornersNewton[:,:,0])
        reproducedImageCornerFile = self.modelFilename.replace("knot_coeff", "naive_vs_newton", 1)
        np.savetxt(reproducedImageCornerFile, np.concatenate((imageCornersNaive[:,:,0], imageCornersNewton[:,:,0]), axis=1), \
                fmt=['%.9f', '%.9f', '%d', '%.9f', '%.9f', '%d'])

    def simulate(self, landmarkOutputFilename, reprojectionSigma):
        '''simulate camera observations for frames at all ref state times and plus noise''' 
        print('Reprojection Gaussian noise sigma {:.3f}'.format(reprojectionSigma))
        imageCornerOffsetNorms = list()
        bins = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, \
                3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0]
        for frameId in range(0, self.refStateTimes.shape[0], 1):
             state_time = self.refStateTimes[frameId]
             line_delay = self.lineDelay
             unusedCornerProjections, frameKeypoints, imageCornerOffset = \
                    self.newtonMethodToRsProjection(state_time, line_delay, reprojectionSigma, self.timeOffset)
             imageCornerOffsetNorms += imageCornerOffset
             if frameId % 300 == 0:
                print('  Projected {:4d} target corners for state at {:.9f}'.format(len(frameKeypoints), state_time))
             for lmObs in frameKeypoints:
                 self.landmark2ObservationList[lmObs[0]].append((frameId, lmObs[1], lmObs[2], lmObs[3], lmObs[4]))
        lmStream = open(landmarkOutputFilename, 'w')
        lmStream.write(('#Each line contains a landmark: lm id, lm p_w(homogeneous xyzw), quality, num obs, '
                        '<frameid> <keypoint id, x, y, size>, ... , <frameid> <keypoint id, x, y, size>\n'))
        probe = 0
        homogeneousW = 1.0
        positionQuality = 1.0
        for lm, obsList in sorted(self.landmark2ObservationList.iteritems()):
            assert probe == lm        
            lmStream.write('%d ' % lm + ' '.join(map(str, self.allTargetCorners[lm,:]))+ \
              ' %.3f %.3f %d' % (homogeneousW, positionQuality, len(obsList)))
            for observation in obsList:
                lmStream.write(' '+' '.join(map(str, observation)))
            lmStream.write('\n')
            probe += 1
        lmStream.close()
        print '  Written landmark observations to', landmarkOutputFilename
        print '  Histogram of the offset in the infinity norm due to line delay and noise'
        counts, newBins, patches = plt.hist(imageCornerOffsetNorms, bins)
        print('  counts:\t{}\nbins:\t{}'.format(counts, newBins))        
        plt.title('Distribution of offset in the infinity norm due to line delay and noise')        
        plt.show()

        # TODO(jhuai): how about simulating IMU measurements? Do we need to account for IMU misalignment models?
        # How is the gravity accounted for?
        

    def naiveMethodToRsProjection(self, state_time, line_delay, verbose=False):
        """
        This method is not proved theoretically to converge, but it performs as precise as 
        Newton's method empirically, though slower.
        return: 
            1. Projected image corners according to a rolling shutter model, NX3X1 array.
        """
        imageCornerProjected= list()
        if verbose:
            print 'Naive method for state time %.9f' % state_time
        for iota in range(self.simulatedObs.getTotalTargetPoint()):
            # get the initial observation
            sm_T_w_c, isValid = getCameraPoseAt(state_time, self.simulPoseSplineDv, T_b_c=self.T_imu_c0)
            if not isValid:
                continue
            lastImagePoint = self.simulatedObs.projectATargetPoint(self.camGeometry, sm_T_w_c, iota) #3x1, generated with the GS model
            if lastImagePoint[2, 0] == 0.0:
                continue
            numIter = 0  
            aborted = False 
            if verbose:
                print 'lmId', iota, 'iter', numIter, 'image coords', lastImagePoint.T  
            if np.absolute(line_delay) < 1e-8:
                imageCornerProjected.append(lastImagePoint)
                continue         
            while numIter < 8:
                currTime = (lastImagePoint[1, 0] - self.imageHeight/2)*line_delay + state_time            
                sm_T_w_cx, isValid = getCameraPoseAt(currTime, self.simulPoseSplineDv, T_b_c=self.T_imu_c0)
                imagePoint = self.simulatedObs.projectATargetPoint(self.camGeometry, sm_T_w_cx, iota)
                if not isValid or imagePoint[2, 0] == 0.0:
                    aborted = True
                    break
                delta = np.absolute(lastImagePoint[1,0] - imagePoint[1,0])            
                numIter += 1
                if verbose:  
                    print 'lmId', iota, 'iter', numIter, 'image coords', imagePoint.T                    
                lastImagePoint = imagePoint
                if delta < 1e-3:
                    break
            if verbose:
                print
            if not aborted:
                imageCornerProjected.append(lastImagePoint)
        return np.array(imageCornerProjected)
    
    def newtonMethodToRsProjection(self, state_time, line_delay, reprojectionSigma = 1.0, time_offset=0.0, verbose = False):
        """
        return: 
            1. Projected image corners according to a rolling shutter model, NX3X1 array.
            2. Projected image corners plus gaussian noise of sigmas, a list of tuples, 
                each tuple (landmark index, keypoint index, pt.x, pt.y, keypoint size)
            3. The inf norm of the offset between the noisy measurement and projected 
                measurement according to a global shutter model.
        """
        imageCornerProjected = list() # image keypoints free of noise effect
        imageCornerProjectedOffset = list()
        frameKeypoints = list()
        kpId = 0
        if verbose:
            print 'Newton method for state time %.9f with time offset %.9f' % (state_time, time_offset)
        if state_time + time_offset <= self.simulPoseSplineDv.spline().t_min() or state_time + time_offset >= self.simulPoseSplineDv.spline().t_max():
            print "Warn: %.9f time out of range [%.9f, %.9f] in newton method Rs simulation" % (state_time + time_offset, self.simulPoseSplineDv.spline().t_min(), self.simulPoseSplineDv.spline().t_max())
            return np.array([[[]]]), frameKeypoints, list()
       
        for iota in range(self.simulatedObs.getTotalTargetPoint()): 
            sm_T_w_c, isValid = getCameraPoseAt(state_time + time_offset, self.simulPoseSplineDv, T_b_c=self.T_imu_c0)       
            lastImagePoint = self.simulatedObs.projectATargetPoint(self.camGeometry, sm_T_w_c, iota) #3x1, this is GS camera model
            if not isValid or lastImagePoint[2, 0] == 0.0:
                continue
            numIter = 0
            aborted = False
            if verbose:
                print 'lmId', iota, 'iter', numIter, 'image coords', lastImagePoint.T
            if np.absolute(line_delay) < 1e-8:
                imageCornerProjected.append(lastImagePoint)          
                xnoise = gauss(0.0, reprojectionSigma)
                ynoise = gauss(0.0, reprojectionSigma)
                imageCornerProjectedOffset.append(max(np.abs(xnoise), np.abs(ynoise)))
                frameKeypoints.append((iota, kpId, lastImagePoint[0, 0] + xnoise, lastImagePoint[1, 0] + ynoise, 12))
                kpId += 1        
                continue      
            # solve y=g(y) where y is the vertical projection in pixels
            initialImagePoint = lastImagePoint
            while numIter < 8:
                # now we have y_0, i.e., lastImagePoint[1, 0], complete the iteration by computing y_1

                # compute g(y_0)
                currTime = (lastImagePoint[1, 0] - self.imageHeight/2)*line_delay + state_time + time_offset
                sm_T_w_cx, isValid = getCameraPoseAt(currTime, self.simulPoseSplineDv, T_b_c=self.T_imu_c0)
                imagePoint0 = self.simulatedObs.projectATargetPoint(self.camGeometry, sm_T_w_cx, iota)
                if not isValid or imagePoint0[2, 0] == 0.0:
                    aborted = True
                    break
                   
                # compute Jacobian of g(y) relative to y at y_0
                eps = 1
                currTime = (lastImagePoint[1, 0] + eps - self.imageHeight/2)*line_delay + state_time + time_offset
                sm_T_w_cx, isValid = getCameraPoseAt(currTime, self.simulPoseSplineDv, T_b_c=self.T_imu_c0)
                imagePoint1 = self.simulatedObs.projectATargetPoint(self.camGeometry, sm_T_w_cx, iota)
                if not isValid or imagePoint1[2, 0] == 0.0:
                    aborted = True
                    break
                jacob = (imagePoint1[1, 0] - imagePoint0[1, 0])/eps
            
                # compute y_1
                lastImagePoint[0, 0] = imagePoint0[0, 0]    
                delta = imagePoint0[1, 0] - lastImagePoint[1, 0]        
                lastImagePoint[1, 0] = lastImagePoint[1, 0] - (imagePoint0[1, 0] - lastImagePoint[1, 0])/(jacob - 1)
                numIter += 1
                if verbose:  
                    print 'lmId', iota, 'iter', numIter, 'image coords', imagePoint0.T                    
                if np.absolute(delta) < 1e-4:
                    break

            if not aborted:
                # TODO: randomly mask some landmark observations
                imageCornerProjected.append(imagePoint0)                
                xnoise = gauss(0.0, reprojectionSigma)
                ynoise = gauss(0.0, reprojectionSigma)
                frameKeypoints.append((iota, kpId, imagePoint0[0, 0] + xnoise, imagePoint0[1, 0] + ynoise, 12))
                imageCornerProjectedOffset.append(max(np.abs(initialImagePoint[0, 0] - imagePoint0[0, 0] - xnoise), \
                  np.abs(initialImagePoint[1, 0] - imagePoint0[1, 0] - ynoise)))
                kpId += 1
        assert kpId == len(imageCornerProjected)
        return np.array(imageCornerProjected), frameKeypoints, imageCornerProjectedOffset
