import copy
import os
from random import gauss
import sys
import numpy as np
import matplotlib.pyplot as plt

import aslam_backend as aopt
import aslam_cv as acv
import aslam_cameras_april as acv_april
import sm
import kalibr_common as kc
import kalibr_errorterms as ket

import BSplineIO

def getCameraPoseAt(timeScalar, poseSplineDv, T_b_c):
    timeOffsetPadding = 0.0
    dv = aopt.Scalar(timeScalar)
    timeExpression = dv.toExpression()

    if timeScalar <= poseSplineDv.spline().t_min() or timeScalar >= poseSplineDv.spline().t_max():
        # print("getCameraPose: {:.9f} time out of range [{:.9f}, {:.9f}]".format( \
        #     timeScalar, poseSplineDv.spline().t_min(), poseSplineDv.spline().t_max()))
        return sm.Transformation(), False

    T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
    sm_T_w_c = sm.Transformation(T_w_b.toTransformationMatrix())*T_b_c  
    return sm_T_w_c, True


def printExtraCameraDetails(camConfig):
    resolution = camConfig.getResolution()
    print('  Camera resolution: {}'.format(resolution))
    imageNoise = camConfig.getImageNoise()
    print('  Image noise std dev: {}'.format(imageNoise))
    lineDelay = camConfig.getLineDelayNanos()
    print("  Line delay: {} ns".format(lineDelay))
    updateRate = camConfig.getUpdateRate()
    print("  Update rate: {} Hz".format(updateRate))


def printExtraImuDetails(imuConfig):
    initialGyroBias = imuConfig.getInitialGyroBias()
    print('  Initial gyro bias: {}'.format(initialGyroBias))
    initialAccBias = imuConfig.getInitialAccBias()
    print('  Initial accelerometer bias: {}'.format(initialAccBias))
    gravityInTarget = imuConfig.getGravityInTarget()
    print('  Gravity in target: {}'.format(gravityInTarget))


class RsVisualInertialMeasViaBSplineSimulator(object):
    '''
    simulate visual(rolling shutter) inertial measurements with provided
    BSpline models representing realistic motion and IMU biases.
    '''
    def __init__(self, args):
        self.pose_file = args.pose_file
        self.poseSplineDv = BSplineIO.loadPoseBSpline(args.pose_file)
        self.gyroBiasSplineDv = BSplineIO.loadBSpline(args.gyro_bias_file)
        self.accBiasSplineDv = BSplineIO.loadBSpline(args.acc_bias_file)
        self.showOnScreen = not args.dontShowReport

        print("Camera chain from {}".format(args.chain_yaml))
        chain = kc.CameraChainParameters(args.chain_yaml)        
        camNr = 0
        self.T_c0_imu = chain.getExtrinsicsImuToCam(camNr)  
        self.T_imu_c0 = self.T_c0_imu.inverse()
        camConfig = chain.getCameraParameters(camNr)
        camConfig.printDetails()
        printExtraCameraDetails(camConfig)

        camera = kc.AslamCamera.fromParameters(camConfig)
        self.camGeometry = camera.geometry
        self.cameraConfig = camConfig

        targetConfig = kc.CalibrationTargetParameters(args.target_yaml)
        print("Target used in the simulation:")
        targetConfig.printDetails()
        self.targetObservation = None
        self.allTargetCorners = None
        self.setupCalibrationTarget(targetConfig, showExtraction=False, showReproj=False, imageStepping=False)
        self.imageHeight = self.cameraConfig.getResolution()[1]
        self.timeOffset = chain.getTimeshiftCamImu(camNr)

        print("IMU configuration:")
        self.imuConfig = kc.ImuParameters(args.imu_yaml)
        self.imuConfig.printDetails()
        printExtraImuDetails(self.imuConfig)

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

        self.targetObservation = acv.GridCalibrationTargetObservation(grid)
        self.allTargetCorners = self.targetObservation.getAllCornersTargetFrame() # nx3
        assert self.allTargetCorners.shape[0] == self.targetObservation.getTotalTargetPoint()

    def checkNaiveVsNewtonRsProjection(self):
        for frameId in range(1):
            state_time = self.refStateTimes[0]
            line_delay = float(self.cameraConfig.getLineDelayNanos()) * 1e-9
            imageCornersNaive = self.naiveMethodToRsProjection(state_time, line_delay, False)
            imageCornersNewton, unusedKeypoints, _ = \
                self.newtonMethodToRsProjection(state_time, line_delay, 1.0, False)

        assert np.allclose(imageCornersNaive[:, :, 0], imageCornersNewton[:, :, 0])
        reproducedImageCornerFile = self.pose_file.replace("pose", "naive_vs_newton", 1)
        np.savetxt(reproducedImageCornerFile,
                   np.concatenate((imageCornersNaive[:, :, 0], imageCornersNewton[:, :, 0]), axis=1), \
                   fmt=['%.9f', '%.9f', '%d', '%.9f', '%.9f', '%d'])

    def __generateSampleTimes(self, tmin, tmax, rate):
        timeList = list()
        interval = 1.0 / rate
        t = tmin + 1e-6
        while t < tmax:
            timeList.append(t)
            t += interval
        return timeList

    def __generateStateTimes(self, rate, timePadding):
        tmin = max(self.poseSplineDv.spline().t_min(), self.gyroBiasSplineDv.spline().t_min(), \
                self.accBiasSplineDv.spline().t_min()) + timePadding
        tmax = min(self.poseSplineDv.spline().t_max(), self.gyroBiasSplineDv.spline().t_max(), \
                self.accBiasSplineDv.spline().t_max()) - timePadding
        return self.__generateSampleTimes(tmin, tmax, rate)

    def simulateImuData(self, trueImuTimes):
        """simulate inertial measurements at true epochs without time offset.
        Imu biases are added. White noise, and random walk are also optional.
        """
        q_i_b_prior = np.array([0., 0., 0., 1.])
        q_i_b_Dv = aopt.RotationQuaternionDv(q_i_b_prior)
        r_b_Dv = aopt.EuclideanPointDv(np.array([0., 0., 0.]))

        # gravity in target example: np.array([0.0, 9.81, 0.0])
        gravity = self.imuConfig.getGravityInTarget()
        gravityDv = aopt.EuclideanDirection(np.array(self.imuConfig.getGravityInTarget()).T)
        gravityExpression = gravityDv.toExpression()

        omegaDummy = np.zeros((3, 1))
        alphaDummy = np.zeros((3, 1))
        weightDummy = 1.0

        imuData = np.zeros((len(trueImuTimes), 6))
        imuBiases = np.zeros((len(trueImuTimes), 6))

        gyroSpline = self.gyroBiasSplineDv.spline()
        accSpline = self.accBiasSplineDv.spline()

        gyroNoiseDiscrete, gyroNoise, gyroWalk = self.imuConfig.getGyroStatistics()
        accNoiseDiscrete, accNoise, accWalk = self.imuConfig.getAccelerometerStatistics()
        Rgyro = np.eye(3) * gyroNoiseDiscrete * gyroNoiseDiscrete
        Raccel = np.eye(3) * accNoiseDiscrete * accNoiseDiscrete
        omegaInvR = np.linalg.inv(Rgyro)
        alphaInvR = np.linalg.inv(Raccel)

        # TODO(jhuai): add IMU noise.
        for index, tk in enumerate(trueImuTimes):
            # GyroscopeError(measurement, invR, angularVelocity, bias)
            w_b = self.poseSplineDv.angularVelocityBodyFrame(tk)
            b_i = self.gyroBiasSplineDv.toEuclideanExpression(tk,0)
            C_i_b = q_i_b_Dv.toExpression()
            w = C_i_b * w_b
            gerr = ket.EuclideanError(omegaDummy, omegaInvR * weightDummy, w + b_i)
            omega = gerr.getPredictedMeasurement()
            gyroBias = gyroSpline.eval(tk)

            # check
            gerr2 = ket.EuclideanError(omegaDummy, omegaInvR * weightDummy, w)
            omega2 = gerr2.getPredictedMeasurement()
            assert np.linalg.norm(omega - omega2 - gyroBias) < 1e-8

            C_b_w = self.poseSplineDv.orientation(tk).inverse()
            a_w = self.poseSplineDv.linearAcceleration(tk)
            b_i = self.accBiasSplineDv.toEuclideanExpression(tk,0)
            w_b = self.poseSplineDv.angularVelocityBodyFrame(tk)
            w_dot_b = self.poseSplineDv.angularAccelerationBodyFrame(tk)
            C_i_b = q_i_b_Dv.toExpression()
            r_b = r_b_Dv.toExpression()
            a = C_i_b * (C_b_w * (a_w - gravityExpression) + \
                            w_dot_b.cross(r_b) + w_b.cross(w_b.cross(r_b)))
            aerr = ket.EuclideanError(alphaDummy, alphaInvR * weightDummy, a + b_i)
            alpha = aerr.getPredictedMeasurement()
            accBias = accSpline.eval(tk)

            imuData[index, :3] = omega
            imuData[index, 3:] = alpha
            imuBiases[index, :3] = gyroBias
            imuBiases[index, 3:] = accBias
        return trueImuTimes, imuData, imuBiases

    def simulateCameraObservations(self, trueFrameTimes, outputDir):
        '''simulate camera observations for frames at all ref state times and plus noise'''
        # simulate camera observations, save to vertices, tracks, observations, landmarks per maplab csv format.
        # https://github.com/ethz-asl/maplab/wiki/CSV-Dataset-Format
        # Descriptors are not needed. In tracks, track_id can be set to -1 as it is not used for now.
        # Timestamps of vertices and tracks should be in camera clock. 
        # The timestamp in tracks.csv for each keypoint is the timestamp for the observing frame.

        # simulate RS observations at state times, but the camera timestamps are shifted by offset.
        imageCornerOffsetNorms = list()
        bins = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, \
                3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0]
        imageNoise = self.cameraConfig.getImageNoise()
        frameKeypointList = list()
        landmark_observations = dict()
        for iota in range(self.targetObservation.getTotalTargetPoint()):
            landmark_observations[iota]=list()
        cameraIndex = 0 # assume only one camera is used.
        cameraTimeOffset = self.timeOffset
        for vertexId, frameTime in enumerate(trueFrameTimes):
            _, noisyKeypoints, keypointOffsets = \
                self.newtonMethodToRsProjection(frameTime,
                                                float(self.cameraConfig.getLineDelayNanos()) * 1e-9,
                                                imageNoise)
            imageCornerOffsetNorms += keypointOffsets
            if vertexId % 300 == 0:
                print('  Projected {:d} target landmarks for state at {:.9f}'.format(len(noisyKeypoints), frameTime))
            for keypoint in noisyKeypoints:
                landmark_observations[keypoint[0]].append(
                    (vertexId, keypoint[1], keypoint[2], keypoint[3], keypoint[4]))
            frameKeypointList.append(noisyKeypoints)

        observationCsv = os.path.join(outputDir, "observations.csv")
        with open(observationCsv, 'w') as stream:
            header = ', '.join(["vertex index", "frame index", "keypoint index", "landmark index"])
            stream.write('{}\n'.format(header))
            probe = 0
            for landmarkId, observationList in sorted(landmark_observations.iteritems()):
                assert probe == landmarkId
                for observation in observationList:
                    stream.write('{}, {}, {}, {}\n'.format(observation[0], cameraIndex, observation[1], landmarkId))
                probe += 1

        trackCsv = os.path.join(outputDir, "tracks.csv")
        with open(trackCsv, 'w') as stream:
            header = ', '.join(
                ["timestamp [ns]", "vertex index", "frame index", "keypoint index", "keypoint measurement 0 [px]",
                 "keypoint measurement 1 [px]", "keypoint measurement uncertainty", "keypoint scale",
                 "keypoint track id"])
            stream.write('{}\n'.format(header))
            for vertexId, frameKeypoints in enumerate(frameKeypointList):
                for keypoint in frameKeypoints:
                    timeString = BSplineIO.toNanosecondString(trueFrameTimes[vertexId] - cameraTimeOffset)
                    stream.write('{}, {}, {}, {:d}, {:.5f}, {:.5f}, {}, {}, {}\n'.format(
                        timeString, vertexId, cameraIndex, keypoint[1], keypoint[2], keypoint[3], imageNoise,
                        keypoint[4], -1))

        print('  Written landmark observations to {}'.format(observationCsv))
        print('  Histogram of norm of the offset due to line delay and noise')
        counts, newBins, patches = plt.hist(imageCornerOffsetNorms, bins)
        print('  counts:{}\n  bins:{}'.format(counts, newBins))
        plt.title('Distribution of norm of the offsets due to line delay and noise')
        if self.showOnScreen:
            plt.show()

    def simulate(self, outputDir):
        landmarkCsv = os.path.join(outputDir, "landmarks.csv")
        print("Saving landmarks to {}...".format(landmarkCsv))
        with open(landmarkCsv, 'w') as stream:
            header = ', '.join(["landmark index", "landmark position x [m]",
                                "landmark position y [m]", "landmark position z [m]"])
            stream.write('{}\n'.format(header))
            for index, row in enumerate(self.allTargetCorners):
                stream.write("{}, {}, {}, {}\n".format(index, row[0], row[1], row[2]))
        timePadding = 2.0 / self.cameraConfig.getUpdateRate()
        trueFrameTimes = self.__generateStateTimes(self.cameraConfig.getUpdateRate(), timePadding)

        print('Simulating states...')
        print("  Camera frame true start time {:.9f} and true finish time {:.9f}".format( \
                trueFrameTimes[0], trueFrameTimes[-1]))
        vertexCsv = os.path.join(outputDir, "vertices.csv")
        with open(vertexCsv, 'w') as vertexStream:
            BSplineIO.saveStates(trueFrameTimes, self.poseSplineDv, self.gyroBiasSplineDv.spline(), \
                    self.accBiasSplineDv.spline(), self.timeOffset, vertexStream)
            print("  Written simulated states to {}".format(vertexCsv))

        print("Simulating IMU data...")
        imuTimePadding = 2.0 / self.cameraConfig.getUpdateRate()
        trueImuTimes = self.__generateStateTimes(self.imuConfig.getUpdateRate(), imuTimePadding)
        imuTimes, imuData, imuBiases = self.simulateImuData(trueImuTimes)
        imuCsv = os.path.join(outputDir, "imu.csv")
        with open(imuCsv, "w") as stream:
            header = ', '.join(["timestamp [ns]", "acc x [m/s^2]", "acc y [m/s^2]", "acc z [m/s^2]",
                                   "gyro x [rad/s]", "gyro y [rad/s]", "gyro z [rad/s]", "bias acc x [m/s^2]",
                                   "bias acc y [m/s^2]", "bias acc z [m/s^2]", "bias gyro x [rad/s]",
                                   "bias gyro y [rad/s]", "bias gyro z [rad/s]"])
            stream.write('{}\n'.format(header))
            for index, time in enumerate(imuTimes):
                dataString = ', '.join(map(str, imuData[index, :]))
                biasString = ', '.join(map(str, imuBiases[index, :]))
                stream.write("{}, {}, {}\n".format(BSplineIO.toNanosecondString(time), dataString, biasString))

        print("Simulating camera observations...")
        self.simulateCameraObservations(trueFrameTimes, outputDir)

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
        for iota in range(self.targetObservation.getTotalTargetPoint()):
            # get the initial observation
            sm_T_w_c, isValid = getCameraPoseAt(state_time, self.poseSplineDv, self.T_imu_c0)
            if not isValid:
                continue
            lastImagePoint = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_c, iota) # 3x1.
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
                sm_T_w_cx, isValid = getCameraPoseAt(currTime, self.poseSplineDv, self.T_imu_c0)
                imagePoint = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_cx, iota)
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

    def newtonMethodToRsProjection(self, state_time, line_delay, reprojectionSigma = 1.0, verbose = False):
        """
        params:
            state_time: camera mid exposure timestamp without time offset or rolling shutter effect.
            For landmark i observed at vertical coordinate v_i in frame j, the state_time, t_j_imu
            satisfies t_j_imu = t_j_cam + t_offset.
            With the rolling shutter, we have t_j_cam + t_offset + (v_i - 0.5 * h) * t_line = t_j_i.

        return:
            1. Projected landmarks in image according to a rolling shutter model, NX3X1 array.
            2. Projected landmarks in image plus gaussian noise, a list of tuples,
                each tuple (landmark index, keypoint index, pt.x, pt.y, keypoint size)
            3. The norm of the offset between the noisy measurement and projected
                measurement according to a global shutter model.
        """
        imageCornerProjected = list() # image keypoints free of noise effect
        imageCornerProjectedOffset = list()
        frameKeypoints = list()
        kpId = 0
        if verbose:
            print('Newton method for state time {:.9f}'.format(state_time))
        if state_time <= self.poseSplineDv.spline().t_min() or state_time >= self.poseSplineDv.spline().t_max():
            print("RS projection warn: {:.9f} time out of range [{:.9f}, {:.9f}] in newton method Rs simulation".
                format(state_time, self.poseSplineDv.spline().t_min(), self.poseSplineDv.spline().t_max()))
            return np.array([[[]]]), frameKeypoints, list()

        numOutOfBound = 0
        numFailedProjection = 0
        numLandmarks = self.targetObservation.getTotalTargetPoint()
        for iota in range(numLandmarks):
            sm_T_w_c, isValid = getCameraPoseAt(state_time, self.poseSplineDv, self.T_imu_c0)

            lastImagePoint = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_c, iota) # 3x1.
            if not isValid:
                numOutOfBound += 1
                continue
            if lastImagePoint[2, 0] == 0.0:
                numFailedProjection += 1
                continue
            numIter = 0
            aborted = False
            if verbose:
                print('lmId {} iter {} image coords {}'.format(iota, numIter, lastImagePoint.T))
            if np.absolute(line_delay) < 1e-8:
                imageCornerProjected.append(lastImagePoint)          
                xnoise = gauss(0.0, reprojectionSigma)
                ynoise = gauss(0.0, reprojectionSigma)
                imageCornerProjectedOffset.append(np.linalg.norm([xnoise, ynoise]))
                frameKeypoints.append((iota, kpId, lastImagePoint[0, 0] + xnoise, lastImagePoint[1, 0] + ynoise, 12))
                kpId += 1
                continue
            # solve y=g(y) where y is the vertical projection in pixels
            initialImagePoint = copy.deepcopy(lastImagePoint)
            while numIter < 6:
                # now we have y_0, i.e., lastImagePoint[1, 0], complete the iteration by computing y_1

                # compute g(y_0)
                currTime = (lastImagePoint[1, 0] - self.imageHeight/2)*line_delay + state_time
                sm_T_w_cx, isValid = getCameraPoseAt(currTime, self.poseSplineDv, self.T_imu_c0)

                imagePoint0 = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_cx, iota)
                if not isValid:
                    numOutOfBound += 1
                    aborted = True
                    break
                if imagePoint0[2, 0] == 0.0:
                    numFailedProjection += 1
                    aborted = True
                    break
                # compute Jacobian of g(y) relative to y at y_0
                eps = 1
                currTime = (lastImagePoint[1, 0] + eps - self.imageHeight/2)*line_delay + state_time
                sm_T_w_cx, isValid = getCameraPoseAt(currTime, self.poseSplineDv, self.T_imu_c0)

                imagePoint1 = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_cx, iota)
                if not isValid:
                    numOutOfBound += 1
                    aborted = True
                    break
                if imagePoint0[2, 0] == 0.0:
                    numFailedProjection += 1
                    aborted = True
                    break
                jacob = (imagePoint1[1, 0] - imagePoint0[1, 0])/eps

                # compute y_1
                lastImagePoint[0, 0] = imagePoint0[0, 0]    
                delta = imagePoint0[1, 0] - lastImagePoint[1, 0]        
                lastImagePoint[1, 0] = lastImagePoint[1, 0] - (imagePoint0[1, 0] - lastImagePoint[1, 0])/(jacob - 1)
                numIter += 1
                if verbose:
                    print('lmId {} iter {} image coords {}'.format(iota, numIter, imagePoint0.T))
                if np.absolute(delta) < 1e-4:
                    break

            if not aborted:
                # TODO: randomly mask some landmark observations
                imageCornerProjected.append(imagePoint0)
                xnoise = gauss(0.0, reprojectionSigma)
                ynoise = gauss(0.0, reprojectionSigma)
                frameKeypoints.append((iota, kpId, imagePoint0[0, 0] + xnoise, imagePoint0[1, 0] + ynoise, 12))
                imageCornerProjectedOffset.append(np.linalg.norm([initialImagePoint[0, 0] - imagePoint0[0, 0] - xnoise,
                                                                  initialImagePoint[1, 0] - imagePoint0[1, 0] - ynoise]))
                kpId += 1
        if numOutOfBound > 0 or numFailedProjection > numLandmarks / 2:
            print("  For frame at {:.6f} s, {} out of time bound landmarks, {} failed to project landmarks".format( \
                state_time, numOutOfBound, numFailedProjection))

        assert kpId == len(imageCornerProjected)
        return np.array(imageCornerProjected), frameKeypoints, imageCornerProjectedOffset
