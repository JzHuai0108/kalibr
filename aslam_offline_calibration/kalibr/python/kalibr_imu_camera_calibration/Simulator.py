import copy
import math
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


def addNoiseToImuReadings(imuMeasurements, imuParameters):
    """
    :param imuParameters: imuConfig
    :param times: time of each IMU reading in seconds.
    :param imuMeasurements: numpy array, N x 6, accel and gyro data.
    :return: 
    """
    trueBiases = copy.deepcopy(imuMeasurements)
    noisyImuMeasurements = copy.deepcopy(imuMeasurements)
    bgk = imuParameters.getInitialGyroBias()
    bak = imuParameters.getInitialAccBias()
    gyroNoiseDiscrete, gyroWalk, gyroNoise = imuParameters.getGyroStatistics()
    accNoiseDiscrete, accWalk, accNoise = imuParameters.getAccelerometerStatistics()
    sqrtRate = math.sqrt(imuParameters.getUpdateRate())
    sqrtDeltaT = 1.0 / sqrtRate

    for index, reading in enumerate(imuMeasurements):
        trueBiases[index, :3] = bak
        trueBiases[index, 3:] = bgk
        # eq 50, Oliver Woodman, An introduction to inertial navigation
        noisyImuMeasurements[index, :3] = imuMeasurements[index, :3] + bak + np.random.normal(0, accNoiseDiscrete, 3)
        noisyImuMeasurements[index, 3:] = imuMeasurements[index, 3:] + bgk + np.random.normal(0, gyroNoiseDiscrete, 3)
        # eq 51, Oliver Woodman, An introduction to inertial navigation,
        # we do not divide sqrtDeltaT by sqrtT because sigma_gw_c is bias white noise density
        # for bias random walk (BRW) whereas eq 51 uses bias instability (BS) having the
        # same unit as the IMU measurements. also see eq 9 therein.
        bgk += np.random.normal(0, gyroWalk * sqrtDeltaT, 3)
        bak += np.random.normal(0, accWalk * sqrtDeltaT, 3)
    return noisyImuMeasurements, trueBiases


class RsCameraSimulator(object):
    def __init__(self, args):
        self.pose_file = args.pose_file
        self.poseSplineDv = BSplineIO.selectiveLoadPoseBSpline(args.pose_file)
        self.showOnScreen = not args.dontShowReport

        print("Camera chain from {}".format(args.chain_yaml))
        chain = kc.CameraChainParameters(args.chain_yaml)
        camNr = 0
        camConfig = chain.getCameraParameters(camNr)
        camConfig.printDetails()
        printExtraCameraDetails(camConfig)
        self.T_imu_cam = sm.Transformation()
        self.timeOffset = 0

        camera = kc.AslamCamera.fromParameters(camConfig)
        self.camGeometry = camera.geometry
        self.cameraConfig = camConfig

        targetConfig = kc.CalibrationTargetParameters(args.target_yaml)
        print("Target used in the simulation:")
        targetConfig.printDetails()
        self.targetObservation = None
        self.allTargetCorners = None
        self.setupCalibrationTarget(targetConfig, showExtraction=False, showReproj=False, imageStepping=False)
        self.imageWidth = self.cameraConfig.getResolution()[0]
        self.imageHeight = self.cameraConfig.getResolution()[1]

    def setupCalibrationTarget(self, targetConfig, showExtraction=False, showReproj=False, imageStepping=False):
        '''copied from IccCamera class'''
        # load the calibration target configuration
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
            options.minTagsForValidObs = int(np.max([targetParams['tagRows'], targetParams['tagCols']]) + 1)

            grid = acv_april.GridCalibrationTargetAprilgrid(targetParams['tagRows'],
                                                            targetParams['tagCols'],
                                                            targetParams['tagSize'],
                                                            targetParams['tagSpacing'],
                                                            options)
        else:
            raise RuntimeError("Unknown calibration target.")

        options = acv.GridDetectorOptions()
        options.imageStepping = imageStepping
        options.plotCornerReprojection = showReproj
        options.filterCornerOutliers = True

        self.targetObservation = acv.GridCalibrationTargetObservation(grid)
        self.allTargetCorners = self.targetObservation.getAllCornersTargetFrame()  # nx3
        assert self.allTargetCorners.shape[0] == self.targetObservation.getTotalTargetPoint()

    def generateSampleTimes(self, tmin, tmax, rate):
        timeList = list()
        interval = 1.0 / rate
        t = tmin
        while t < tmax:
            timeList.append(t)
            t += interval
        return timeList

    def generateStateTimes(self, rate, timePadding):
        tmin = self.poseSplineDv.spline().t_min() + timePadding
        tmax = self.poseSplineDv.spline().t_max() - timePadding
        return self.generateSampleTimes(tmin, tmax, rate)

    def checkNaiveVsNewtonRsProjection(self, outputDir):
        timePadding = 2.5 / self.cameraConfig.getUpdateRate()
        trueFrameTimes = self.generateStateTimes(self.cameraConfig.getUpdateRate(), timePadding)
        state_time = trueFrameTimes[0]
        line_delay = float(self.cameraConfig.getLineDelayNanos()) * 1e-9
        imageCornersNaive = self.naiveMethodToRsProjection(state_time, line_delay, self.T_imu_cam, False)
        imageCornersNewton, unusedKeypoints, _ = \
            self.newtonMethodToRsProjection(state_time, line_delay, self.T_imu_cam, 1.0, False)
        if imageCornersNaive.shape[0] == 0 or imageCornersNewton.shape[0] == 0:
            print("None successfully projected landmarks!")
        elif imageCornersNaive.shape[0] == imageCornersNewton.shape[0]:
            print("naive shape {} newton shape {}".format(imageCornersNaive.shape, imageCornersNewton.shape))
            assert np.allclose(imageCornersNaive[:, :, 0], imageCornersNewton[:, :, 0])
            imageCornerFile = os.path.join(outputDir, "naive_vs_newton_rs_check.txt")
            np.savetxt(imageCornerFile, np.concatenate((imageCornersNaive[:, :, 0], imageCornersNewton[:, :, 0]), axis=1),
                       fmt=['%.9f', '%.9f', '%d', '%.9f', '%.9f', '%d'])
        else:
            print('#Naive projections {} #Newton projections {}'.format(imageCornersNaive.shape[0],
                                                                        imageCornersNewton.shape[0]))
            lastIndex = min(imageCornersNaive.shape[0], imageCornersNewton.shape[0], 20)
            if lastIndex > 0:
                print('First {} rows of naive and newton projections:\n{}'.format(
                    lastIndex,
                    np.concatenate((imageCornersNaive[:lastIndex, :, 0], imageCornersNewton[:lastIndex, :, 0]),
                                   axis=1)))

    def naiveMethodToRsProjection(self, state_time, line_delay, T_imu_cam, verbose=False):
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
            sm_T_w_c, validPose = getCameraPoseAt(state_time, self.poseSplineDv, T_imu_cam)
            if not validPose:
                continue
            validProjection, lastImagePoint = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_c, iota) # 3x1.
            if not validProjection:
                continue
            numIter = 0
            aborted = False
            if verbose:
                print 'lmId', iota, 'iter', numIter, 'image coords', lastImagePoint.T
            if np.absolute(line_delay) < 1e-8:
                imageCornerProjected.append(lastImagePoint)
                continue
            while numIter < 8:
                currTime = (lastImagePoint[1, 0] - self.imageHeight * 0.5) * line_delay + state_time
                sm_T_w_cx, validPose = getCameraPoseAt(currTime, self.poseSplineDv, T_imu_cam)
                validProjection, imagePoint = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_cx, iota)
                if not validPose or not validProjection:
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

    def newtonMethodToRsProjection(self, state_time, line_delay, T_imu_cam, reprojectionSigma = 1.0, verbose = False):
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
            sm_T_w_c, validPose = getCameraPoseAt(state_time, self.poseSplineDv, T_imu_cam)
            validProjection, lastImagePoint = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_c, iota) # 3x1.
            if not validPose:
                numOutOfBound += 1
                continue
            if not validProjection:
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
                currTime = (lastImagePoint[1, 0] - self.imageHeight * 0.5) * line_delay + state_time
                sm_T_w_cx, validPose = getCameraPoseAt(currTime, self.poseSplineDv, T_imu_cam)

                validProjection, imagePoint0 = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_cx, iota)
                if not validPose:
                    numOutOfBound += 1
                    aborted = True
                    break
                if not validProjection:
                    numFailedProjection += 1
                    aborted = True
                    break
                # compute Jacobian of g(y) relative to y at y_0
                eps = 1
                currTime = (lastImagePoint[1, 0] + eps - self.imageHeight * 0.5) * line_delay + state_time
                sm_T_w_cx, validPose = getCameraPoseAt(currTime, self.poseSplineDv, T_imu_cam)

                validProjection, imagePoint1 = self.targetObservation.projectATargetPoint(self.camGeometry, sm_T_w_cx, iota)
                if not validPose:
                    numOutOfBound += 1
                    aborted = True
                    break
                if not validProjection:
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
                imageCornerProjected.append(imagePoint0)
                xnoise = gauss(0.0, reprojectionSigma)
                ynoise = gauss(0.0, reprojectionSigma)
                noisyPoint = [noisyValue(imagePoint0[0, 0], self.imageWidth, xnoise),
                              noisyValue(imagePoint0[1, 0], self.imageHeight, ynoise)]
                frameKeypoints.append((iota, kpId, noisyPoint[0], noisyPoint[1], 12))
                imageCornerProjectedOffset.append(np.linalg.norm([initialImagePoint[0, 0] - noisyPoint[0],
                                                                  initialImagePoint[1, 0] - noisyPoint[1]]))
                kpId += 1
        if numOutOfBound > 0 or numFailedProjection > numLandmarks / 2:
            print("  For frame at {:.6f} s, {} out of time bound landmarks, {} failed to project landmarks".format( \
                state_time, numOutOfBound, numFailedProjection))

        assert kpId == len(imageCornerProjected)
        return np.array(imageCornerProjected), frameKeypoints, imageCornerProjectedOffset

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
                                                self.T_imu_cam,
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
                    timeString = BSplineIO.secondToNanosecondString(trueFrameTimes[vertexId] - cameraTimeOffset)
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

    def simulateLandmarks(self, outputDir):
        landmarkCsv = os.path.join(outputDir, "landmarks.csv")
        print("Saving landmarks to {}...".format(landmarkCsv))
        with open(landmarkCsv, 'w') as stream:
            header = ', '.join(["landmark index", "landmark position x [m]",
                                "landmark position y [m]", "landmark position z [m]"])
            stream.write('{}\n'.format(header))
            for index, row in enumerate(self.allTargetCorners):
                stream.write("{}, {}, {}, {}\n".format(index, row[0], row[1], row[2]))

    def computeCameraRate(self):
        lineDelay = self.cameraConfig.getLineDelayNanos()
        maxFrameRate = math.floor(1e9 / ((lineDelay + 1000) * self.imageHeight))
        cameraRate = min(maxFrameRate, self.cameraConfig.getUpdateRate())
        return cameraRate

    def simulateStates(self, outputDir):
        cameraRate = self.computeCameraRate()
        timePadding = 2.5 / cameraRate
        trueFrameTimes = self.generateStateTimes(cameraRate, timePadding)

        print('Simulating states at camera rate {}...'.format(cameraRate))
        print("  Camera frame true start time {:.9f} and true finish time {:.9f}".format(
            trueFrameTimes[0], trueFrameTimes[-1]))
        vertexCsv = os.path.join(outputDir, "vertices.csv")
        with open(vertexCsv, 'w') as vertexStream:
            BSplineIO.saveCameraStates(trueFrameTimes, self.poseSplineDv, vertexStream)
            print("  Written simulated states to {}".format(vertexCsv))
        return trueFrameTimes

    def simulate(self, outputDir):
        self.simulateLandmarks(outputDir)
        trueFrameTimes = self.simulateStates(outputDir)
        print("Simulating camera observations...")
        self.simulateCameraObservations(trueFrameTimes, outputDir)

class RsCameraImuSimulator(RsCameraSimulator):
    '''
    simulate visual(rolling shutter) inertial measurements with provided
    BSpline models representing realistic motion and IMU biases.
    '''
    def __init__(self, args):
        super(RsCameraImuSimulator, self).__init__(args)
        self.biasFromSplines = args.biasFromSplines
        self.gyroBiasSplineDv = None
        self.accBiasSplineDv = None
        if self.biasFromSplines:
            self.gyroBiasSplineDv = BSplineIO.loadBSpline(args.gyro_bias_file)
            self.accBiasSplineDv = BSplineIO.loadBSpline(args.acc_bias_file)

        chain = kc.CameraChainParameters(args.chain_yaml)
        camNr = 0
        T_cam_imu = chain.getExtrinsicsImuToCam(camNr)
        self.T_imu_cam = T_cam_imu.inverse()
        self.timeOffset = chain.getTimeshiftCamImu(camNr)

        print("IMU configuration:")
        self.imuConfig = kc.ImuParameters(args.imu_yaml)
        self.imuConfig.printDetails()
        printExtraImuDetails(self.imuConfig)


    def generateStateTimes(self, rate, timePadding):
        if self.gyroBiasSplineDv:
            tmin = max(self.poseSplineDv.spline().t_min(), self.gyroBiasSplineDv.spline().t_min(),
                    self.accBiasSplineDv.spline().t_min()) + timePadding
            tmax = min(self.poseSplineDv.spline().t_max(), self.gyroBiasSplineDv.spline().t_max(),
                    self.accBiasSplineDv.spline().t_max()) - timePadding
        else:
            tmin = self.poseSplineDv.spline().t_min() + timePadding
            tmax = self.poseSplineDv.spline().t_max() - timePadding
        return self.generateSampleTimes(tmin, tmax, rate)

    def simulateImuDataAtTimes(self, trueImuTimes):
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

        gyroNoiseDiscrete, gyroNoise, gyroWalk = self.imuConfig.getGyroStatistics()
        accNoiseDiscrete, accNoise, accWalk = self.imuConfig.getAccelerometerStatistics()
        Rgyro = np.eye(3) * gyroNoiseDiscrete * gyroNoiseDiscrete
        Raccel = np.eye(3) * accNoiseDiscrete * accNoiseDiscrete
        omegaInvR = np.linalg.inv(Rgyro)
        alphaInvR = np.linalg.inv(Raccel)

        for index, tk in enumerate(trueImuTimes):
            # GyroscopeError(measurement, invR, angularVelocity, bias)
            w_b = self.poseSplineDv.angularVelocityBodyFrame(tk)
            C_i_b = q_i_b_Dv.toExpression()
            w = C_i_b * w_b
            if self.biasFromSplines:
                b_i = self.gyroBiasSplineDv.toEuclideanExpression(tk, 0)
                gerr = ket.EuclideanError(omegaDummy, omegaInvR * weightDummy, w + b_i)
                omega = gerr.getPredictedMeasurement()
                gyroSpline = self.gyroBiasSplineDv.spline()
                gyroBias = gyroSpline.eval(tk)
            else:
                gerr = ket.EuclideanError(omegaDummy, omegaInvR * weightDummy, w)
                omega = gerr.getPredictedMeasurement()
                gyroBias = np.array([0, 0, 0])

            C_b_w = self.poseSplineDv.orientation(tk).inverse()
            a_w = self.poseSplineDv.linearAcceleration(tk)
            w_b = self.poseSplineDv.angularVelocityBodyFrame(tk)
            w_dot_b = self.poseSplineDv.angularAccelerationBodyFrame(tk)
            C_i_b = q_i_b_Dv.toExpression()
            r_b = r_b_Dv.toExpression()
            a = C_i_b * (C_b_w * (a_w - gravityExpression) + \
                            w_dot_b.cross(r_b) + w_b.cross(w_b.cross(r_b)))
            if self.biasFromSplines:
                b_i = self.accBiasSplineDv.toEuclideanExpression(tk, 0)
                aerr = ket.EuclideanError(alphaDummy, alphaInvR * weightDummy, a + b_i)
                alpha = aerr.getPredictedMeasurement()
                accSpline = self.accBiasSplineDv.spline()
                accBias = accSpline.eval(tk)
            else:
                aerr = ket.EuclideanError(alphaDummy, alphaInvR * weightDummy, a)
                alpha = aerr.getPredictedMeasurement()
                accBias = np.array([0, 0, 0])

            imuData[index, :3] = alpha
            imuData[index, 3:] = omega
            imuBiases[index, :3] = accBias
            imuBiases[index, 3:] = gyroBias

        if not self.biasFromSplines:
            print("  Add noise to IMU readings...")
            noisyImuData, imuBiases = addNoiseToImuReadings(imuData, self.imuConfig)
        else:
            noisyImuData = imuData

        return trueImuTimes, noisyImuData, imuBiases

    def simulateImuData(self, outputDir):
        print("Simulating IMU data...")
        cameraRate = self.computeCameraRate()
        imuTimePadding = 2.0 / cameraRate
        trueImuTimes = self.generateStateTimes(self.imuConfig.getUpdateRate(), imuTimePadding)
        imuTimes, imuData, imuBiases = self.simulateImuDataAtTimes(trueImuTimes)
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
                stream.write("{}, {}, {}\n".format(BSplineIO.secondToNanosecondString(time), dataString, biasString))

    def simulateStates(self, outputDir):
        cameraRate = self.computeCameraRate()
        timePadding = 2.5 / cameraRate
        trueFrameTimes = self.generateStateTimes(cameraRate, timePadding)

        print('Simulating states at camera rate {}...'.format(cameraRate))
        print("  Camera frame true start time {:.9f} and true finish time {:.9f}".format(
            trueFrameTimes[0], trueFrameTimes[-1]))
        vertexCsv = os.path.join(outputDir, "vertices.csv")
        with open(vertexCsv, 'w') as vertexStream:
            if self.gyroBiasSplineDv:
                BSplineIO.saveStates(trueFrameTimes, self.poseSplineDv, self.gyroBiasSplineDv.spline(),
                                     self.accBiasSplineDv.spline(), self.timeOffset, vertexStream)
            else:
                BSplineIO.saveStates(trueFrameTimes, self.poseSplineDv, None, None, self.timeOffset, vertexStream)
            print("  Written simulated states to {}".format(vertexCsv))
        return trueFrameTimes

    def simulate(self, outputDir):
        super(RsCameraImuSimulator, self).simulate(outputDir)
        self.simulateImuData(outputDir)


def noisyValue(x, upperbound, noise):
    if x <= 1 or x >= upperbound - 1:
        noisyx = x
    elif x + noise < 0:
        noisyx = x - noise
    elif x + noise > upperbound:
        noisyx = x - noise
    else:
        noisyx = x + noise
    return noisyx
