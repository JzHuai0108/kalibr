"""
Simulate camera and IMU data given pose and bias bsplines.
This module only supports saving simulation results for one camera
limited by the maplab vertex csv file.
"""

import copy
import math
import os
from random import gauss

import numpy as np
import matplotlib.pyplot as plt

import aslam_backend as aopt
import aslam_cv as acv
import aslam_cameras_april as acv_april
import sm
import kalibr_common as kc
import kalibr_errorterms as ket

from . import BSplineIO

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
        self.chain = kc.CameraChainParameters(args.chain_yaml)
        self.T_imu_cam_list = []
        self.timeOffsetList = []
        self.camGeometryList = []
        numCameras = self.chain.numCameras()
        for i in range(numCameras):
            camConfig = self.chain.getCameraParameters(i)
            print("Camera {}:".format(i))
            camConfig.printDetails()
            printExtraCameraDetails(camConfig)
            # These parameters are set to default values assuming no IMU is present.
            self.T_imu_cam_list.append(sm.Transformation())
            self.timeOffsetList.append(0)
            camera = kc.AslamCamera.fromParameters(camConfig)
            self.camGeometryList.append(camera.geometry)

        targetConfig = kc.CalibrationTargetParameters(args.target_yaml)
        print("Target used in the simulation:")
        targetConfig.printDetails()
        self.targetObservation = None
        self.allTargetCorners = None
        self.setupCalibrationTarget(targetConfig, showExtraction=False, showReproj=False, imageStepping=False)

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
        camId = 0  # only check one camera
        camConfig = self.chain.getCameraParameters(camId)
        timePadding = 2.5 / camConfig.getUpdateRate()
        trueFrameTimes = self.generateStateTimes(camConfig.getUpdateRate(), timePadding)
        state_time = trueFrameTimes[0]
        line_delay = float(camConfig.getLineDelayNanos()) * 1e-9
        resolution = camConfig.getResolution()
        imageCornersNaive = self.naiveMethodToRsProjection(self.camGeometryList[camId], state_time, line_delay,
                                                           self.T_imu_cam_list[camId], resolution, False)
        imageCornersNewton, unusedKeypoints, _ = \
            self.newtonMethodToRsProjection(self.camGeometryList[camId], state_time, line_delay,
                                            self.T_imu_cam_list[camId], resolution, 1.0, False)
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

    def naiveMethodToRsProjection(self, camGeometry, state_time, line_delay, T_imu_cam, resolution, verbose=False):
        """
        This method is not proved theoretically to converge, but it performs as precise as
        Newton's method empirically, though slower.
        return:
            1. Projected image corners according to a rolling shutter model, NX3X1 array.
        """
        imageCornerProjected= list()
        if verbose:
            print('Naive method for state time %.9f' % state_time)
        imageHeight = resolution[1]
        for iota in range(self.targetObservation.getTotalTargetPoint()):
            # get the initial observation
            sm_T_w_c, validPose = getCameraPoseAt(state_time, self.poseSplineDv, T_imu_cam)
            if not validPose:
                continue
            validProjection, lastImagePoint = self.targetObservation.projectATargetPoint(camGeometry, sm_T_w_c, iota, False) # 3x1.
            if not validProjection:
                continue
            numIter = 0
            aborted = False
            if verbose:
                print('lmId', iota, 'iter', numIter, 'image coords', lastImagePoint.T)
            if np.absolute(line_delay) < 1e-8:
                imageCornerProjected.append(lastImagePoint)
                continue
            while numIter < 8:
                currTime = (lastImagePoint[1, 0] - imageHeight * 0.5) * line_delay + state_time
                sm_T_w_cx, validPose = getCameraPoseAt(currTime, self.poseSplineDv, T_imu_cam)
                validProjection, imagePoint = self.targetObservation.projectATargetPoint(camGeometry, sm_T_w_cx, iota, False)
                if not validPose or not validProjection:
                    aborted = True
                    break
                delta = np.absolute(lastImagePoint[1,0] - imagePoint[1,0])
                numIter += 1
                if verbose:
                    print('lmId', iota, 'iter', numIter, 'image coords', imagePoint.T)
                lastImagePoint = imagePoint
                if delta < 1e-3:
                    break
            if verbose:
                print
            if not aborted:
                imageCornerProjected.append(lastImagePoint)
        return np.array(imageCornerProjected)

    def newtonMethodToRsProjection(self, camGeometry, state_time, line_delay, T_imu_cam, resolution,
                                   reprojectionSigma = 1.0, verbose = False):
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
        imageWidth = resolution[0]
        imageHeight = resolution[1]
        for iota in range(numLandmarks):
            sm_T_w_c, validPose = getCameraPoseAt(state_time, self.poseSplineDv, T_imu_cam)
            validProjection, lastImagePoint = self.targetObservation.projectATargetPoint(camGeometry, sm_T_w_c, iota, False) # 3x1.
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
                frameKeypoints.append((iota, kpId, lastImagePoint[0, 0] + xnoise, lastImagePoint[1, 0] + ynoise, reprojectionSigma, 12, -1))
                kpId += 1
                continue
            # solve y=g(y) where y is the vertical projection in pixels
            initialImagePoint = copy.deepcopy(lastImagePoint)
            while numIter < 6:
                # now we have y_0, i.e., lastImagePoint[1, 0], complete the iteration by computing y_1

                # compute g(y_0)
                currTime = (lastImagePoint[1, 0] - imageHeight * 0.5) * line_delay + state_time
                sm_T_w_cx, validPose = getCameraPoseAt(currTime, self.poseSplineDv, T_imu_cam)

                validProjection, imagePoint0 = self.targetObservation.projectATargetPoint(camGeometry, sm_T_w_cx, iota, False)
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
                currTime = (lastImagePoint[1, 0] + eps - imageHeight * 0.5) * line_delay + state_time
                sm_T_w_cx, validPose = getCameraPoseAt(currTime, self.poseSplineDv, T_imu_cam)

                validProjection, imagePoint1 = self.targetObservation.projectATargetPoint(camGeometry, sm_T_w_cx, iota, False)
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
                noisyPoint = [noisyValue(imagePoint0[0, 0], imageWidth, xnoise),
                              noisyValue(imagePoint0[1, 0], imageHeight, ynoise)]
                frameKeypoints.append((iota, kpId, noisyPoint[0], noisyPoint[1], reprojectionSigma, 12, -1))
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

        landmarkObservations = dict()
        for iota in range(self.targetObservation.getTotalTargetPoint()):
            landmarkObservations[iota] = list()
        frameKeypointMap = dict()
        for vertexId, state_time in enumerate(trueFrameTimes):
            frameKeypointMap[vertexId] = list()

        for cameraIndex in range(self.chain.numCameras()):
            cameraTimeOffset = self.chain.getTimeshiftCamImu(cameraIndex)
            T_imu_cam = self.T_imu_cam_list[cameraIndex]
            camConfig = self.chain.getCameraParameters(cameraIndex)
            resolution = camConfig.getResolution()
            imageNoise = camConfig.getImageNoise()
            linedelaysecs = float(camConfig.getLineDelayNanos()) * 1e-9
            for vertexId, frameTime in enumerate(trueFrameTimes):
                rawFrameTime = frameTime - cameraTimeOffset
                _, noisyKeypoints, keypointOffsets = self.newtonMethodToRsProjection(
                    self.camGeometryList[cameraIndex], frameTime, linedelaysecs, T_imu_cam, resolution, imageNoise)
                imageCornerOffsetNorms += keypointOffsets
                if vertexId % 300 == 0:
                    print('  Projected {:d} target landmarks for state at {:.9f}'.format(len(noisyKeypoints), frameTime))
                for keypoint in noisyKeypoints:
                    landmarkObservations[keypoint[0]].append((vertexId, cameraIndex, keypoint[1]))
                frameKeypointMap[vertexId].append((acv.Time(rawFrameTime), cameraIndex, noisyKeypoints))

        observationCsv = os.path.join(outputDir, "observations.csv")
        kc.VimapCsvWriter.saveObservations(landmarkObservations, observationCsv)

        trackCsv = os.path.join(outputDir, "tracks.csv")
        kc.VimapCsvWriter.saveTracks(frameKeypointMap, trackCsv)

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

    def computeCameraRate(self, camId):
        camConfig = self.chain.getCameraParameters(camId)
        lineDelay = camConfig.getLineDelayNanos()
        imageHeight = camConfig.getResolution()[1]
        maxFrameRate = math.floor(1e9 / ((lineDelay + 1000) * imageHeight))
        cameraRate = min(maxFrameRate, camConfig.getUpdateRate())
        return cameraRate

    def simulateStates(self, outputDir):
        cameraRate = self.computeCameraRate(0)
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

    def saveCameraIntrinsics(self, outputDir):
        """
        save the intrinsics for only one camera.
        :param outputDir:
        :return:
        """
        yamlfile = os.path.join(outputDir, "camchain.yaml")
        mycamchain = kc.CameraChainParameters(yamlfile, createYaml=True)
        numCams = self.chain.numCameras()
        for i in range(numCams):
            camParams = self.chain.getCameraParameters(i)
            mycamchain.addCameraAtEnd(camParams)
        mycamchain.writeYaml()

    def simulate(self, outputDir):
        self.simulateLandmarks(outputDir)
        trueFrameTimes = self.simulateStates(outputDir)
        print("Simulating camera observations...")
        self.simulateCameraObservations(trueFrameTimes, outputDir)
        self.saveCameraIntrinsics(outputDir)

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

        self.T_imu_cam_list = []
        self.timeOffsetList = []
        for camId in range(self.chain.numCameras()):
            T_cam_imu = self.chain.getExtrinsicsImuToCam(camId)
            self.T_imu_cam_list.append(T_cam_imu.inverse())
            self.timeOffsetList.append(self.chain.getTimeshiftCamImu(camId))

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

        gravity = self.imuConfig.getGravityInTarget()
        gravityDv = aopt.EuclideanDirection(np.array(gravity))
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
        cameraRate = self.computeCameraRate(0)
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
        """
        save poses for every image at timestamp in camera clock.
        """
        cameraRate = self.computeCameraRate(0)
        timePadding = 2.5 / cameraRate
        trueFrameTimes = self.generateStateTimes(cameraRate, timePadding)

        print('Simulating states at camera rate {}...'.format(cameraRate))
        print("  Camera frame true start time {:.9f} and true finish time {:.9f}".format(
            trueFrameTimes[0], trueFrameTimes[-1]))
        vertexCsv = os.path.join(outputDir, "vertices.csv")
        with open(vertexCsv, 'w') as vertexStream:
            if self.gyroBiasSplineDv:
                BSplineIO.saveStates(trueFrameTimes, self.poseSplineDv, self.gyroBiasSplineDv.spline(),
                                     self.accBiasSplineDv.spline(), self.timeOffsetList[0], vertexStream)
            else:
                BSplineIO.saveStates(trueFrameTimes, self.poseSplineDv, None, None, self.timeOffsetList[0], vertexStream)
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
