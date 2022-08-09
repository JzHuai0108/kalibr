#encoding:UTF-8
from __future__ import print_function
import time

import sm
import aslam_backend as aopt
import aslam_cv_backend as acvb
import aslam_cv as acv
import aslam_splines as asp
import incremental_calibration as inc
from kalibr_common import ConfigReader as cr
import bsplines
import kalibr_common as kc
import numpy as np
import multiprocessing
import sys
import gc
import math
from .ReprojectionErrorKnotSequenceUpdateStrategy import *
from .RsPlot import plotSpline
from .RsPlot import plotSplineValues
import pylab as pl
import pdb
from kalibr_imu_camera_calibration import BSplineIO
from kalibr_imu_camera_calibration import IccSensors as sens

# make numpy print prettier
np.set_printoptions(suppress=True)

CALIBRATION_GROUP_ID = 0
TRANSFORMATION_GROUP_ID = 1
LANDMARK_GROUP_ID = 2
HELPER_GROUP_ID = 3

class RsCalibratorConfiguration(object):
    deltaX = 1e-8
    """Stopping criterion for the optimizer on x"""

    deltaJ = 1e-4
    """Stopping criterion for the optimizer on J"""

    maxNumberOfIterations = 20
    """Maximum number of iterations of the batch optimizer"""

    maxKnotPlacementIterations = 10
    """Maximum number of iterations to take in the adaptive knot-placement step"""

    adaptiveKnotPlacement = True
    """Whether to enable adaptive knot placement"""

    knotUpdateStrategy = ReprojectionErrorKnotSequenceUpdateStrategy
    """The adaptive knot placement strategy to use"""

    timeOffsetConstantSparsityPattern = 0.08
    """A time offset to pad the blocks generated in the hessian/jacobian to ensure a constant symbolic representation
    of the batch estimation problem, even when a change in the shutter timing shifts the capture time to another
    spline segment.
    """

    inverseFeatureCovariance = 1/0.26
    """The inverse covariance of the feature detector. Used to standardize the error terms."""

    estimateParameters = {'shutter': True, 'intrinsics': True, 'distortion': True,
                          'timeOffset': False, 'pose': True, 'landmarks': False}
    """Which parameters to estimate. Dictionary with shutter, intrinsics, distortion, pose, landmarks as bool"""

    splineOrder = 4
    """Order of the spline to use for ct-parametrization"""

    timeOffsetPadding = 0.05
    """Time offset to add to the beginning and end of the spline to ensure we remain
    in-bounds while estimating time-depending parameters that shift the spline.
    """

    numberOfKnots = None
    """Set to an integer to start with a fixed number of uniformly distributed knots on the spline."""

    W = None
    """6x6 diagonal matrix with a weak motion prior"""

    framerate = 30
    """The approximate framerate of the camera. Required as approximate threshold in adaptive
    knot placement and for initializing a knot sequence if no number of knots is given.
    """

    chain_yaml = None
    """Camera system configuration yaml. If provided, it will be used to initialize the camera projection and distortion parameters!"""

    saveSamplePoses = False

    reprojectFrameIndex = 0

    recoverCov = False

    def validate(self, isRollingShutter):
        """Validate the configuration."""
        # only rolling shutters can be estimated
        if (not isRollingShutter):
            self.estimateParameters['shutter'] = False
            self.adaptiveKnotPlacement = False


class ImuDataDescription(object):
    def __init__(self, bagfiles, bag_from_to, perform_sync, constant_bias):
        self.bagfile = bagfiles
        self.bag_from_to = bag_from_to
        self.perform_synchronization = perform_sync
        self.constant_bias = constant_bias


class RsCalibrator(object):

    __observations = None
    """Store the list of observations"""

    __cameraGeometry = None
    """The geometry container of which the calibration is performed."""

    __camera = None
    """The camera geometry itself."""

    __camera_dv = None
    """The camera design variable"""

    __cameraModelFactory = None
    """Factory object that can create a typed objects for a camera (error terms, frames, design variables etc)"""

    __poseSpline = None
    """The spline describing the pose of the camera"""

    __poseSpline_dv = None
    """The design variable representation of the pose spline of the camera"""

    __config = None
    """Configuration container \see RsCalibratorConfiguration"""

    __frames = []
    """All frames observed"""

    __reprojection_errors = []
    """Reprojection errors of the latest optimizer iteration"""

    __std_camera = None

    __ImuList = None

    __imageHeight = 0

    def calibrate(self,
        cameraGeometry,
        observations,
        config
    ):
        """
        A Motion regularization term is added with low a priori knowledge to avoid
        diverging parts in the spline of too many knots are selected/provided or if
        no image information is available for long sequences and to regularize the
        last few frames (which typically contain no image information but need to have
        knots to /close/ the spline).

        Kwargs:
            cameraGeometry (kcc.CameraGeometry): a camera geometry object with an initialized target
            observations ([]: The list of observation \see extractCornersFromDataset
            config (RsCalibratorConfiguration): calibration configuration
        """

        ## set internal objects
        self.__observations = observations
        self.__cameraGeometry = cameraGeometry
        self.__cameraModelFactory = cameraGeometry.model
        self.__camera_dv = cameraGeometry.dv
        self.__camera = cameraGeometry.geometry
        self.__config = config
        self.__config.validate(self.__isRollingShutter())

        # obtain initial guesses for extrinsics and intrinsics
        if (not self.__generateIntrinsicsInitialGuess()):
            sm.logError("Could not generate initial guess.")

        # obtain the extrinsic initial guess for every observation
        self.__generateExtrinsicsInitialGuess()

        # set the value for the motion prior term or uses the defaults
        W = self.__getMotionModelPriorOrDefault()

        times = [observation.time().toSec() for observation in self.__observations]
        times = np.sort(times)
        deltaTimes = np.diff(times)
        self.__config.framerate = 1.0 / np.median(deltaTimes)

        self.__poseSpline = self.__generateInitialSpline(
            self.__config.splineOrder,
            self.__config.timeOffsetPadding,
            self.__config.numberOfKnots,
            self.__config.framerate
        )

        # build estimator problem
        optimisation_problem = self.__buildOptimizationProblem(W)

        status = self.__runOptimization(
            optimisation_problem,
            self.__config.deltaJ,
            self.__config.deltaX,
            self.__config.maxNumberOfIterations
        )

        # continue with knot replacement
        if self.__config.adaptiveKnotPlacement:
            knotUpdateStrategy = self.__config.knotUpdateStrategy(self.__config.framerate)

            for iteration in range(self.__config.maxKnotPlacementIterations):

                # generate the new knots list
                [knots, requiresUpdate] = knotUpdateStrategy.generateKnotList(
                    self.__reprojection_errors,
                    self.__poseSpline_dv.spline()
                )
                # if no new knotlist was generated, we are done.
                if (not requiresUpdate):
                    break;

                # otherwise update the spline dv and rebuild the problem
                self.__poseSpline = knotUpdateStrategy.getUpdatedSpline(self.__poseSpline_dv.spline(), knots, self.__config.splineOrder)

                optimisation_problem = self.__buildOptimizationProblem(W)
                status = self.__runOptimization(
                    optimisation_problem,
                    self.__config.deltaJ,
                    self.__config.deltaX,
                    self.__config.maxNumberOfIterations
                )

        if self.__config.saveSamplePoses:
            self.__saveBSplinePoses()

        self.__printResults()
        self.__saveParametersYaml()
        if status and self.__config.recoverCov:
            self.recoverCovariance(optimisation_problem)

    def __saveBSplinePoses(self):
            bspline = self.__poseSpline_dv.spline()
            rate = 100
            interval = 1.0 / rate
            padding = 2
            samplePoseTimes = np.arange(bspline.t_min() + padding, bspline.t_max() - padding, interval)
            refPoseStream = open("sampled_poses.txt", 'w')
            print("%poses {} Hz from the RS calibrator B-splines: time (sec), T_w_c (txyz, qxyzw).".format(rate), file=refPoseStream)
            BSplineIO.sampleAndSaveBSplinePoses(samplePoseTimes, self.__poseSpline_dv, stream=refPoseStream)
            refPoseStream.close()


    def __generateExtrinsicsInitialGuess(self):
        """Estimate the pose of the camera with a PnP solver. Call after initializing the intrinsics"""
        # estimate and set T_c in the observations
        for idx, observation in enumerate(self.__observations):
            (success, T_t_c) = self.__camera.estimateTransformation(observation)
            if (success):
                observation.set_T_t_c(T_t_c)
            else:
                sm.logWarn("Could not estimate T_t_c for observation at index {0}".format(idx))

        return

    def __generateIntrinsicsInitialGuess(self):
        """
        Get an initial guess for the camera geometry (intrinsics, distortion). Distortion is typically left as 0,0,0,0.
        The parameters of the geometryModel are updated in place.
        """
        if self.__config.chain_yaml:
            camchain = kc.CameraChainParameters(self.__config.chain_yaml)
            camConfig = camchain.getCameraParameters(0)
            aslamCamera = kc.AslamCamera.fromParameters(camConfig)
            resolution = camConfig.getResolution()
            self.__imageHeight = resolution[1]
            self.__camera = aslamCamera.geometry
            self.__camera_dv = self.__cameraModelFactory.designVariable(self.__camera)
            status = True
            print('External projection and distortion parameters {}'.format(
                self.__camera.getParameters(True, True, True).T))
        else:
            if self.__isRollingShutter():
                self.__imageHeight = self.__observations[0].imRows()
                self.__camera.shutter().setParameters(np.array([1.0 / self.__config.framerate / float(sensorRows)]))
            status = self.__camera.initializeIntrinsics(self.__observations)
            print('Initial projection and distortion parameters {}'.format(self.__camera.getParameters(True, True, True).T))
        return status

    def __getMotionModelPriorOrDefault(self):
        """Get the motion model prior or the default value"""
        W = self.__config.W
        if W is None:
            W = np.eye(6)
            W[:3,:3] *= 1e-3
            W[3:,3:] *= 1
            W *= 1e-2
        return W

    def __generateInitialSpline(self, splineOrder, timeOffsetPadding, numberOfKnots = None, framerate = None):
        poseSpline = bsplines.BSplinePose(splineOrder, sm.RotationVector())

        # Get the observation times.
        times = np.array([observation.time().toSec() for observation in self.__observations ])
        # get the pose values of the initial transformations at observation time
        curve = np.matrix([ poseSpline.transformationToCurveValue( observation.T_t_c().T() ) for observation in self.__observations]).T
        # make sure all values are well defined
        if np.isnan(curve).any():
            raise RuntimeError("Nans in curve values")
            sys.exit(0)
        # Add 2 seconds on either end to allow the spline to slide during optimization
        times = np.hstack((times[0] - (timeOffsetPadding * 2.0), times, times[-1] + (timeOffsetPadding * 2.0)))
        curve = np.hstack((curve[:,0], curve, curve[:,-1]))

        self.__ensureContinuousRotationVectors(curve)

        seconds = times[-1] - times[0]

        # fixed number of knots
        if (numberOfKnots is not None):
            knots = numberOfKnots
        # otherwise with framerate estimate
        else:
            knots = int(round(seconds * framerate/3))

        print("")
        print("Initializing a pose spline with %d knots (%f knots per second over %f seconds)" % ( knots, knots/seconds, seconds))
        poseSpline.initPoseSplineSparse(times, curve, knots, 1e-4)
        return poseSpline


    def loadImu(self, imu_yaml, imu_model, imuDescription):
        if imu_yaml is None:
            self.__ImuList = None
            return
        imus = list()
        imuConfig = kc.ImuParameters(imu_yaml)
        imuConfig.printDetails()

        if imu_model != "calibrated":
            raise Exception("Only calibrated IMU model is supported for simplicity!")
        imus.append(sens.IccImu(imuConfig, imuDescription, isReferenceImu=True, estimateTimedelay=False))
        self.__ImuList = imus


    def __addImuErrors(self, problem):
        poseSpline = self.__poseSpline_dv.spline()
        for imu in self.__ImuList:
            splineOrder = 4
            biasKnotsPerSecond = 5
            imu.initBiasSplines(poseSpline, splineOrder, biasKnotsPerSecond)

        # estimate gravity in the world coordinate frame as the mean specific force.
        R_i_c = np.identity(3)
        if self.__config.chain_yaml:
            camchain = kc.CameraChainParameters(self.__config.chain_yaml)
            T_cam_imu = camchain.getExtrinsicsImuToCam(0)
            T_imu_cam = T_cam_imu.inverse()
            R_i_c = sm.quat2r(T_imu_cam.q())

        a_w = []
        for im in self.__ImuList[0].imuData:
            tk = im.stamp.toSec()
            if tk > poseSpline.t_min() and tk < poseSpline.t_max():
                a_w.append(np.dot(poseSpline.orientation(tk), np.dot(R_i_c, - im.alpha)))
        if len(a_w) == 0:
            print("poseSpline t_min {}, t_max {}.".format(poseSpline.t_min(), poseSpline.t_max()))
            print("IMU data t_min {}, t_max {}.".format(self.__ImuList[0].imuData[0].stamp.toSec(),
                                                        self.__ImuList[0].imuData[-1].stamp.toSec()))
            raise IOError('No corresponding IMU data were found.')
        
        mean_a_w = np.mean(np.asarray(a_w).T, axis=1)
        # A rough gravity magnitude is OK for RS camera calibration with loose IMU constraints.
        gravity_w = mean_a_w / np.linalg.norm(mean_a_w) * 9.80655
        print("Gravity was intialized to {} [m/s^2]".format(gravity_w))

        # Add the calibration target orientation design variable. (expressed as gravity vector in target frame)
        self.gravityDv = aopt.EuclideanDirection(gravity_w)
        self.gravityExpression = self.gravityDv.toExpression()
        self.gravityDv.setActive(True)
        problem.addDesignVariable(self.gravityDv, HELPER_GROUP_ID)

        for imu in self.__ImuList:
            imu.addDesignVariables(problem)

        huberAccel = -1
        huberGyro = -1
        gyroNoiseScale = 1.0
        accelNoiseScale = 1.0
        for imu in self.__ImuList:
            imu.addAccelerometerErrorTerms(problem, self.__poseSpline_dv, self.gravityExpression, mSigma=huberAccel, accelNoiseScale=accelNoiseScale)
            imu.addGyroscopeErrorTerms(problem, self.__poseSpline_dv, mSigma=huberGyro, gyroNoiseScale=gyroNoiseScale, g_w=self.gravityExpression)
            imu.addBiasMotionTerms(problem)

    def __buildOptimizationProblem(self, W):
        """Build the optimisation problem"""
        problem = inc.CalibrationOptimizationProblem()

        # Initialize all design variables.
        self.__initPoseDesignVariables(problem)

        #####
        ## build error terms and add to problem

        # store all frames
        self.__frames = []
        self.__reprojection_errors = []

        # This code assumes that the order of the landmarks in the observations
        # is invariant across all observations. At least for the chessboards it is true.

        #####
        # add all the landmarks once
        landmarks = []
        landmarks_expr = []
        target = self.__cameraGeometry.ctarget.detector.target()
        for i in range(0, target.size()):
            landmark_w_dv = aopt.HomogeneousPointDv(sm.toHomogeneous(target.point(i)))
            landmark_w_dv.setActive(self.__config.estimateParameters['landmarks'])
            landmarks.append(landmark_w_dv)
            landmarks_expr.append(landmark_w_dv.toExpression())
            problem.addDesignVariable(landmark_w_dv, LANDMARK_GROUP_ID)

        #####
        # activate design variables
        self.__camera_dv.setActive(
            self.__config.estimateParameters['intrinsics'],
            self.__config.estimateParameters['distortion'],
            self.__config.estimateParameters['shutter']
        )

        #####
        # Add design variables

        # add the camera design variables last for optimal sparsity patterns
        problem.addDesignVariable(self.__camera_dv.shutterDesignVariable(), CALIBRATION_GROUP_ID)
        problem.addDesignVariable(self.__camera_dv.projectionDesignVariable(), CALIBRATION_GROUP_ID)
        problem.addDesignVariable(self.__camera_dv.distortionDesignVariable(), CALIBRATION_GROUP_ID)

        #####
        # Regularization term / motion prior
        if self.__ImuList:
            self.__addImuErrors(problem)
            W *= 1e-2
        # Add the time delay design variable.
        self.cameraTimeToImuTimeDv = aopt.Scalar(0.0)
        self.cameraTimeToImuTimeDv.setActive(self.__config.estimateParameters['timeOffset'])
        problem.addDesignVariable(self.cameraTimeToImuTimeDv, CALIBRATION_GROUP_ID)


        motionError = asp.BSplineMotionError(self.__poseSpline_dv, W)
        problem.addErrorTerm(motionError)

        dummyPoint = np.array([0, self.__imageHeight / 2])
        centerRowTemporalOffset = self.__camera_dv.temporalOffset(dummyPoint)

        #####
        # add a reprojection error for every corner of each observation
        for observation in self.__observations:
            # only process successful observations of a pattern
            if (observation.hasSuccessfulObservation()):
                # add a frame
                frame = self.__cameraModelFactory.frameType()
                frame.setGeometry(self.__camera)
                frame.setTime(observation.time())
                self.__frames.append(frame)

                #####
                # add an error term for every observed corner
                corner_id_list = observation.getCornersIdx()
                for index, point in enumerate(observation.getCornersImageFrame()):
                    # keypoint time offset by line delay as expression type
                    keypoint_time = self.cameraTimeToImuTimeDv.toExpression() + \
                                    self.__camera_dv.keypointTime(frame.time(), point) - centerRowTemporalOffset

                    # from target to world transformation.
                    T_w_t = self.__poseSpline_dv.transformationAtTime(
                        keypoint_time,
                        self.__config.timeOffsetConstantSparsityPattern,
                        self.__config.timeOffsetConstantSparsityPattern
                    )
                    T_t_w = T_w_t.inverse()

                    # transform target point to camera frame
                    p_t = T_t_w * landmarks_expr[corner_id_list[index]]

                    # create the keypoint
                    keypoint = acv.Keypoint2()
                    keypoint.setMeasurement(point)
                    inverseFeatureCovariance = self.__config.inverseFeatureCovariance
                    keypoint.setInverseMeasurementCovariance(np.eye(len(point)) * inverseFeatureCovariance)
                    frame.addKeypoint(keypoint)

                    # create reprojection error
                    reprojection_error = self.__buildErrorTerm(
                        frame,
                        index,
                        p_t,
                        self.__camera_dv,
                        self.__poseSpline_dv
                    )
                    self.__reprojection_errors.append(reprojection_error)
                    problem.addErrorTerm(reprojection_error)

        return problem

    def getReprojectedCorners(self, frameIndex):
        """
        Reproject detected corners of a frame with the RS and the GS.
        :param frameIndex: the index of fhe image frame within the used images for reprojection.
        :return: A NX6 array. Each row corresponds to an observed landmark,
        The columns correspond to reprojected RS point, reprojected GS point, and the detected corner.
        """
        observation = self.__observations[frameIndex]
        print("Reprojecting corners for frame {} at time {}.".format(frameIndex, observation.time().toSec()))
        imageCornerPoints = np.array(observation.getCornersImageFrame())  # Nx2

        # Rolling shutter projections
        # add all the landmarks once
        landmarks = []
        landmarks_expr = []
        target = self.__cameraGeometry.ctarget.detector.target()
        for i in range(0, target.size()):
            landmark_w_dv = aopt.HomogeneousPointDv(sm.toHomogeneous(target.point(i)))
            landmark_w_dv.setActive(self.__config.estimateParameters['landmarks'])
            landmarks.append(landmark_w_dv)
            landmarks_expr.append(landmark_w_dv.toExpression())

        frame = self.__cameraModelFactory.frameType()
        frame.setGeometry(self.__camera)
        frame.setTime(observation.time())
        # build an error term for every observed corner
        corner_id_list = observation.getCornersIdx()
        predictedMeasurements = list()
        for index, point in enumerate(observation.getCornersImageFrame()):
            # keypoint time offset by line delay as expression type
            keypoint_time = self.__camera_dv.keypointTime(frame.time(), point)

            # from target to world transformation.
            T_w_t = self.__poseSpline_dv.transformationAtTime(
                keypoint_time,
                self.__config.timeOffsetConstantSparsityPattern,
                self.__config.timeOffsetConstantSparsityPattern
            )
            T_t_w = T_w_t.inverse()

            # transform target point to camera frame
            p_t = T_t_w * landmarks_expr[corner_id_list[index]]

            # create the keypoint
            keypoint_index = frame.numKeypoints()
            keypoint = acv.Keypoint2()
            keypoint.setMeasurement(point)
            inverseFeatureCovariance = self.__config.inverseFeatureCovariance
            keypoint.setInverseMeasurementCovariance(np.eye(len(point)) * inverseFeatureCovariance)
            frame.addKeypoint(keypoint)

            rerr = self.__cameraModelFactory.reprojectionError(frame, keypoint_index, p_t, self.__camera_dv)
            rerr.evaluateError()
            predictedMeas = imageCornerPoints[index, :].T - rerr.error()
            predictedMeasurements.append(predictedMeas)
        rsImageCornerProjected = np.array(predictedMeasurements)  # NX2

        # Global shutter projections
        # Build a transformation expression for the time.
        cameraTimeToImuTimeDv = aopt.Scalar(0.0)

        proj = self.__camera_dv.projectionDesignVariable().value()
        sensorRows = proj.rv()
        halfSensorRows = sensorRows / 2

        lineDelay = self.__camera_dv.shutterDesignVariable().value().lineDelay()
        frameTime = cameraTimeToImuTimeDv.toExpression() + observation.time().toSec() + halfSensorRows * lineDelay
        frameTimeScalar = frameTime.toScalar()
        # as we are applying an initial time shift outside the optimization so
        # we need to make sure that we dont add data outside the spline definition
        if frameTimeScalar <= self.__poseSpline_dv.spline().t_min() or frameTimeScalar >= self.__poseSpline_dv.spline().t_max():
            return np.array([])

        T_w_c = self.__poseSpline_dv.transformationAtTime(
            frameTime,
            self.__config.timeOffsetConstantSparsityPattern,
            self.__config.timeOffsetConstantSparsityPattern
        )
        # Simple approach to reproject landmarks with float numbers.
        observation.set_T_t_c(sm.Transformation(T_w_c.toTransformationMatrix()))
        gsImageCornerProjected = np.array(observation.getCornerReprojection(self.__camera))  # Nx2
        if gsImageCornerProjected.shape[0] == 0: # when simulated observations are used, the GS corners can be empty.
            return np.concatenate((rsImageCornerProjected, imageCornerPoints), axis=1)
        return np.concatenate((rsImageCornerProjected, gsImageCornerProjected, imageCornerPoints), axis=1)

    def __buildErrorTerm(self, frame, keypoint_index, p_t, camera_dv, poseSpline_dv):
        """
        Build an error term that considers the shutter type. A Global Shutter camera gets the standard reprojection error
        a Rolling Shutter gets the adaptive covariance error term that considers the camera motion.
        """
        # it is a global shutter camera -> no covariance error
        if (self.__isRollingShutter()):
            return self.__cameraModelFactory.reprojectionErrorAdaptiveCovariance(
                frame,
                keypoint_index,
                p_t,
                camera_dv,
                poseSpline_dv
            )
        else:
            return self.__cameraModelFactory.reprojectionError(
                frame,
                keypoint_index,
                p_t,
                camera_dv
            )

    def __ensureContinuousRotationVectors(self, curve):
        """
        Ensures that the rotation vector does not flip and enables a continuous trajectory modeling.
        Updates curves in place.
        """
        for i in range(1, curve.shape[1]):
            previousRotationVector = curve[3:6,i-1]
            r = curve[3:6,i]
            angle = np.linalg.norm(r)
            axis = r/angle
            best_r = r
            best_dist = np.linalg.norm( best_r - previousRotationVector)

            for s in range(-3,4):
                aa = axis * (angle + math.pi * 2.0 * s)
                dist = np.linalg.norm( aa - previousRotationVector )
                if dist < best_dist:
                    best_r = aa
                    best_dist = dist
            curve[3:6,i] = best_r

    def __initPoseDesignVariables(self, problem):
        """Get the design variable representation of the pose spline and add them to the problem"""
        # get the design variable
        self.__poseSpline_dv = asp.BSplinePoseDesignVariable(self.__poseSpline)
        # activate all contained dv and add to problem
        for i in range(0, self.__poseSpline_dv.numDesignVariables()):
            dv = self.__poseSpline_dv.designVariable(i)
            dv.setActive(self.__config.estimateParameters['pose'])
            problem.addDesignVariable(dv, TRANSFORMATION_GROUP_ID)

    def __runOptimization(self, problem ,deltaJ, deltaX, maxIt):
        """Run the given optimization problem problem"""

        print("run new optimisation with initial values:")
        self.__printResults()

        # verbose and choldmod solving with schur complement trick
        options = aopt.Optimizer2Options()
        options.verbose = True
        options.linearSolver = aopt.BlockCholeskyLinearSystemSolver()
        options.doSchurComplement = True

        # stopping criteria
        options.maxIterations = maxIt
        options.convergenceDeltaJ = deltaJ
        options.convergenceDeltaX = deltaX

        # use the dogleg trustregion policy
        options.trustRegionPolicy = aopt.DogLegTrustRegionPolicy()

        # create the optimizer
        optimizer = aopt.Optimizer2(options)
        optimizer.setProblem(problem)

        # go for it:
        status = optimizer.optimize()
        if status:
            corners = self.getReprojectedCorners(self.__config.reprojectFrameIndex)
            np.savetxt("reprojected_corners_{}.txt".format(self.__config.reprojectFrameIndex), corners)

        return status

    def recoverCovariance(self, problem):
        """Computing covariance takes so long because of the complex prior BSplineMotionError"""
        #Covariance ordering (=dv ordering)
        #ORDERING:   N=num cams
        #            camera -->  sum(sub) * N
        #                a) shutter    --> 1
        #                b) projection --> omni:5, pinhole: 4
        #                c) distortion --> 4
        tic = time.time()
        estimator = inc.IncrementalEstimator(CALIBRATION_GROUP_ID)
        rval = estimator.addBatch(problem, True)
        est_stds = np.sqrt(estimator.getSigma2Theta().diagonal())
        toc = time.time()
        elapsed = toc - tic
        print("Covariance recovery takes {} secs".format(elapsed))
        #split and store the variance
        std_camera = list()
        offset=0
        nt = self.__camera.geometry.minimalDimensionsShutter() +  \
            self.__camera.geometry.minimalDimensionsProjection() +  \
            self.__camera.geometry.minimalDimensionsDistortion()      
        std_camera.extend(est_stds[offset:offset+nt].flatten().tolist())
        offset = offset+nt
        self.__std_camera = std_camera
        print('std_camera: {}'.format(std_camera))

    def __isRollingShutter(self):
        return self.__cameraModelFactory.shutterType == acv.RollingShutter

    def __printResults(self):
        shutter = self.__camera_dv.shutterDesignVariable().value()
        proj = self.__camera_dv.projectionDesignVariable().value()
        dist = self.__camera_dv.distortionDesignVariable().value()
        dt = self.cameraTimeToImuTimeDv.toScalar()
        print("")
        if not self.__std_camera:
            if (self.__isRollingShutter()):
                print("LineDelay: {}".format(shutter.lineDelay()))
            p = proj.getParameters().flatten()
            print("Intrinsics: {}".format(p))
            d = dist.getParameters().flatten()
            print("Distortion: {}".format(d))
        else:
            if (self.__isRollingShutter()):
                print("LineDelay: {} +/- {}".format(shutter.lineDelay(), self.__std_camera[0]))
            p = proj.getParameters().flatten()
            print("Intrinsics: {} +/- {}".format(p, self.__std_camera[1:1+p.shape[0]]))
            d = dist.getParameters().flatten()
            print("Distortion: {} +/- {}".format(d, self.__std_camera[1+p.shape[0]:]))
        print("timeshift_cam_imu: {}".format(dt))

    def __saveParametersYaml(self):
        # Create new config file
        try:
            regulartag = self.__cameraGeometry.dataset.bagfile.translate({ord(c):None for c in "<>:/\|?*"})
        except TypeError:
            regulartag = self.__cameraGeometry.dataset.bagfile.translate(None, "<>:/\|?*")
        bagtag = regulartag.replace('.bag', '', 1)
        resultFile = "camchain-" + bagtag + ".yaml"
        chain = cr.CameraChainParameters(resultFile, createYaml=True)
        camParams = cr.CameraParameters(resultFile, createYaml=True)
        camParams.setRosTopic(self.__cameraGeometry.dataset.topic)

        # Intrinsics
        cameraModels = {
            # Rolling shutter
            acvb.DistortedPinholeRs: 'pinhole',
            acvb.EquidistantPinholeRs: 'pinhole',
            acvb.DistortedOmniRs: 'omni',
            # Global shutter
            acvb.DistortedPinhole: 'pinhole',
            acvb.EquidistantPinhole: 'pinhole',
            acvb.DistortedOmni: 'omni'}
        cameraModel = cameraModels[self.__cameraGeometry.model]
        proj = self.__camera_dv.projectionDesignVariable().value()
        camParams.setIntrinsics(cameraModel, proj.getParameters().flatten())
        camParams.setResolution([proj.ru(), proj.rv()])

        # Distortion
        distortionModels = {
            # Rolling shutter
            acvb.DistortedPinholeRs: 'radtan',
            acvb.EquidistantPinholeRs: 'equidistant',
            acvb.DistortedOmniRs: 'radtan',
            # Global shutter
            acvb.DistortedPinhole: 'radtan',
            acvb.EquidistantPinhole: 'equidistant',
            acvb.DistortedOmni: 'radtan'}
        distortionModel = distortionModels[self.__cameraGeometry.model]
        dist = self.__camera_dv.distortionDesignVariable().value()
        camParams.setDistortion(distortionModel, dist.getParameters().flatten())

        # Shutter
        shutter = self.__camera_dv.shutterDesignVariable().value()
        camParams.setLineDelay(shutter.lineDelay())

        chain.addCameraAtEnd(camParams)
        chain.writeYaml()