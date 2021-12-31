import time

import aslam_backend as aopt
import aslam_splines as asp
from . import IccUtil as util
import incremental_calibration as inc
import kalibr_common as kc
import sm

import gc
import numpy as np
import multiprocessing
import sys

# make numpy print prettier
np.set_printoptions(suppress=True)

CALIBRATION_GROUP_ID = 0
HELPER_GROUP_ID = 1

class IccCalibratorConfiguration(object):
    estimateParameters = {'shutter': False, 'intrinsics': False, 'distortion': False, 'timeOffset': False,
                          'chainExtrinsics': False, 'gravityLength': False, 'pose': True, 'landmarks': False}
    initialGravityEstimate = np.array([0.0,9.81,0.0])


class IccCalibrator(object):
    CALIBRATION_GROUP_ID = CALIBRATION_GROUP_ID
    HELPER_GROUP_ID = HELPER_GROUP_ID
    def __init__(self, config):
        self.ImuList = []
        self.__config = config
        self.CameraChain = None
        self.poseDv = None
        self.gravityDv = None
        self.gravityExpression = None
        self.problem = None
        self.optimizer = None
        self.optimizationResult = None

    def getEstimateParameters(self):
        return self.__config.estimateParameters

    def initDesignVariables(self, problem, poseSpline):
        # Initialize the system pose spline (always attached to imu0) 
        self.poseDv = asp.BSplinePoseDesignVariable( poseSpline )
        self.addSplineDesignVariables(problem, self.poseDv)

        # Add the calibration target orientation design variable. (expressed as gravity vector in target frame)
        if self.__config.estimateParameters['gravityLength']:
            self.gravityDv = aopt.EuclideanPointDv( self.__config.initialGravityEstimate )
        else:
            self.gravityDv = aopt.EuclideanDirection( self.__config.initialGravityEstimate )
        self.gravityExpression = self.gravityDv.toExpression()  
        self.gravityDv.setActive( True )
        problem.addDesignVariable(self.gravityDv, HELPER_GROUP_ID)
        
        #Add all DVs for all IMUs
        for imu in self.ImuList:
            imu.addDesignVariables( problem )
        
        #Add all DVs for the camera chain    
        self.CameraChain.addDesignVariables( problem, self.__config.estimateParameters )

    def addPoseMotionTerms(self, problem, tv, rv):
        wt = 1.0/tv;
        wr = 1.0/rv
        W = np.diag([wt,wt,wt,wr,wr,wr])
        asp.addMotionErrorTerms(problem, self.poseDv, W, errorOrder)
        
    #add camera to sensor list (create list if necessary)
    def registerCamChain(self, sensor):
        self.CameraChain = sensor

    def registerImu(self, sensor):
        self.ImuList.append( sensor )
            
    def buildProblem( self, 
                      splineOrder=6,
                      biasSplineOrder=4,
                      poseKnotsPerSecond=70, 
                      biasKnotsPerSecond=70, 
                      doPoseMotionError=False, 
                      mrTranslationVariance=1e6,
                      mrRotationVariance=1e5,
                      doBiasMotionError=True,
                      blakeZisserCam=-1,
                      huberAccel=-1,
                      huberGyro=-1,
                      maxIterations=20,
                      gyroNoiseScale=1.0,
                      accelNoiseScale=1.0,
                      timeOffsetPadding=0.02,
                      timeOffsetConstantSparsityPattern=0.08,
                      verbose=False  ):

        print("\tSpline order: %d" % (splineOrder))
        print("\tBias spline order: %d" % (biasSplineOrder))
        print("\tPose knots per second: %d" % (poseKnotsPerSecond))
        print("\tDo pose motion regularization: %s" % (doPoseMotionError))
        print("\t\txddot translation variance: %f" % (mrTranslationVariance))
        print("\t\txddot rotation variance: %f" % (mrRotationVariance))
        print("\tBias knots per second: %d" % (biasKnotsPerSecond))
        print("\tDo bias motion regularization: %s" % (doBiasMotionError))
        print("\tBlake-Zisserman on reprojection errors %s" % blakeZisserCam)
        print("\tAcceleration Huber width (sigma): %f" % (huberAccel))
        print("\tGyroscope Huber width (sigma): %f" % (huberGyro))
        print("\tDo time calibration: %s" % (self.__config.estimateParameters['timeOffset']))
        print("\tMax iterations: %d" % (maxIterations))
        print("\tTime offset padding: %f" % (timeOffsetPadding))

        for id, cam in enumerate(self.CameraChain.camList):
            print("\tCamera {} uses rolling shutter model? {}, line delay {} (sec).".format(
                id, cam.isRollingShutter(), cam.getLineDelaySeconds()))
            cam.computeCameraPoses()

        ############################################
        ## initialize camera chain
        ############################################
        #estimate the timeshift for all cameras to the main imu
        if self.__config.estimateParameters['timeOffset']:
            for cam in self.CameraChain.camList:
                cam.findTimeshiftCameraImuPrior(self.ImuList[0], verbose)

        #obtain orientation prior between main imu and camera chain (if no external input provided)
        #and initial estimate for the direction of gravity
        self.CameraChain.findOrientationPriorCameraChainToImu(self.ImuList[0])
        estimatedGravity = self.CameraChain.getEstimatedGravity()
        self.__config.initialGravityEstimate = estimatedGravity

        ############################################
        ## init optimization problem
        ############################################
        #initialize a pose spline using the camera poses in the camera chain
        poseSpline = self.CameraChain.initializePoseSplineFromCameraChain(splineOrder, poseKnotsPerSecond, timeOffsetPadding)
        
        # Initialize bias splines for all IMUs
        for imu in self.ImuList:
            imu.initBiasSplines(poseSpline, biasSplineOrder, biasKnotsPerSecond)
        
        # Now I can build the problem
        problem = inc.CalibrationOptimizationProblem()

        # Initialize all design variables.
        self.initDesignVariables(problem, poseSpline)
        
        ############################################
        ## add error terms
        ############################################
        #Add calibration target reprojection error terms for all camera in chain
        self.CameraChain.addCameraChainErrorTerms(problem, self.poseDv, blakeZissermanDf=blakeZisserCam,
                                                  timeOffsetConstantSparsityPattern=timeOffsetConstantSparsityPattern)
        
        # Initialize IMU error terms.
        for imu in self.ImuList:
            imu.addAccelerometerErrorTerms(problem, self.poseDv, self.gravityExpression, mSigma=huberAccel, accelNoiseScale=accelNoiseScale)
            imu.addGyroscopeErrorTerms(problem, self.poseDv, mSigma=huberGyro, gyroNoiseScale=gyroNoiseScale, g_w=self.gravityExpression)

            # Add the bias motion terms.
            imu.addBiasMotionTerms(problem)
            
        # Add the pose motion terms.
        if doPoseMotionError:
            self.addPoseMotionTerms(problem, mrTranslationVariance, mrRotationVariance)

        self.problem = problem


    def optimize(self, options=None, maxIterations=30, recoverCov=False):

        if options is None:
            options = aopt.Optimizer2Options()
            options.verbose = True
            options.doLevenbergMarquardt = True
            options.levenbergMarquardtLambdaInit = 10.0
            options.nThreads = max(1,multiprocessing.cpu_count()-1)
            options.convergenceDeltaX = 1e-5
            options.convergenceDeltaJ = 1e-2
            options.maxIterations = maxIterations
            options.trustRegionPolicy = aopt.LevenbergMarquardtTrustRegionPolicy(options.levenbergMarquardtLambdaInit)
            options.linearSolver = aopt.BlockCholeskyLinearSystemSolver()

        #run the optimization
        self.optimizer = aopt.Optimizer2(options)
        self.optimizer.setProblem(self.problem)

        optimizationFailed=False
        try: 
            retval = self.optimizer.optimize()
            self.optimizationResult = retval
            if retval.linearSolverFailure:
                optimizationFailed = True
        except:
            optimizationFailed = True

        if optimizationFailed:
            sm.logError("Optimization failed!")
            # raise RuntimeError("Optimization failed!")

        #free some memory
        del self.optimizer
        gc.collect()
        if recoverCov:
            self.recoverCovariance()
        

    def recoverCovariance(self):
        #Covariance ordering (=dv ordering)
        #ORDERING:   N=num cams
        #            1. transformation imu-cam0 --> 6
        #            2. camera time2imu --> 1*numCams (only if enabled)
        print("Recovering covariance is problematic because evaluation of Jacobians of\n"
              "BSpline MotionErrors for IMU biases is not implemented. Despite these\n"
              "void Jacobians, the computation for covariance takes too long!")
        tic = time.time()
        estimator = inc.IncrementalEstimator(CALIBRATION_GROUP_ID)
        rval = estimator.addBatch(self.problem, True)    
        est_stds = np.sqrt(estimator.getSigma2Theta().diagonal())
        toc = time.time()
        elapsed = toc - tic
        print("Covariance recovery takes {} secs".format(elapsed))

        #split and store the variance
        self.std_trafo_ic = np.array(est_stds[0:6])
        self.std_times = np.array(est_stds[6:])

        cam_std_start_index = 0
        for cam_id, cam in enumerate(self.CameraChain.camList):
            num_associated_stds = cam.associateVariableStds(
                est_stds, cam_std_start_index, self.__config.estimateParameters, cam_id)
            cam_std_start_index += num_associated_stds
    
    def saveImuSetParametersYaml(self, resultFile):
        imuSetConfig = kc.ImuSetParameters(resultFile, True)
        for imu in self.ImuList:
            imuConfig = imu.getImuConfig()
            imuConfig.setGravityInTarget(self.gravityDv.toEuclidean())
            imuSetConfig.addImuParameters(imu_parameters=imuConfig)
        imuSetConfig.writeYaml(resultFile)

    def saveCamChainParametersYaml(self, resultFile):    
        chain = self.CameraChain.chainConfig
        nCams = len(self.CameraChain.camList)
    
        # Calibration results
        for camNr in range(0,nCams):
            #cam-cam baselines           
            if camNr > 0:
                T_cB_cA, baseline = self.CameraChain.getResultBaseline(camNr-1, camNr)
                chain.setExtrinsicsLastCamToHere(camNr, T_cB_cA)

            #imu-cam trafos
            T_ci = self.CameraChain.getResultTrafoImuToCam(camNr)
            chain.setExtrinsicsImuToCam(camNr, T_ci)

            if self.__config.estimateParameters['timeOffset']:
                #imu to cam timeshift
                timeshift = float(self.CameraChain.getResultTimeShift(camNr))
                chain.setTimeshiftCamImu(camNr, timeshift)

            if self.__config.estimateParameters['shutter']:
                lineDelay = self.CameraChain.getResultLineDelay(camNr)
                chain.setLineDelay(camNr, int(lineDelay * 1e9))

            if self.__config.estimateParameters['intrinsics']:
                model, coeffs = chain.getIntrinsics(camNr)
                projection = self.CameraChain.getResultProjection(camNr)
                chain.setIntrinsics(camNr, model, projection)

            if self.__config.estimateParameters['distortion']:
                model, coeffs = chain.getDistortion(camNr)
                distortion = self.CameraChain.getResultDistortion(camNr)
                chain.setDistortion(camNr, model, distortion)

        try:
            chain.writeYaml(resultFile)
        except:
            print("ERROR: Could not write parameters to file: {0}\n".format(resultFile))
    
    def computeResidualStatistics(self):
        """
        return
        stats: a dict including JFinal, camera noise std, and IMU noise parameters.
        JFinal is the final cost function value which considers the measurement noise weighting but not the MEstimator weighting.
        """
        noiseParameterNames = ["accelerometer_noise_density", "accelerometer_random_walk", 
            "gyroscope_noise_density", "gyroscope_random_walk"]
        
        stats = dict()
        stats['JFinal'] = self.optimizationResult.JFinal

        for cidx, cam in enumerate(self.CameraChain.camList):
            if len(cam.allReprojectionErrors)>0:
                rawErrors = np.array([rerr.error() for reprojectionErrors in cam.allReprojectionErrors 
                        for rerr in reprojectionErrors])

                camNoise = np.std(rawErrors, 0, ddof=1)
                camName = 'cam{}'.format(cidx)
                stats[camName] = dict()
                stats[camName]['image_noise_std_dev'] = camNoise
                cov = np.matmul(rawErrors.transpose(), rawErrors) / rawErrors.shape[0]
                stats[camName]['image_noise_cov'] = cov
                print("Reprojection error (cam{0}) [px]: mean {1}, median {2}, std: {3}, cov: {4}, #terms: {5}".format(
                        cidx, np.mean(rawErrors, 0), np.median(rawErrors, 0), camNoise, cov, len(rawErrors)))
            else:
                print("Reprojection error (cam{0}) [px]:     no corners".format(cidx))

        for iidx, imu in enumerate(self.ImuList):
            f = imu.imuConfig.getUpdateRate()
            rootf = np.sqrt(f)
            rootdt = 1.0 / rootf

            # compute noise stats
            eGyro = np.array([ e.error() for e in imu.gyroErrors ]) # N x 3
            gyroNoiseDiscrete = np.std(eGyro, 0, ddof=1)
            print("Gyroscope error (imu{0}) [rad/s]: mean {1}, median {2}, std: {3}".format(
                    iidx, np.mean(eGyro, 0), np.median(eGyro, 0), gyroNoiseDiscrete))
            eAccel = np.array([ e.error() for e in imu.accelErrors ])
            accelNoiseDiscrete = np.std(eAccel, 0, ddof=1)
            print("Accelerometer error (imu{0}) [m/s^2]: mean {1}, median {2}, std: {3}".format(
                    iidx, np.mean(eAccel, 0), np.median(eAccel, 0), accelNoiseDiscrete))

            imuName = 'imu{}'.format(iidx)
            stats[imuName] = dict()
            stats[imuName]["accelerometer_noise_density"] = accelNoiseDiscrete * rootdt
            stats[imuName]["gyroscope_noise_density"] = gyroNoiseDiscrete * rootdt

            accelNoiseCovDiscrete = np.matmul(eAccel.transpose(), eAccel) / eAccel.shape[0]
            gyroNoiseCovDiscrete = np.matmul(eGyro.transpose(), eGyro) / eGyro.shape[0]
            stats[imuName]["accelerometer_noise_cov"] = accelNoiseCovDiscrete / f
            stats[imuName]["gyroscope_noise_cov"] = gyroNoiseCovDiscrete / f

            # compute bias random walk stats by sampling the bias splines
            padding = 1.0 # remove padding from both ends to avert ripple effect.
            clampstart = self.poseDv.spline().t_min() + padding
            clampend = self.poseDv.spline().t_max() - padding

            samplingFactor = [0.1, 0.3, 1, 3, 10]
            accWalkList = np.zeros((len(samplingFactor), 4))
            gyroWalkList = np.zeros((len(samplingFactor), 4))
            for fid, factor in enumerate(samplingFactor):
                biasSamplingInterval = factor / f
                imuTimes = np.arange(clampstart, clampend, biasSamplingInterval)

                gyroBiasList = []
                accBiasList = []
                for time in imuTimes:
                    gyro_bias = imu.evaluateGyroBias(time)
                    acc_bias = imu.evaluateAccelerometerBias(time)
                    gyroBiasList.append(gyro_bias)
                    accBiasList.append(acc_bias)
                gyroBiasDiff = np.diff(gyroBiasList, axis=0)
                accBiasDiff = np.diff(accBiasList, axis=0)

                gyroWalkDiscrete = np.std(gyroBiasDiff, 0, ddof=1)
                accWalkDiscrete = np.std(accBiasDiff, 0, ddof=1)

                gyroWalk = gyroWalkDiscrete * rootf
                accWalk = accWalkDiscrete * rootf

                gyroWalkList[fid][0] = biasSamplingInterval
                gyroWalkList[fid][1:] = gyroWalk
                accWalkList[fid][0] = biasSamplingInterval
                accWalkList[fid][1:] = accWalk

            stats[imuName]["accelerometer_random_walk"] = accWalkList
            stats[imuName]["gyroscope_random_walk"] = gyroWalkList
        return stats

    @staticmethod
    def addSplineDesignVariables(problem, dvc, setActive=True, group_id=HELPER_GROUP_ID):
        for i in range(0,dvc.numDesignVariables()):
            dv = dvc.designVariable(i)
            dv.setActive(setActive)
            problem.addDesignVariable(dv, group_id)
