import aslam_backend as aopt
import aslam_splines as asp
import IccUtil as util
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

def addSplineDesignVariables(problem, dvc, setActive=True, group_id=HELPER_GROUP_ID):
    for i in range(0,dvc.numDesignVariables()):
        dv = dvc.designVariable(i)
        dv.setActive(setActive)
        problem.addDesignVariable(dv, group_id)

class IccCalibrator(object):
    def __init__(self, config):
        self.ImuList = []
        self.__config = config
        self.CameraChain = None
        self.poseDv = None
        self.gravityDv = None
        self.gravityExpression = None
        self.noTimeCalibration = False
        self.problem = None
        self.optimizer = None

    def initDesignVariables(self, problem, poseSpline):
        # Initialize the system pose spline (always attached to imu0) 
        self.poseDv = asp.BSplinePoseDesignVariable( poseSpline )
        addSplineDesignVariables(problem, self.poseDv)

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
                      verbose=False  ):

        print "\tSpline order: %d" % (splineOrder)
        print "\tPose knots per second: %d" % (poseKnotsPerSecond)
        print "\tDo pose motion regularization: %s" % (doPoseMotionError)
        print "\t\txddot translation variance: %f" % (mrTranslationVariance)
        print "\t\txddot rotation variance: %f" % (mrRotationVariance)
        print "\tBias knots per second: %d" % (biasKnotsPerSecond)
        print "\tDo bias motion regularization: %s" % (doBiasMotionError)
        print "\tBlake-Zisserman on reprojection errors %s" % blakeZisserCam
        print "\tAcceleration Huber width (sigma): %f" % (huberAccel)
        print "\tGyroscope Huber width (sigma): %f" % (huberGyro)
        print "\tDo time calibration: %s" % (self.__config.estimateParameters['timeOffset'])
        print "\tMax iterations: %d" % (maxIterations)
        print "\tTime offset padding: %f" % (timeOffsetPadding)

        for cam in self.CameraChain.camList:
            cam.generateIntrinsicsInitialGuess(self.__config.estimateParameters['intrinsics'],
                                               self.__config.estimateParameters['distortion'],
                                               self.__config.estimateParameters['shutter'])
            cam.computeCameraPoses()

        ############################################
        ## initialize camera chain
        ############################################
        #estimate the timeshift for all cameras to the main imu
        self.noTimeCalibration = not self.__config.estimateParameters['timeOffset']
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
            imu.initBiasSplines(poseSpline, splineOrder, biasKnotsPerSecond)
        
        # Now I can build the problem
        problem = inc.CalibrationOptimizationProblem()

        # Initialize all design variables.
        self.initDesignVariables(problem, poseSpline)
        
        ############################################
        ## add error terms
        ############################################
        #Add calibration target reprojection error terms for all camera in chain
        self.CameraChain.addCameraChainErrorTerms(problem, self.poseDv, blakeZissermanDf=blakeZisserCam, timeOffsetPadding=timeOffsetPadding)
        
        # Initialize IMU error terms.
        for imu in self.ImuList:
            imu.addAccelerometerErrorTerms(problem, self.poseDv, self.gravityExpression, mSigma=huberAccel, accelNoiseScale=accelNoiseScale)
            imu.addGyroscopeErrorTerms(problem, self.poseDv, mSigma=huberGyro, gyroNoiseScale=gyroNoiseScale, g_w=self.gravityExpression)

            # Add the bias motion terms.
            if doBiasMotionError:
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
            if retval.linearSolverFailure:
                optimizationFailed = True
        except:
            optimizationFailed = True

        if optimizationFailed:
            sm.logError("Optimization failed!")
            raise RuntimeError("Optimization failed!")
        
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
        
        print "Recovering covariance..."
        estimator = inc.IncrementalEstimator(CALIBRATION_GROUP_ID)
        rval = estimator.addBatch(self.problem, True)    
        est_stds = np.sqrt(estimator.getSigma2Theta().diagonal())
        
        #split and store the variance
        self.std_trafo_ic = np.array(est_stds[0:6])
        self.std_times = np.array(est_stds[6:])
    
    def saveImuSetParametersYaml(self, resultFile):
        imuSetConfig = kc.ImuSetParameters(resultFile, True)
        for imu in self.ImuList:
            imuConfig = imu.getImuConfig()
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
            print "ERROR: Could not write parameters to file: {0}\n".format(resultFile)
            
