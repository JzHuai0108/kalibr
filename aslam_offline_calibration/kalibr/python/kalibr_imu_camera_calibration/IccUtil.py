from sm import PlotCollection
import IccPlots as plots
import matplotlib.pyplot as plt #plt.hist
import aslam_backend as aopt
import aslam_cv as acv
import aslam_cameras_april as acv_april
import sm
import bsplines
import aslam_splines as asp
import kalibr_common as kc

import numpy as np
import pylab as pl
import sys
import subprocess
import yaml
import time
from matplotlib.backends.backend_pdf import PdfPages
import StringIO
import matplotlib.patches as patches
from random import gauss # simulate reprojection errors

def printErrorStatistics(cself, dest=sys.stdout):
    # Reprojection errors
    print >> dest, "Normalized Residuals\n----------------------------"
    for cidx, cam in enumerate(cself.CameraChain.camList):
        if len(cam.allReprojectionErrors)>0:
            e2 = np.array([ np.sqrt(rerr.evaluateError()) for reprojectionErrors in cam.allReprojectionErrors for rerr in reprojectionErrors])
            print >> dest, "Reprojection error (cam{0}):     mean {1}, median {2}, std: {3}".format(cidx, np.mean(e2), np.median(e2), np.std(e2) )
        else:
            print >> dest, "Reprojection error (cam{0}):     no corners".format(cidx)
    
    for iidx, imu in enumerate(cself.ImuList):
        # Gyro errors
        e2 = np.array([ np.sqrt(e.evaluateError()) for e in imu.gyroErrors ])
        print >> dest, "Gyroscope error (imu{0}):        mean {1}, median {2}, std: {3}".format(iidx, np.mean(e2), np.median(e2), np.std(e2))
        # Accelerometer errors
        e2 = np.array([ np.sqrt(e.evaluateError()) for e in imu.accelErrors ])
        print >> dest, "Accelerometer error (imu{0}):    mean {1}, median {2}, std: {3}".format(iidx, np.mean(e2), np.median(e2), np.std(e2))

    print >> dest, ""
    print >> dest, "Residuals\n----------------------------"
    for cidx, cam in enumerate(cself.CameraChain.camList):
        if len(cam.allReprojectionErrors)>0:
            e2 = np.array([ np.linalg.norm(rerr.error()) for reprojectionErrors in cam.allReprojectionErrors for rerr in reprojectionErrors])
            print >> dest, "Reprojection error (cam{0}) [px]:     mean {1}, median {2}, std: {3}".format(cidx, np.mean(e2), np.median(e2), np.std(e2) )
        else:
            print >> dest, "Reprojection error (cam{0}) [px]:     no corners".format(cidx)
    
    for iidx, imu in enumerate(cself.ImuList):
        # Gyro errors
        e2 = np.array([ np.linalg.norm(e.error()) for e in imu.gyroErrors ])
        print >> dest, "Gyroscope error (imu{0}) [rad/s]:     mean {1}, median {2}, std: {3}".format(iidx, np.mean(e2), np.median(e2), np.std(e2))
        # Accelerometer errors
        e2 = np.array([ np.linalg.norm(e.error()) for e in imu.accelErrors ])
        print >> dest, "Accelerometer error (imu{0}) [m/s^2]: mean {1}, median {2}, std: {3}".format(iidx, np.mean(e2), np.median(e2), np.std(e2))

def printGravity(cself):
    print
    print "Gravity vector: (in target coordinates): [m/s^2]"
    print cself.gravityDv.toEuclidean()

def printResults(cself, withCov=False):
    nCams = len(cself.CameraChain.camList)
    for camNr in range(0,nCams):
        T_cam_b = cself.CameraChain.getResultTrafoImuToCam(camNr)

        print
        print "Transformation T_cam{0}_imu0 (imu0 to cam{0}, T_ci): ".format(camNr)
        if withCov and camNr==0:
            print "\t quaternion: ", T_cam_b.q(), " +- ", cself.std_trafo_ic[0:3]
            print "\t translation: ", T_cam_b.t(), " +- ", cself.std_trafo_ic[3:]
        print T_cam_b.T()
        
        if not cself.noTimeCalibration:
            print
            print "cam{0} to imu0 time: [s] (t_imu = t_cam + shift)".format(camNr)
            print cself.CameraChain.getResultTimeShift(camNr),
            
            if withCov:
                print " +- ", cself.std_times[camNr]
            else:
                print

    print
    for (imuNr, imu) in enumerate(cself.ImuList):
        print "IMU{0}:\n".format(imuNr), "----------------------------"
        imu.getImuConfig().printDetails()
            
def printBaselines(self):
    #print all baselines in the camera chain
    if nCams > 1:
        for camNr in range(0,nCams-1):
            T, baseline = cself.CameraChain.getResultBaseline(camNr, camNr+1)
            
            if cself.CameraChain.camList[camNr+1].T_extrinsic_fixed:
                isFixed = "(fixed to external data)"
            else:
                isFixed = ""
            
            print
            print "Baseline (cam{0} to cam{1}): [m] {2}".format(camNr, camNr+1, isFixed)
            print T.T()
            print baseline, "[m]"
    
   
def generateReport(cself, filename="report.pdf", showOnScreen=True):  
    figs = list()
    plotter = PlotCollection.PlotCollection("Calibration report")
    offset = 3010
    
    #Output calibration results in text form.
    sstream = StringIO.StringIO()
    printResultTxt(cself, sstream)
    
    text = [line for line in StringIO.StringIO(sstream.getvalue())]
    linesPerPage = 40
    
    while True:
        fig = pl.figure(offset)
        offset += 1

        left, width = .05, 1.
        bottom, height = -.05, 1.
        right = left + width
        top = bottom + height
        
        ax = fig.add_axes([.0, .0, 1., 1.])
        # axes coordinates are 0,0 is bottom left and 1,1 is upper right
        p = patches.Rectangle((left, bottom), width, height, fill=False, transform=ax.transAxes, \
                                 clip_on=False, edgecolor="none")
        ax.add_patch(p)
        pl.axis('off')

        printText = lambda t: ax.text(left, top, t, fontsize=8, \
                                     horizontalalignment='left', verticalalignment='top',\
                                     transform=ax.transAxes)
        
        if len(text) > linesPerPage:
            printText("".join(text[0:linesPerPage]))
            figs.append(fig)
            text = text[linesPerPage:]
        else:
            printText("".join(text[0:]))
            figs.append(fig)
            break
    
    #plot imu stuff (if we have imus)
    for iidx, imu in enumerate(cself.ImuList):
        f = pl.figure(offset+iidx)
        plots.plotAccelerations(cself, iidx, fno=f.number, noShow=True)
        plotter.add_figure("imu{0}: accelerations".format(iidx), f)
        figs.append(f)
        offset += len(cself.ImuList)

        f = pl.figure(offset+iidx)
        plots.plotAccelErrorPerAxis(cself, iidx, fno=f.number, noShow=True)
        plotter.add_figure("imu{0}: acceleration error".format(iidx), f)
        figs.append(f)
        offset += len(cself.ImuList)

        f = pl.figure(offset+iidx)
        plots.plotAccelBias(cself, iidx, fno=f.number, noShow=True)
        plotter.add_figure("imu{0}: accelerometer bias".format(iidx), f)
        figs.append(f)
        offset += len(cself.ImuList)

        f = pl.figure(offset+iidx)
        plots.plotAngularVelocities(cself, iidx, fno=f.number, noShow=True)
        plotter.add_figure("imu{0}: angular velocities".format(iidx), f)
        figs.append(f)
        offset += len(cself.ImuList)

        f = pl.figure(offset+iidx)
        plots.plotGyroErrorPerAxis(cself, iidx, fno=f.number, noShow=True)
        plotter.add_figure("imu{0}: angular velocity error".format(iidx), f)
        figs.append(f)
        offset += len(cself.ImuList)

        f = pl.figure(offset+iidx)
        plots.plotAngularVelocityBias(cself, iidx, fno=f.number, noShow=True)
        plotter.add_figure("imu{0}: gyroscope bias".format(iidx), f)
        figs.append(f)
        offset += len(cself.ImuList)

    #plot cam stuff
    if cself.CameraChain:        
        for cidx, cam in enumerate(cself.CameraChain.camList):
            f = pl.figure(offset+cidx)
            title="cam{0}: reprojection errors".format(cidx);
            plots.plotReprojectionScatter(cself, cidx, fno=f.number, noShow=True, title=title)
            plotter.add_figure(title, f)
            figs.append(f)
            offset += len(cself.CameraChain.camList)

    #write to pdf
    pdf=PdfPages(filename)
    for fig in figs:
        pdf.savefig(fig)
    pdf.close()

    if showOnScreen:
        plotter.show()

def saveBspline(cself, filename="bspline.txt"):    
    
    idx = 0
    imu = cself.ImuList[idx]    
    poseSplineDv = cself.poseDv

    imuTimes = np.array([im.stamp.toSec() + imu.timeOffset for im in imu.imuData \
                      if im.stamp.toSec() + imu.timeOffset > poseSplineDv.spline().t_min() \
                      and im.stamp.toSec() + imu.timeOffset < poseSplineDv.spline().t_max() ])

    timeList = list()
    for obs in cself.CameraChain.camList[0].targetObservations:
        # Build a transformation expression for the time.
        frameTime = cself.CameraChain.camList[0].cameraTimeToImuTimeDv.toExpression() + obs.time().toSec() + cself.CameraChain.camList[0].timeshiftCamToImuPrior
        frameTimeScalar = frameTime.toScalar()
        if frameTimeScalar > imuTimes[0] and frameTimeScalar < imuTimes[-1]:
            timeList.append(frameTimeScalar)
    refStateTimes = np.array(timeList)

 
    refPoseStream = open(filename.replace("bspline", "ref_pose",1), 'w')
    print >> refPoseStream, "%poses generated at the IMU rate from the B-spline: time, T_w_b(txyz, qxyzw)"  
    saveBsplineRefPose(imuTimes, poseSplineDv, stream=refPoseStream)
    refPoseStream.close()
    
    refStateFile = filename.replace("bspline", "ref_state",1)
    refStateStream = open(refStateFile, 'w')
    print >> sys.stdout, "  Saving system states at camera rate generated from B-spline to", refStateFile   
    saveBsplineRefStates(refStateTimes, poseSplineDv, cself, refStateStream)
    refStateStream.close()
   
    refImuFile = filename.replace("bspline", "ref_imu_meas",1)
    saveBsplineRefImuMeas(cself, refImuFile) 
  
    modelFile = filename.replace("bspline", "knotCoeffT",1) 
    saveBsplineModel(cself, modelFile)
    
    if 0:
        camNr = 0
        imageCornerPoints = cself.CameraChain.getCornersImageSample(poseSplineDv, 0.0).T
        targetCornerPoints = cself.CameraChain.getCornersTargetSample(camNr).T
        sampleImageCorners = filename.replace("bspline", "sampleImageCorners",1)
        sampleTargetCorners = filename.replace("bspline", "sampleTargetCorners",1)
        np.savetxt(sampleImageCorners,imageCornerPoints, fmt=['%.5f', '%.5f', '%.5f', '%.5f', '%.5f', '%.5f'])
        np.savetxt(sampleTargetCorners,targetCornerPoints, fmt=['%.5f', '%.5f', '%.5f'])
    return modelFile

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

class RsVisualInertialMeasViaBsplineSimulator(object):
    '''simulate visual(rolling shutter) inertial measurements with a predefined BSpline model representing a realistic motion'''
    def __init__(self, filename="knotCoeffT.txt", targetYaml= "april_6x6.yaml", \
      chainYaml='camchain.yaml', time_offset=0.0, time_readout=30e-3):
        self.modelFilename = filename
        self.simulPoseSplineDv = loadBsplineModel(filename)        

        
        refStateFile = filename.replace("knotCoeffT", "ref_state", 1)
        fullRefStates = loadArrayWithHeader(refStateFile)

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
            print >> trimmedStateStream, '%d %s %.9f' % (iota-startIndex, iota-startIndex, fullRefStates[iota, 2]), \
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
        print "which target is used in the simulation?"
        targetConfig.printDetails()
        self.setupCalibrationTarget(targetConfig, showExtraction=False, showReproj=False, imageStepping=False)
        self.imageHeight = resolution[1]  
        self.timeOffset = time_offset 
        self.lineDelay = time_readout/self.imageHeight
        print "line delay ", self.lineDelay
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
        print "all target corners shape", self.allTargetCorners.shape   
        
        assert self.allTargetCorners.shape[0] == self.simulatedObs.getTotalTargetPoint()
   
    
    def simulateOneFrame(self):
        for frameId in range(1):
             state_time = self.refStateTimes[0]
             line_delay = self.lineDelay             
             imageCornerProjectedArray = self.naiveMethodToRsProjection(state_time, line_delay, True)
             imageCornerProjectedArray2, unusedKeypoints, unusedImageOffset = self.newtonMethodToRsProjection(state_time, line_delay, 1.0, 0.0, True)        

        imageCornerOffsetNorms = list()
        bins = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0]
        imageCornerOffsetNorms += unusedImageOffset
        counts, newBins, patches = plt.hist(imageCornerOffsetNorms, bins)
        print '  counts', counts, '\n  bins\t', newBins
        plt.title('Distribution of offset infinity norm due to line delay and noise')
        plt.show()
       
        reproducedImageCornerFile = self.modelFilename.replace("knotCoeffT", "sampleImageCornersByNaiveAndNewton", 1)
        print "  ReproducedImageCorner shape ", imageCornerProjectedArray.shape, ' should be (144, 3, 1)'
        np.savetxt(reproducedImageCornerFile, np.concatenate((imageCornerProjectedArray[:,:,0], imageCornerProjectedArray2[:,:,0]), axis=1), fmt=['%.9f', '%.9f', '%d', '%.9f', '%.9f', '%d'])

    def simulate(self, landmarkOutputFilename, reprojectionSigma):
        '''simulate camera observations for frames at all ref state times and plus noise''' 
        print 'Reprojection Gaussian noise sigma %.5f' % reprojectionSigma
        print 'Image time advance relative to Imu clock', self.timeOffset
        
        imageCornerOffsetNorms = list()
        bins = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0]
        for frameId in range(0, self.refStateTimes.shape[0], 1):
             state_time = self.refStateTimes[frameId]
             line_delay = self.lineDelay
             unusedCornerProjections, frameKeypoints, imageCornerOffset = self.newtonMethodToRsProjection(state_time, line_delay, reprojectionSigma, self.timeOffset)
             imageCornerOffsetNorms += imageCornerOffset
             print '  Projected %4d' % len(frameKeypoints), 'target corners for state at %.9f' % state_time 
             for lmObs in frameKeypoints: #each lmObs (lmId, keypoint Id, x, y, size)
                 self.landmark2ObservationList[lmObs[0]].append((frameId, lmObs[1], lmObs[2], lmObs[3], lmObs[4]))
        lmStream = open(landmarkOutputFilename, 'w')
        lmStream.write('%each line contains a landmark: lm id, lm p_w(homogeneous xyzw), quality, num obs, <frameid> <keypoint id, x, y, size>, ... , <frameid> <keypoint id, x, y, size>\r\n')
        probe = 0
        homogeneousW = 1.0
        positionQuality = 1.0
        for lm, obsList in sorted(self.landmark2ObservationList.iteritems()):
            assert probe == lm            
            lmStream.write('%d ' % lm + ' '.join(map(str, self.allTargetCorners[lm,:]))+ \
              ' %.3f %.3f %d' % (homogeneousW, positionQuality, len(obsList)))
            for observation in obsList:
                lmStream.write(' '+' '.join(map(str, observation)))
            lmStream.write('\r\n')
            probe += 1
        lmStream.close()
        print '  Written landmark observations to', landmarkOutputFilename
        print '  To use the simulated dataset, you need to pull out the gravity vector from report-...-dynamic.pdf and plug it into the customized bundle adjustment calibrator'
        print '  Histogram of the infinity norm the offset due to line delay and noise'
        counts, newBins, patches = plt.hist(imageCornerOffsetNorms, bins)
        print '  counts', counts, '\n  bins\t', newBins        
        plt.title('Distribution of offset infinity norm due to line delay and noise')        
        plt.show()
        print

    def naiveMethodToRsProjection(self, state_time, line_delay, verbose=False):   
        '''this method is not proved theoretically to converge, but it performs as precise as Newton's method empirically, though slower'''
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
        imageCornerProjected = list() #image keypoints free of noise effect
        imageCornerProjectedOffset = list() # the pixel offset due to line delay and noise effect
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
            if verbose:
                print
            if not aborted:
                #TODO: randomly mask some landmark observations
                imageCornerProjected.append(imagePoint0)                
                xnoise = gauss(0.0, reprojectionSigma)
                ynoise = gauss(0.0, reprojectionSigma)
                frameKeypoints.append((iota, kpId, imagePoint0[0, 0] + xnoise, imagePoint0[1, 0] + ynoise, 12))
                imageCornerProjectedOffset.append(max(np.abs(initialImagePoint[0, 0] - imagePoint0[0, 0] - xnoise), \
                  np.abs(initialImagePoint[1, 0] - imagePoint0[1, 0] - ynoise)))
                kpId += 1        
        assert kpId == len(imageCornerProjected)
        return np.array(imageCornerProjected), frameKeypoints, imageCornerProjectedOffset
    

def saveBsplineRefPose(times, poseSplineDv, stream=sys.stdout, T_b_c=sm.Transformation()):
    
    timeOffsetPadding = 0.0  
    for timeScalar in times:    
        dv = aopt.Scalar(timeScalar)
        timeExpression = dv.toExpression()
        
        if timeScalar <= poseSplineDv.spline().t_min() or timeScalar >= poseSplineDv.spline().t_max():
            print >> sys.stdout, "Warn: time out of range "
            continue       
        T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
        sm_T_w_c = sm.Transformation(T_w_b.toTransformationMatrix())*T_b_c
        # quatInv used here to convert kalibr's JPL quaternion to Halmilton quaternion
        print >> stream, '%.9f' % timeScalar, ' '.join(map(str,sm_T_w_c.t())), ' '.join(map(str, sm.quatInv(sm_T_w_c.q())))
      
def saveBsplineRefStates(times, poseSplineDv, cself, stream=sys.stdout, T_b_c=sm.Transformation()):
    '''save the system states at the camera reference times in a customized format. system states are (T_w_b(txyz, qxyzw), v_w, bg, ba)'''
    
    stream.write('%Each line contains the state at a frame\'s reference time: all in metric units unless specified otherwise\n')
    stream.write('%frame id, frame id in source, timestamp, T_WB[xyz, qxyzw], v_W, bg, ba, isKeyFrame[0 or 1]\n')

    idx = 0
    imu = cself.ImuList[idx]    
   
    gyroBias = imu.gyroBiasDv.spline()   
    accBias = imu.accelBiasDv.spline()

    print '\t\t\tstart time\t\tfinish time'
    print 'poseSpline\t%.9f\t%.9f' % (poseSplineDv.spline().t_min(), poseSplineDv.spline().t_max())
    print 'gyroBias\t%.9f\t%.9f' % (gyroBias.t_min(), gyroBias.t_max())
    print 'accBias\t\t%.9f\t%.9f' % (accBias.t_min(), accBias.t_max())
    print 'imu.timeOffset\t%.9f' % imu.timeOffset
        
    timeOffsetPadding = 0.0
    frameId = 0  
    for timeScalar in times:    
        dv = aopt.Scalar(timeScalar)
        timeExpression = dv.toExpression()        
        if timeScalar <= poseSplineDv.spline().t_min() or timeScalar >= poseSplineDv.spline().t_max() or \
           timeScalar <= gyroBias.t_min() or timeScalar >= gyroBias.t_max() or \
           timeScalar <= accBias.t_min() or timeScalar >= accBias.t_max():
            print "Warn: time out of range in generating ref states"
            continue
        T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
        sm_T_w_b = sm.Transformation(T_w_b.toTransformationMatrix())
        v_w = poseSplineDv.linearVelocity(timeScalar).toEuclidean()
        
        gyro_bias = gyroBias.evalD(timeScalar,0)   
        acc_bias = accBias.evalD(timeScalar,0)
        # quatInv used here to convert kalibr's JPL quaternion to Halmilton quaternion
        print >> stream, frameId, frameId, '%.9f' % timeScalar, ' '.join(map(str,sm_T_w_b.t())), \
        ' '.join(map(str, sm.quatInv(sm_T_w_b.q()))), ' '.join(map(str,v_w)), \
        ' '.join(map(str, gyro_bias)), ' '.join(map(str, acc_bias)), 1
        frameId += 1

def loadArrayWithHeader(arrayFile):
    with open(arrayFile) as f:
        lines = (line for line in f if not (line.startswith('#') or line.startswith('%')))        
        return np.loadtxt(lines, delimiter=' ', skiprows=0)

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

def saveBsplineRefImuMeas(cself, filename):
    print >> sys.stdout, "  Saving IMU measurements generated from B-spline (time wxyz, axyz bg, ba in metric units) to", filename

    idx = 0
    imu = cself.ImuList[idx]    
    poseSplineDv = cself.poseDv
    times = np.array([im.stamp.toSec() + imu.timeOffset for im in imu.imuData \
                      if im.stamp.toSec() + imu.timeOffset > poseSplineDv.spline().t_min() \
                      and im.stamp.toSec() + imu.timeOffset < poseSplineDv.spline().t_max() ])

    predictedAng_body =  np.array([err.getPredictedMeasurement() for err in imu.gyroErrors])    
    predicetedAccel_body =  np.array([err.getPredictedMeasurement() for err in imu.accelErrors])

    gyroBias = imu.gyroBiasDv.spline()    
    gyro_bias_spline = np.array([gyroBias.evalD(t,0) for t in times])
    
    accBias = imu.accelBiasDv.spline()    
    acc_bias_spline = np.array([accBias.evalD(t,0) for t in times])

    print >> sys.stdout, 'Epitome of predicted inertial measurements'
    print >> sys.stdout, "\t#times", times.shape
    print >> sys.stdout, "\t#gyro", predictedAng_body.shape
    print >> sys.stdout, "\t#accel", predicetedAccel_body.shape
    print >> sys.stdout, "\t#gyro bias", gyro_bias_spline.shape
    print >> sys.stdout, "\t#accel bias", acc_bias_spline.shape
    
    whole=np.concatenate((np.array([times]).T, predictedAng_body, predicetedAccel_body, gyro_bias_spline, acc_bias_spline),axis=1)
    np.savetxt(filename,whole, fmt=['%.9f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f'])

def saveBsplineModel(cself, filename):    
    poseSpline = cself.poseDv.spline()    
    print >> sys.stdout,"  splineOrder", poseSpline.splineOrder(), "should be 6"    
    poseSpline.savePoseSplineToFile(filename)   
    print >> sys.stdout, "  saved B spline model to ", filename

def loadBsplineModel(knotCoeffFile):
    splineOrder = 6
    poseSpline = bsplines.BSplinePose(splineOrder, sm.RotationVector() )
    poseSpline.initPoseSplineFromFile(knotCoeffFile)
    print "Initialized a pose spline with ", poseSpline.knots().size, "knots and coeff shape", poseSpline.coefficients().shape
    poseDv = asp.BSplinePoseDesignVariable( poseSpline )
    return poseDv

def saveResultTxt(cself, filename='cam_imu_result.txt'):
    f = open(filename, 'w')
    printResultTxt(cself, stream=f)

def printResultTxt(cself, stream=sys.stdout):
    
    print >> stream, "Calibration results"
    print >> stream, "==================="   
    printErrorStatistics(cself, stream)
  
    # Calibration results
    nCams = len(cself.CameraChain.camList)
    for camNr in range(0,nCams):
        T = cself.CameraChain.getResultTrafoImuToCam(camNr)
        print >> stream, ""
        print >> stream, "Transformation (cam{0}):".format(camNr)
        print >> stream, "-----------------------"
        print >> stream, "T_ci:  (imu0 to cam{0}): ".format(camNr)   
        print >> stream, T.T()
        print >> stream, ""
        print >> stream, "T_ic:  (cam{0} to imu0): ".format(camNr)   
        print >> stream, T.inverse().T()
    
        # Time
        print >> stream, ""
        print >> stream, "timeshift cam{0} to imu0: [s] (t_imu = t_cam + shift)".format(camNr)
        print >> stream, cself.CameraChain.getResultTimeShift(camNr)
        print >> stream, ""

    #print all baselines in the camera chain
    if nCams > 1:
        print >> stream, "Baselines:"
        print >> stream, "----------"

        for camNr in range(0,nCams-1):
            T, baseline = cself.CameraChain.getResultBaseline(camNr, camNr+1)
            print >> stream, "Baseline (cam{0} to cam{1}): ".format(camNr, camNr+1)
            print >> stream, T.T()
            print >> stream, "baseline norm: ", baseline,  "[m]"
            print >> stream, ""
    
    # Gravity
    print >> stream, ""
    print >> stream, "Gravity vector in target coords: [m/s^2]"
    print >> stream, cself.gravityDv.toEuclidean()
    
    print >> stream, ""
    print >> stream, ""
    print >> stream, "Calibration configuration"
    print >> stream, "========================="
    print >> stream, ""

    for camNr, cam in enumerate( cself.CameraChain.camList ):
        print >> stream, "cam{0}".format(camNr)
        print >> stream, "-----"
        cam.camConfig.printDetails(stream)
        cam.targetConfig.printDetails(stream)
        print >> stream, ""
    
	print >> stream, ""
    print >> stream, ""
    print >> stream, "IMU configuration"
    print >> stream, "================="
    print >> stream, ""
    for (imuNr, imu) in enumerate(cself.ImuList):
        print >> stream, "IMU{0}:\n".format(imuNr), "----------------------------"
        imu.getImuConfig().printDetails(stream)
        print >> stream, ""
