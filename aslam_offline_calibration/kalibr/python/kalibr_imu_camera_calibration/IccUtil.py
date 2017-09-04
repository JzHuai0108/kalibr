from sm import PlotCollection
import IccPlots as plots

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
    timeList = list()
    for obs in cself.CameraChain.camList[0].targetObservations:
        # Build a transformation expression for the time.
        frameTime = cself.CameraChain.camList[0].cameraTimeToImuTimeDv.toExpression() + obs.time().toSec() + cself.CameraChain.camList[0].timeshiftCamToImuPrior
        frameTimeScalar = frameTime.toScalar()
        timeList.append(frameTimeScalar)
    camTimes = np.array(timeList)

    imuTimes = np.array([im.stamp.toSec() + imu.timeOffset for im in imu.imuData \
                      if im.stamp.toSec() + imu.timeOffset > poseSplineDv.spline().t_min() \
                      and im.stamp.toSec() + imu.timeOffset < poseSplineDv.spline().t_max() ])

    refPoseStream = open(filename.replace("bspline", "ref_pose",1), 'w')
    print >> refPoseStream, "%B spline poses: time, T_w_b(txyz, qxyzw) estimated by Kalibr calibrator"  
    saveBsplineRefPose(imuTimes, poseSplineDv, stream=refPoseStream)
    refPoseStream.close()
    refTimeFile = filename.replace("bspline", "ref_camera_time",1)
    np.savetxt(refTimeFile, camTimes.T, fmt='%.9f')

    refImuFile = filename.replace("bspline", "ref_imu_meas",1)
    saveBsplineRefImuMeas(cself, refImuFile) 
  
    modelFile = filename.replace("bspline", "knotCoeffT",1) 
    saveBsplineModel(cself, modelFile)
    
    camNr = 0
    imageCornerPoints = cself.CameraChain.getCornersImageSample(poseSplineDv, 0.0).T
    targetCornerPoints = cself.CameraChain.getCornersTargetSample(camNr).T
    sampleImageCorners = filename.replace("bspline", "sampleImageCorners",1)
    sampleTargetCorners = filename.replace("bspline", "sampleTargetCorners",1)
    np.savetxt(sampleImageCorners,imageCornerPoints, fmt=['%.5f', '%.5f', '%.5f', '%.5f', '%.5f', '%.5f'])
    np.savetxt(sampleTargetCorners,targetCornerPoints, fmt=['%.5f', '%.5f', '%.5f'])
    return modelFile

def loadBspline(filename="knotCoeffT.txt", targetYaml= "april_6x6.yaml"):
    poseSplineDv2 = loadBsplineModel(filename)   
    refTimeFile = filename.replace("knotCoeffT", "ref_camera_time", 1)
    camTimes = np.loadtxt(refTimeFile)
    print >> sys.stdout, "sampling starting time ", '%.9f' % camTimes[0]
    refPoseFile = filename.replace("knotCoeffT", "ref_posex", 1)
    refPoseStream2 = open(refPoseFile, 'w')
    print >> sys.stdout, "got a transformation sample ", poseSplineDv2.spline().transformation(poseSplineDv2.spline().t_min() + 1.0)
    
    chainYaml = filename.replace("knotCoeffT", "camchain", 1)
    chainYaml = chainYaml.replace(".txt", ".yaml", 1)
    print "Reading camera chain:", chainYaml
    chain = kc.CameraChainParameters(chainYaml)  
    chain.printDetails()
    camNr=0
    T_c0_imu = chain.getExtrinsicsImuToCam(camNr)  
    print >> refPoseStream2, "%B spline poses: time, T_w_c(txyz, qxyzw) estimated by Kalibr calibrator"  
    saveBsplineRefPose(camTimes.T, poseSplineDv2, stream=refPoseStream2, T_b_c=T_c0_imu.inverse())
    refPoseStream2.close()

    # simulate some camera observations
    #TODO: generalize the following by referencing 
    #grid = setupCalibrationTarget(self, targetConfig, showExtraction=False, showReproj=False, imageStepping=False):
    
    #load calibration target configuration     
    targetConfig = kc.CalibrationTargetParameters(targetYaml)
    targetParams = targetConfig.getTargetParams()
    print "Initializing calibration target:"
    targetConfig.printDetails()
    options = acv_april.AprilgridOptions()
    options.showExtractionVideo = False
    options.minTagsForValidObs = int( np.max( [targetParams['tagRows'], targetParams['tagCols']] ) + 1 )
            
    grid = acv_april.GridCalibrationTargetAprilgrid(targetParams['tagRows'],
                                                            targetParams['tagCols'], 
                                                            targetParams['tagSize'], 
                                                            targetParams['tagSpacing'], 
                                                            options)
    simulatedObs = acv.GridCalibrationTargetObservation(grid) #this line should work though 
    # in the constructor GridCalibrationTargetObservation(GridCalibrationTarget::Ptr) is
    # inconsistent with GridCalibrationTargetBase::Ptr
    
    camConfig = chain.getCameraParameters(camNr)
    camera = kc.AslamCamera.fromParameters( camConfig)

    # get the GS camera observations
    sm_T_w_c= getCameraPoseAt(camTimes[0], poseSplineDv2, T_b_c=T_c0_imu.inverse())
    print >> sys.stdout, '%.9f' % camTimes[0], ' '.join(map(str,sm_T_w_c.t())), ' '.join(map(str,sm_T_w_c.q()))
    simulatedObs.set_T_t_c(sm_T_w_c)
    v2Point = np.array([[0],[0]])
    imageCornerProjected= list()
    print "total corners ", simulatedObs.getTotalTargetPoint()
    print "image rows ", simulatedObs.imRows()
    
    lineDelay = 20e-3/480
    for iota in range(0, simulatedObs.getTotalTargetPoint(), 1):
        # get the initial observation  
        lastImagePoint = simulatedObs.projectATargetPoint(camera.geometry, sm_T_w_c, iota) #3x1, this is GS camera model
        print "image point ", lastImagePoint.T 
        numIter = 0       
        while numIter < 5:
            currTime = lastImagePoint[1, 0]*lineDelay + camTimes[0]
            print "currtime ", '%.9f' % currTime
            sm_T_w_cx= getCameraPoseAt(currTime, poseSplineDv2, T_b_c=T_c0_imu.inverse())
            imagePoint = simulatedObs.projectATargetPoint(camera.geometry, sm_T_w_cx, iota)
            print "image point ", imagePoint.T
            
            if(np.linalg.norm(lastImagePoint - imagePoint) < 1e-3):
                break
            numIter += 1
            lastImagePoint = imagePoint
        print 
    	imageCornerProjected.append(imagePoint)       
        
    imageCornerProjectedArray = np.array(imageCornerProjected)
    reproducedImageCornerFile = filename.replace("knotCoeffT", "regenerated_corners", 1)
    print "reproducedImageCorner shape ", imageCornerProjectedArray.shape
    np.savetxt(reproducedImageCornerFile, imageCornerProjectedArray, fmt=['%.9f', '%.9f', '%d'])

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
        print >> stream, '%.9f' % timeScalar, ' '.join(map(str,sm_T_w_c.t())), ' '.join(map(str,sm_T_w_c.q()))
      
def getCameraPoseAt(timeScalar, poseSplineDv, T_b_c=sm.Transformation()):
    timeOffsetPadding = 0.0  
      
    dv = aopt.Scalar(timeScalar)
    timeExpression = dv.toExpression()
        
    if timeScalar <= poseSplineDv.spline().t_min() or timeScalar >= poseSplineDv.spline().t_max():
        print >> sys.stdout, "Warn: time out of range "
        return sm.Transformation()
       
    T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
    sm_T_w_c = sm.Transformation(T_w_b.toTransformationMatrix())*T_b_c  
    return sm_T_w_c

def saveBsplineRefImuMeas(cself, filename):
    print >> sys.stdout, "saving Bspline IMU measurements time (axyz, wxyz in metric units) estimated by Kalibr calibrator"
    print >> sys.stdout, "to ", filename

    idx = 0
    imu = cself.ImuList[idx]    
    poseSplineDv = cself.poseDv
    times = np.array([[im.stamp.toSec() + imu.timeOffset for im in imu.imuData \
                      if im.stamp.toSec() + imu.timeOffset > poseSplineDv.spline().t_min() \
                      and im.stamp.toSec() + imu.timeOffset < poseSplineDv.spline().t_max() ]])

    predictedAng_body =  np.array([err.getPredictedMeasurement() for err in imu.gyroErrors])    
    predicetedAccel_body =  np.array([err.getPredictedMeasurement() for err in imu.accelErrors])
    print >> sys.stdout, "times ", times.shape
    print >> sys.stdout, "ang ", predictedAng_body.shape
    print >> sys.stdout, "accel ", predicetedAccel_body.shape    
    whole=np.concatenate((times.T, predicetedAccel_body, predictedAng_body),axis=1)
    np.savetxt(filename,whole, fmt=['%.9f', '%.5f', '%.5f', '%.5f', '%.5f', '%.5f', '%.5f'])

def saveBsplineModel(cself, filename):    
    poseSpline = cself.poseDv.spline()    
    print >> sys.stdout,"splineOrder", poseSpline.splineOrder(), "should be 6"    
    poseSpline.savePoseSplineToFile(filename)   
    print >> sys.stdout, "saved B spline model to ", filename

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
