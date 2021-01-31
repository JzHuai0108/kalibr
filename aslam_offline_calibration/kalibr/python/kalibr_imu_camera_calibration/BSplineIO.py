
import sys

import numpy as np

import aslam_backend as aopt
import aslam_splines as asp
import bsplines
import sm


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
    poseSpline.savePoseSplineToFile(filename)   
    print("  saved B spline model of order {} to {}".format(poseSpline.splineOrder(), filename))


def saveBspline(cself, bagtag):
    idx = 0
    imu = cself.ImuList[idx]    
    poseSplineDv = cself.poseDv

    # imuTimes = np.array([im.stamp.toSec() + imu.timeOffset for im in imu.imuData \
    #                   if im.stamp.toSec() + imu.timeOffset > poseSplineDv.spline().t_min() \
    #                   and im.stamp.toSec() + imu.timeOffset < poseSplineDv.spline().t_max()])

    # refPoseStream = open("bspline_poses.txt", 'w')
    # print >> refPoseStream, "%poses generated at the IMU rate from the B-spline: time, T_w_b(txyz, qxyzw)"  
    # saveBsplineRefPose(imuTimes, poseSplineDv, stream=refPoseStream)
    # refPoseStream.close()

    # timeList = list()
    # for obs in cself.CameraChain.camList[0].targetObservations:
    #     # Build a transformation expression for the time.
    #     frameTime = cself.CameraChain.camList[0].cameraTimeToImuTimeDv.toExpression() + obs.time().toSec() + cself.CameraChain.camList[0].timeshiftCamToImuPrior
    #     frameTimeScalar = frameTime.toScalar()
    #     if frameTimeScalar > imuTimes[0] and frameTimeScalar < imuTimes[-1]:
    #         timeList.append(frameTimeScalar)
    # refStateTimes = np.array(timeList)

    # refStateFile = "bspline_states.txt"
    # refStateStream = open(refStateFile, 'w')
    # print >> sys.stdout, "  Saving system states at camera rate generated from B-spline to", refStateFile   
    # saveBsplineRefStates(refStateTimes, poseSplineDv, cself, refStateStream)
    # refStateStream.close()

    # refImuFile = "bspline_imu_meas.txt"
    # saveBsplineRefImuMeas(cself, refImuFile)
  
    modelFile = "bspline_knot_coeff.txt" 
    saveBsplineModel(cself, modelFile)

    cameraIndex = 0
    frameIndex = 0
    imageCornerPoints = cself.CameraChain.getCornersImageSample(poseSplineDv, 0.0, cameraIndex, frameIndex)
    targetCornerPoints = cself.CameraChain.getCornersTargetSample(cameraIndex, frameIndex)
    sampleImageCorners = "bspline_image_corners_{}_{}.txt".format(cameraIndex, frameIndex)
    sampleTargetCorners = "bspline_target_corners_{}_{}.txt".format(cameraIndex, frameIndex)
    np.savetxt(sampleImageCorners,imageCornerPoints, fmt=['%.5f', '%.5f', '%.5f', '%.5f', '%.5f', '%.5f'])
    np.savetxt(sampleTargetCorners,targetCornerPoints, fmt=['%.5f', '%.5f', '%.5f'])


def loadArrayWithHeader(arrayFile):
    with open(arrayFile) as f:
        lines = (line for line in f if not (line.startswith('#') or line.startswith('%')))        
        return np.loadtxt(lines, delimiter=' ', skiprows=0)


def loadBsplineModel(knotCoeffFile):
    splineOrder = 6
    poseSpline = bsplines.BSplinePose(splineOrder, sm.RotationVector() )
    poseSpline.initPoseSplineFromFile(knotCoeffFile)
    print("Initialized a pose spline with {} knots and coefficients {}.".format( \
            poseSpline.knots().size, poseSpline.coefficients().shape))
    poseDv = asp.BSplinePoseDesignVariable( poseSpline )
    return poseDv
