
import sys

import numpy as np

import aslam_backend as aopt
import aslam_splines as asp
import bsplines
import sm


def saveBSplineRefPose(times, poseSplineDv, stream=sys.stdout, T_b_c=sm.Transformation()):
    timeOffsetPadding = 0.0  
    for time in times:    
        dv = aopt.Scalar(time)
        timeExpression = dv.toExpression()
        
        if time <= poseSplineDv.spline().t_min() or time >= poseSplineDv.spline().t_max():
            print >> sys.stdout, "Warn: time out of range "
            continue       
        T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
        sm_T_w_c = sm.Transformation(T_w_b.toTransformationMatrix())*T_b_c
        # quatInv used here to convert kalibr's JPL quaternion to Halmilton quaternion
        print >> stream, '%.9f' % time, ' '.join(map(str,sm_T_w_c.t())), ' '.join(map(str, sm.quatInv(sm_T_w_c.q())))


def sampleBSplines(stateTimes, poseSplineDv, gyroBiasSpline, accBiasSpline, timeOffset):
    """
    return:
        1. frameIds.
        2. saved stamps, saved stamp + time offset = state stamp.
        3. states, each state T_WB[xyz, qxyzw], v_W, bg, ba.
    """
    states = np.zeros((len(stateTimes),16))
    measuredTimes = np.zeros(len(stateTimes))
    frameIds = np.zeros(len(stateTimes), dtype=np.int32)

    tmin = max(poseSplineDv.spline().t_min(), gyroBiasSpline.t_min(), accBiasSpline.t_min())
    tmax = min(poseSplineDv.spline().t_max(), gyroBiasSpline.t_max(), accBiasSpline.t_max())

    timeOffsetPadding = 0.0
    frameId = 0
    for time in stateTimes:
        if time <= tmin or time >= tmax:
            print "Warn: time out of range in generating a state"
            continue
        dv = aopt.Scalar(time)
        timeExpression = dv.toExpression()  
        T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
        sm_T_w_b = sm.Transformation(T_w_b.toTransformationMatrix())
        v_w = poseSplineDv.linearVelocity(time).toEuclidean()
        gyro_bias = gyroBiasSpline.eval(time)
        acc_bias = accBiasSpline.eval(time)
        frameIds[frameId] = frameId
        measuredTimes[frameId] = time - timeOffset
        states[frameId, 0:3] = sm_T_w_b.t()
        # quatInv converts JPL quaternion to Halmilton quaternion (x,y,z,w).
        states[frameId, 3:7] = sm.quatInv(sm_T_w_b.q())
        states[frameId, 7:10] = v_w
        states[frameId, 10:13] = acc_bias
        states[frameId, 13:16] = gyro_bias
        frameId += 1
    return frameIds, measuredTimes, states


def toNanosecondString(time):
    return "{}{:09d}".format(int(time), int((time-int(time)) * 1e9))


def saveStates(times, poseSplineDv, gyroBiasSpline, accBiasSpline, timeOffset = 0, stream = sys.stdout):
    '''save the system states at times.'''
    stream.write('vertex index, timestamp [ns], position x [m], position y [m], position z [m], '
                 'quaternion x, quaternion y, quaternion z, quaternion w, velocity x [m/s], '
                 'velocity y [m/s], velocity z [m/s], acc bias x [m/s^2], acc bias y [m/s^2], '
                 'acc bias z [m/s^2], gyro bias x [rad/s], gyro bias y [rad/s], gyro bias z [rad/s]\n')
    frameIds, measuredTimes, states = sampleBSplines( \
            times, poseSplineDv, gyroBiasSpline, accBiasSpline, timeOffset)
    for index, row in enumerate(states):
        msg = ', '.join(map(str, row))
        stream.write("{:d}, {}, {}\n".format(frameIds[index], toNanosecondString(measuredTimes[index]), msg))


def saveBSplineRefImuMeas(cself, filename):
    print >> sys.stdout, "  Saving IMU measurements generated from B-spline to", filename

    idx = 0
    imu = cself.ImuList[idx]    
    poseSplineDv = cself.poseDv
    print("  imuData begin at {:.6f} end at {:.6f} imu time offset {}".format(imu.imuData[0].stamp.toSec(), imu.imuData[-1].stamp.toSec(), imu.timeOffset))
    times = np.array([im.stamp.toSec() + imu.timeOffset for im in imu.imuData \
                      if im.stamp.toSec() + imu.timeOffset > poseSplineDv.spline().t_min() \
                      and im.stamp.toSec() + imu.timeOffset < poseSplineDv.spline().t_max() ])

    predictedAng_body =  np.array([err.getPredictedMeasurement() for err in imu.gyroErrors])    
    predicetedAccel_body =  np.array([err.getPredictedMeasurement() for err in imu.accelErrors])

    gyroBias = imu.gyroBiasDv.spline()    
    gyro_bias_spline = np.array([gyroBias.evalD(t,0) for t in times])
    
    accBias = imu.accelBiasDv.spline()    
    acc_bias_spline = np.array([accBias.evalD(t,0) for t in times])

    print >> sys.stdout, '\tEpitome of predicted inertial measurements'
    print >> sys.stdout, "\t#times", times.shape
    print >> sys.stdout, "\t#gyro", predictedAng_body.shape
    print >> sys.stdout, "\t#accel", predicetedAccel_body.shape
    print >> sys.stdout, "\t#gyro bias", gyro_bias_spline.shape
    print >> sys.stdout, "\t#accel bias", acc_bias_spline.shape
    
    whole=np.concatenate((np.array([times]).T, predictedAng_body, predicetedAccel_body, gyro_bias_spline, acc_bias_spline),axis=1)
    np.savetxt(filename,whole, fmt=['%.9f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f'])


def saveBSpline(cself, outputDir):
    idx = 0
    imu = cself.ImuList[idx]    
    poseSplineDv = cself.poseDv

    computeCheckData = False
    if computeCheckData:
        imuTimes = np.array([im.stamp.toSec() + imu.timeOffset for im in imu.imuData if
                             poseSplineDv.spline().t_min() < im.stamp.toSec() + imu.timeOffset < poseSplineDv.spline().t_max()])
        refPoseStream = open("poses_check.txt", 'w')
        print >> refPoseStream, "%poses generated at the IMU rate from the B-spline: time, T_w_b(txyz, qxyzw)"  
        saveBSplineRefPose(imuTimes, poseSplineDv, stream=refPoseStream)
        refPoseStream.close()

        timeList = list()
        for obs in cself.CameraChain.camList[0].targetObservations:
            frameTime = cself.CameraChain.camList[0].cameraTimeToImuTimeDv.toScalar() + obs.time().toSec() + \
                    cself.CameraChain.camList[0].timeshiftCamToImuPrior
            if frameTime > imuTimes[0] and frameTime < imuTimes[-1]:
                timeList.append(frameTime)
        refStateTimes = np.array(timeList)
        refStateFile = "states_check.txt"
        refStateStream = open(refStateFile, 'w')
        print >> sys.stdout, "  Saving system states at camera rate generated from B-spline to", refStateFile
        imu = cself.ImuList[0]
        gyroBias = imu.gyroBiasDv.spline()
        accBias = imu.accelBiasDv.spline()
        saveStates(refStateTimes, poseSplineDv, gyroBias, accBias, 0.0, refStateStream)
        refStateStream.close()

        refImuFile = "imu_check.txt"
        saveBSplineRefImuMeas(cself, refImuFile)

        landmarks = cself.CameraChain.camList[0].detector.target().points()
        landmarkCsv = "landmarks_check.csv"
        with open(landmarkCsv, 'w') as stream:
            header = ', '.join(["landmark index", "landmark position x [m]",
                                "landmark position y [m]", "landmark position z [m]"])
            stream.write('{}\n'.format(header))
            for index, row in enumerate(landmarks):
                stream.write("{}, {}, {}, {}\n".format(index, row[0], row[1], row[2]))

        # check landmarks observed in an image.
        cameraIndex = 0
        frameIndex = 0
        obs = cself.CameraChain.camList[0].targetObservations[0]
        camTimeOffset = cself.CameraChain.camList[0].cameraTimeToImuTimeDv.toScalar()
        frameTime = camTimeOffset + obs.time().toSec() + \
                    cself.CameraChain.camList[0].timeshiftCamToImuPrior
        print('  Saving landmark observation at {:.6f} time shift prior {} residual time shift {}'.format( \
                frameTime, cself.CameraChain.camList[0].timeshiftCamToImuPrior, camTimeOffset))
        imageCornerPoints = cself.CameraChain.getCornersImageSample(poseSplineDv, 0.0, cameraIndex, frameIndex)
        targetCornerPoints = cself.CameraChain.getCornersTargetSample(cameraIndex, frameIndex)
        sampleImageCorners = "image_corners_{}_{}_check.txt".format(cameraIndex, frameIndex)
        sampleTargetCorners = "landmarks_{}_{}_check.txt".format(cameraIndex, frameIndex)
        np.savetxt(sampleImageCorners,imageCornerPoints, fmt=['%.5f', '%.5f', '%.5f', '%.5f', '%.5f', '%.5f'])
        np.savetxt(sampleTargetCorners,targetCornerPoints, fmt=['%.5f', '%.5f', '%.5f'])

    poseFile = "bspline_pose.txt"
    poseSpline = cself.poseDv.spline()
    poseSpline.saveSplineToFile(poseFile)   
    print("  saved pose B splines of order {} to {}".format(poseSpline.splineOrder(), poseFile))

    imu = cself.ImuList[0]
    gyroBias = imu.gyroBiasDv.spline()   
    accBias = imu.accelBiasDv.spline()

    gyroBiasFile = "bspline_gyro_bias.txt"
    gyroBias.saveSplineToFile(gyroBiasFile)
    print("  saved gyro bias B splines of order {} to {}".format(gyroBias.splineOrder(), gyroBiasFile))

    accBiasFile = "bspline_acc_bias.txt"
    accBias.saveSplineToFile(accBiasFile)    
    print("  saved acc bias B splines of order {} to {}".format(accBias.splineOrder(), accBiasFile))

    print('Saved B splines of start and finish time')
    print('\t\t\t\tstart time\t\tfinish time')
    print('\tposeSpline\t{:.9f}\t{:.9f}'.format(poseSpline.t_min(), poseSpline.t_max()))
    print('\tgyroBias\t{:.9f}\t{:.9f}'.format(gyroBias.t_min(), gyroBias.t_max()))
    print('\taccBias\t\t{:.9f}\t{:.9f}'.format(accBias.t_min(), accBias.t_max()))

def loadArrayWithHeader(arrayFile):
    with open(arrayFile) as f:
        lines = (line for line in f if not (line.startswith('#') or line.startswith('%')))        
        return np.loadtxt(lines, delimiter=' ', skiprows=0)

def getSplineOrder(knotCoeffFile):
    with open(knotCoeffFile, 'r') as stream:
        lineNumber = 0
        for line in stream:
            if lineNumber == 2:
                return int(line.split()[0])
            lineNumber += 1

def loadPoseBSpline(knotCoeffFile):
    splineOrder = getSplineOrder(knotCoeffFile)
    poseSpline = bsplines.BSplinePose(splineOrder, sm.RotationVector())
    poseSpline.initSplineFromFile(knotCoeffFile)
    print("  Initialized a pose spline with {} knots and coefficients {}.".format( \
            poseSpline.knots().size, poseSpline.coefficients().shape))
    poseDv = asp.BSplinePoseDesignVariable(poseSpline)
    return poseDv

def loadBSpline(knotCoeffFile):
    splineOrder = getSplineOrder(knotCoeffFile)
    spline = bsplines.BSpline(splineOrder)
    spline.initSplineFromFile(knotCoeffFile)
    print("  Initialized a Euclidean spline with {} knots and coefficients {}.".format( \
            spline.knots().size, spline.coefficients().shape))
    return asp.EuclideanBSplineDesignVariable(spline)

       

