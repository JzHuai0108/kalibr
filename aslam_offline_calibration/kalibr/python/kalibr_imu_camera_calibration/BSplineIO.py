from __future__ import print_function
import math
import sys

import numpy as np

import aslam_backend as aopt
import aslam_splines as asp
import bsplines
import sm


def sampleAndSaveBSplinePoses(times, poseSplineDv, stream=sys.stdout, T_b_c=sm.Transformation()):
    timeOffsetPadding = 0.0  
    for time in times:
        dv = aopt.Scalar(time)
        timeExpression = dv.toExpression()
        
        if time <= poseSplineDv.spline().t_min() or time >= poseSplineDv.spline().t_max():
            print("Warn: time out of range ", file=sys.stdout)
            continue
        T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
        sm_T_w_c = sm.Transformation(T_w_b.toTransformationMatrix())*T_b_c
        # quatInv used here to convert kalibr's JPL quaternion to Hamilton quaternion
        print('{:.9f}, {}, {}'.format(time, ','.join(map(str,sm_T_w_c.t())),
                                      ','.join(map(str, sm.quatInv(sm_T_w_c.q())))), file=stream)

def sampleBSplinePoses(stateTimes, poseSplineDv):
    """
    return:
        1. frameIds.
        2. saved state timestamps.
        3. states, each state T_WB[xyz, qxyzw], v_W.
    """
    states = np.zeros((len(stateTimes),16))
    measuredTimes = np.zeros(len(stateTimes))
    frameIds = np.zeros(len(stateTimes), dtype=np.int32)

    tmin = poseSplineDv.spline().t_min()
    tmax = poseSplineDv.spline().t_max()

    timeOffsetPadding = 0.0
    frameId = 0
    for time in stateTimes:
        if time <= tmin or time >= tmax:
            print("Warn: time out of range in generating a state")
            continue
        dv = aopt.Scalar(time)
        timeExpression = dv.toExpression()
        T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
        sm_T_w_b = sm.Transformation(T_w_b.toTransformationMatrix())
        v_w = poseSplineDv.linearVelocity(time).toEuclidean()

        frameIds[frameId] = frameId
        measuredTimes[frameId] = time
        states[frameId, 0:3] = sm_T_w_b.t()
        # quatInv converts JPL quaternion to Halmilton quaternion (x,y,z,w).
        states[frameId, 3:7] = sm.quatInv(sm_T_w_b.q())
        states[frameId, 7:10] = v_w
        frameId += 1
    return frameIds, measuredTimes, states


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
    if gyroBiasSpline:
        tmin = max(poseSplineDv.spline().t_min(), gyroBiasSpline.t_min(), accBiasSpline.t_min())
        tmax = min(poseSplineDv.spline().t_max(), gyroBiasSpline.t_max(), accBiasSpline.t_max())
    else:
        tmin = poseSplineDv.spline().t_min()
        tmax = poseSplineDv.spline().t_max()

    timeOffsetPadding = 0.0
    frameId = 0
    for time in stateTimes:
        if time <= tmin or time >= tmax:
            print("Warn: time out of range in generating a state")
            continue
        dv = aopt.Scalar(time)
        timeExpression = dv.toExpression()  
        T_w_b = poseSplineDv.transformationAtTime(timeExpression, timeOffsetPadding, timeOffsetPadding)
        sm_T_w_b = sm.Transformation(T_w_b.toTransformationMatrix())
        v_w = poseSplineDv.linearVelocity(time).toEuclidean()
        if gyroBiasSpline:
            gyro_bias = gyroBiasSpline.eval(time)
            acc_bias = accBiasSpline.eval(time)
        else:
            gyro_bias = np.zeros(3)
            acc_bias = np.zeros(3)
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


def secondToNanosecondString(time):
    return "{}{:09d}".format(int(time), int((time-int(time)) * 1e9))


def saveStates(times, poseSplineDv, gyroBiasSpline, accBiasSpline, timeOffset = 0, stream = sys.stdout):
    '''save the system states at times.'''
    stream.write('vertex index, timestamp [ns], position x [m], position y [m], position z [m], '
                 'quaternion x, quaternion y, quaternion z, quaternion w, velocity x [m/s], '
                 'velocity y [m/s], velocity z [m/s], acc bias x [m/s^2], acc bias y [m/s^2], '
                 'acc bias z [m/s^2], gyro bias x [rad/s], gyro bias y [rad/s], gyro bias z [rad/s]\n')
    frameIds, measuredTimes, states = sampleBSplines(
            times, poseSplineDv, gyroBiasSpline, accBiasSpline, timeOffset)
    for index, row in enumerate(states):
        msg = ', '.join(map(str, row))
        stream.write("{:d}, {}, {}\n".format(frameIds[index], secondToNanosecondString(measuredTimes[index]), msg))

def saveCameraStates(times, poseSplineDv, stream = sys.stdout):
    '''save the system states at times.'''
    stream.write('vertex index, timestamp [ns], position x [m], position y [m], position z [m], '
                 'quaternion x, quaternion y, quaternion z, quaternion w\n')
    frameIds, measuredTimes, states = sampleBSplinePoses(times, poseSplineDv)
    for index, row in enumerate(states):
        msg = ', '.join(map(str, row))
        stream.write("{:d}, {}, {}\n".format(frameIds[index], secondToNanosecondString(measuredTimes[index]), msg))

def saveImuMeasurementsFromPoseBSpline(cself, filename):
    """
    generate IMU data from pose bspline, gyro bspline, and accelerometer bspline.
    :param cself:
    :param filename:
    :return: time, predicted angular rate, predicted accelerometer data, gyro bias, accelerometer bias
    """
    print("  Saving IMU measurements generated from B-spline to", filename, file=sys.stdout)

    idx = 0
    imu = cself.ImuList[idx]
    poseSplineDv = cself.poseDv
    print("  imuData begin at {:.6f} end at {:.6f} imu time offset {}".format(imu.imuData[0].stamp.toSec(), imu.imuData[-1].stamp.toSec(), imu.timeOffset))
    times = np.array([im.stamp.toSec() + imu.timeOffset for im in imu.imuData \
                      if poseSplineDv.spline().t_min() < im.stamp.toSec() + imu.timeOffset < poseSplineDv.spline().t_max()])

    predictedAng_body =  np.array([err.getPredictedMeasurement() for err in imu.gyroErrors])    
    predictedAccel_body =  np.array([err.getPredictedMeasurement() for err in imu.accelErrors])

    gyro_bias_spline = np.array([imu.evaluateGyroBias(t) for t in times])
    acc_bias_spline = np.array([imu.evaluateAccelerometerBias(t) for t in times])

    print('\tEpitome of predicted inertial measurements', file=sys.stdout)
    print("\t#times", times.shape, file=sys.stdout)
    print("\t#gyro", predictedAng_body.shape, file=sys.stdout)
    print("\t#accel", predictedAccel_body.shape, file=sys.stdout)
    print("\t#gyro bias", gyro_bias_spline.shape, file=sys.stdout)
    print("\t#accel bias", acc_bias_spline.shape, file=sys.stdout)

    predictedImu=np.concatenate((np.array([times]).T, predictedAng_body, predictedAccel_body, gyro_bias_spline, acc_bias_spline),axis=1)
    np.savetxt(filename,predictedImu, fmt=['%.9f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f', '%.7f'])
    return predictedImu


def saveBSpline(cself, outputDir):
    idx = 0
    imu = cself.ImuList[idx]    
    poseSplineDv = cself.poseDv

    imuTimes = np.array([im.stamp.toSec() + imu.timeOffset for im in imu.imuData if
                         poseSplineDv.spline().t_min() < im.stamp.toSec() + imu.timeOffset < poseSplineDv.spline().t_max()])
    refPoseStream = open("sampled_poses.txt", 'w')
    print("%poses generated at the IMU rate from the B-spline: time, T_w_b(txyz, qxyzw)", file=refPoseStream)
    sampleAndSaveBSplinePoses(imuTimes, poseSplineDv, stream=refPoseStream)
    refPoseStream.close()

    computeCheckData = False
    if computeCheckData:
        refImuFile = "imu_check.txt"
        saveImuMeasurementsFromPoseBSpline(cself, refImuFile)

        # check landmarks observed in an image.
        cameraIndex = 0
        frameIndex = 0
        obs = cself.CameraChain.camList[0].targetObservations[0]
        camTimeOffset = cself.CameraChain.camList[0].cameraTimeToImuTimeDv.toScalar()
        frameTime = camTimeOffset + obs.time().toSec() + \
                    cself.CameraChain.camList[0].timeshiftCamToImuPrior
        print('  Saving landmark observation at {:.6f} time shift prior {} residual time shift {}'.format( \
                frameTime, cself.CameraChain.camList[0].timeshiftCamToImuPrior, camTimeOffset))
        imageCornerPoints = cself.CameraChain.getReprojectedCorners(poseSplineDv, 0.0, cameraIndex, frameIndex)
        targetCornerPoints = cself.CameraChain.getCornersTargetSample(cameraIndex, frameIndex)
        sampleImageCorners = "image_corners_{}_{}_check.txt".format(cameraIndex, frameIndex)
        sampleTargetCorners = "landmarks_{}_{}_check.txt".format(cameraIndex, frameIndex)
        np.savetxt(sampleImageCorners,imageCornerPoints, fmt=['%.5f', '%.5f', '%.5f', '%.5f', '%.5f', '%.5f'])
        np.savetxt(sampleTargetCorners,targetCornerPoints, fmt=['%.5f', '%.5f', '%.5f'])

    poseFile = "bspline_pose.txt"
    poseSpline = cself.poseDv.spline()
    poseSpline.saveSplineToFile(poseFile)   
    print("  saved pose B splines of order {} to {}".format(poseSpline.splineOrder(), poseFile))

    print('Saved B splines of start and finish time')
    print('\t\t\t\tstart time\t\tfinish time')
    print('\tposeSpline\t{:.9f}\t{:.9f}'.format(poseSpline.t_min(), poseSpline.t_max()))
    imu = cself.ImuList[0]
    if not imu.constantBias:
        gyroBias = imu.gyroBiasDv.spline()
        accBias = imu.accelBiasDv.spline()
        print('\tgyroBias\t{:.9f}\t{:.9f}'.format(gyroBias.t_min(), gyroBias.t_max()))
        print('\taccBias\t\t{:.9f}\t{:.9f}'.format(accBias.t_min(), accBias.t_max()))

        gyroBiasFile = "bspline_gyro_bias.txt"
        gyroBias.saveSplineToFile(gyroBiasFile)
        print("  saved gyro bias B splines of order {} to {}".format(gyroBias.splineOrder(), gyroBiasFile))

        accBiasFile = "bspline_acc_bias.txt"
        accBias.saveSplineToFile(accBiasFile)
        print("  saved acc bias B splines of order {} to {}".format(accBias.splineOrder(), accBiasFile))


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


def generateRandomPoses():
    timeList = []
    smTList = []
    for i in range(4):
        q = np.random.rand(4)
        q = q / np.linalg.norm(q)
        p = np.random.rand(3)
        smTList.append(sm.Transformation(q, p))
        timeList.append(i + 1)
    return timeList, smTList


def savePoses(timeList, smTList, outputfile):
    with open(outputfile, "w") as stream:
        print("%time (sec), T_w_c (txyz, qxyzw)", file=stream)
        for index, T in enumerate(smTList):
            print('{:.9f}, {}, {}'.format(timeList[index], ','.join(map(str, T.t())),
                                          ','.join(map(str, sm.quatInv(T.q())))), file=stream)


def loadPoses(poseFile):
    """
        load poses from a file, each line
        time(sec), T_w_b(txyz, qxyzw)
        return time list, sm Transformation list
    """
    timeList = []
    smTransformationList = []
    with open(poseFile, 'r') as stream:
        for line in stream:
            if line.startswith('%') or line.startswith('#'):
                continue
            segments = line.split(',')
            time = float(segments[0])
            pxyz = np.array(map(float, segments[1:4]))
            qxyzw = np.array(map(float, segments[4:8]))
            qxyzw[:3] = - qxyzw[:3]  # Hamilton to JPL convention.
            T_w_b = sm.Transformation(qxyzw, pxyz)
            timeList.append(time)
            smTransformationList.append(T_w_b)
    return timeList, smTransformationList


def projectPoses(smTransformList, projectionCode, maxTranslation=2.5):
    """project the pose along specific axis, for instance project along x, and make other components take the average."""
    componentList = [sm.fromTEuler(transform.T()) for transform in smTransformList] # tx, ty, tz, theta x, theta y, theta z.
    componentArray = np.array(componentList)
    meanValues = []
    for i in range(6):
        column = componentArray[:, i]
        if i < 3:
            mean = column[(-maxTranslation < column) & (column < maxTranslation)].mean()
        else:
            mean = column.mean()
        meanValues.append(mean)
    meanValues = np.array(meanValues)

    projectDict = {'tx' : 0,
                   'ty' : 1,
                   'tz' : 2,
                   'rx' : 3,
                   'ry' : 4,
                   'rz' : 5}
    variableIndex = projectDict[projectionCode]

    # use average values for not affected components.
    rows = componentArray.shape[0]
    for i in range(6):
        if i == variableIndex:
            for j in range(rows):
                if componentArray[j, i] > maxTranslation or componentArray[j, i] < -maxTranslation:
                    if j == 0:
                        componentArray[j, i]= math.copysign(1, componentArray[j, i]) * maxTranslation
                    else:
                        componentArray[j, i] = componentArray[j-1, i]
        else:
            componentArray[:, i] = meanValues[i]

    return [sm.Transformation(sm.toTEuler(row)) for row in componentArray]


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


def __isBSplineFile(filename):
    with open(filename, 'r') as stream:
        first_line = stream.readline()
        if 'splineOrder' in first_line:
            return True
    return False


def ensureContinuousRotationVectors(curve):
    """
    Ensures that the rotation vector does not flip and enables a continuous trajectory modeling.
    Updates curves in place.
    Copied from RsCalibrator.py
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


def generateInitialSpline(times, smTransformList, splineOrder = 5, timeOffsetPadding = 0.05, numberOfKnots = None):
    """
    Adapted from RsCalibrator.py
    :param times: list of time in seconds
    :param smTransformList: list of sm transforms for each timestamp
    :param splineOrder:
    :param timeOffsetPadding:
    :param numberOfKnots:
    :param framerate:
    :return:
    """
    poseSpline = bsplines.BSplinePose(splineOrder, sm.RotationVector())
    curve = np.matrix([ poseSpline.transformationToCurveValue(transform.T()) for transform in smTransformList]).T
    if np.isnan(curve).any():
        raise RuntimeError("Nans in curve values")
        sys.exit(0)
    # Add padding on either end to allow the spline to slide during optimization.
    times = np.hstack((times[0] - (timeOffsetPadding * 2.0), times, times[-1] + (timeOffsetPadding * 2.0)))
    curve = np.hstack((curve[:,0], curve, curve[:,-1]))
    ensureContinuousRotationVectors(curve)

    seconds = times[-1] - times[0]
    framerate = len(times) / seconds
    if numberOfKnots is not None:
        knots = numberOfKnots
    else:
        knots = int(round(seconds * framerate/3))

    print("Initializing a pose spline with %d knots (%f knots per second over %f seconds)" % ( knots, 100, seconds))
    poseSpline.initPoseSplineSparse(times, curve, knots, 1e-4)
    return poseSpline


def selectiveLoadPoseBSpline(filename):
    isBSpline = __isBSplineFile(filename)
    if isBSpline:
        return loadPoseBSpline(filename)
    else:
        timeList, smTransformList = loadPoses(filename)
        poseSpline = generateInitialSpline(timeList, smTransformList)
        poseDv = asp.BSplinePoseDesignVariable(poseSpline)
        return poseDv
