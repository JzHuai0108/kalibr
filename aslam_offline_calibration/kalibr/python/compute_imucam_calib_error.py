# compute errors of camera-IMU calibration results relative to a reference calibration.
import math
import os
import sys
import numpy as np

import kalibr_common as kc
import sm


def replace_all(text, dic):
    """
    replace the substrings in text with those in dic, if the replacement is insensitive to the order of keys and values.
    https://stackoverflow.com/questions/6116978/how-to-replace-multiple-substrings-of-a-string
    :param text:
    :param dic:
    :return:
    """
    for i, j in dic.items():
        text = text.replace(i, j)
    return text


def parseMeanMedianStd(line):
    index = line.find('mean')
    numberline = line[index + len('mean'):]
    dict = {'median': '', 'std': '', ':': ' ', ',': ' ', '#terms': ''}
    clearline = replace_all(numberline, dict)
    numbers = clearline.split()
    return map(float, numbers)

def parseGravityLine(line):
    dict = {'[' : '', '  ': ' ', ']' : ''}
    clearline = replace_all(line, dict)
    numbers = clearline.split()
    return [float(numbers[0]), float(numbers[1]), float(numbers[2])]


def parseImuCameraCalibrationResult(resulttxt):
    statList = []
    nextgravity = False
    with open(resulttxt, 'r') as stream:
        for line in stream:
            if 'Reprojection error ' in line and '[px]' in line:
                stats = parseMeanMedianStd(line)
                statList.extend(stats)
            if 'Gyroscope error ' in line and '[rad/s]' in line:
                stats = parseMeanMedianStd(line)
                statList.extend(stats)
            if 'Accelerometer error ' in line and '[m/s^2]' in line:
                stats = parseMeanMedianStd(line)
                statList.extend(stats)
            if 'Gravity vector in target coords' in line:
                nextgravity = True
                continue
            if nextgravity:
                stats = parseGravityLine(line)
                statList.extend(stats)
                break
    return statList


def findFileInDir(folder, namekeys):
    """
    find a file with keys in filename under dir, It will not go into subfolders.
    :param folder:
    :param namekeys:
    :return:
    """
    for filename in os.listdir(folder):
        status = True
        for key in namekeys:
            if key not in filename:
                status = False
                break
        if status:
            return os.path.join(folder, filename)


def main():
    if len(sys.argv) < 4:
        print("Usage: {} <calibration result folder> <reference camimu yaml> "
              "<reference IMU yaml> <output csv file in append mode>".format(
            sys.argv[0]))
        sys.exit(1)

    folder = sys.argv[1]
    referenceYaml = sys.argv[2]
    referenceImuYaml = sys.argv[3]
    outputCsv = sys.argv[4]

    if not os.path.isdir(folder):
        print("Calibration result {} does not exist!".format(folder))
        sys.exit(2)
    camimuyaml = findFileInDir(folder, ['camchain-imucam', '.yaml'])
    resulttxt = findFileInDir(folder, ['results-', '.txt'])
    if camimuyaml is None or resulttxt is None:
        print("Failed to find camchain-imucam yaml or results txt under {}".format(folder))
        sys.exit(3)

    estimatedChain = kc.CameraChainParameters(camimuyaml)
    referenceChain = kc.CameraChainParameters(referenceYaml)
    numCams = estimatedChain.numCameras()
    calibErrors = []
    for camNr in range(numCams):
        T_cam_imu = estimatedChain.getExtrinsicsImuToCam(camNr)
        ref_T_cam_imu = referenceChain.getExtrinsicsImuToCam(camNr)
        deltaT = ref_T_cam_imu.inverse() * T_cam_imu
        translationError = np.linalg.norm(deltaT.t()) * 1000
        rotationVector = sm.quat2AxisAngle(deltaT.q())
        rotationError = abs(math.atan(math.tan(np.linalg.norm(rotationVector))) * 180 / math.pi)

        deltaTime = referenceChain.getTimeshiftCamImu(camNr) - estimatedChain.getTimeshiftCamImu(camNr)
        deltaLineDelay = estimatedChain.getLineDelay(camNr) / 1000.0
        calibErrors.append([deltaT.t(), rotationVector, deltaTime, deltaLineDelay])

    statList = parseImuCameraCalibrationResult(resulttxt)
    referenceImu = kc.ImuParameters(referenceImuYaml)
    refgravity = referenceImu.getGravityInTarget()
    refunitgravity = refgravity / np.linalg.norm(refgravity)
    gravity = np.array(statList[-3:])
    unitgravity = gravity / np.linalg.norm(gravity)
    unitgravityerror = refunitgravity - unitgravity

    existingCsv = os.path.isfile(outputCsv)
    with open(outputCsv, 'a') as stream:
        # if not existingCsv:
        #     stream.write("folder, translation_error(mm), rotation error(deg), time offset error(us), line delay (us), "
        #                  "reprojection error (mean, median, std), gyro error (mean, median, std), "
        #                  "accel error (mean, median, std)\n")
        # stream.write("{}, {}, {}, {}, {}, {}\n".format(
        #     folder, translationError, rotationError, deltaTime, deltaLineDelay, ', '.join(map(str, statList))))

        if not existingCsv:
            stream.write("folder, translation_error(cam0), rotation error(cam0), translation_error(cam1), rotation error(cam1),"
                         " time offset error(cam0, us), time offset error(cam1, us), unit gravity error, line delay (us), "
                         "reprojection error (mean, median, std), gyro error (mean, median, std), "
                         "accel error (mean, median, std)\n")
        stream.write("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format(
            folder, ','.join(map(str, calibErrors[0][0])), ','.join(map(str, calibErrors[0][1])),
            ','.join(map(str, calibErrors[1][0])), ','.join(map(str, calibErrors[1][1])),
            calibErrors[0][2], calibErrors[1][2],
            ','.join(map(str, unitgravityerror)),
            deltaLineDelay, ','.join(map(str, statList))))

if __name__ == '__main__':
    main()
