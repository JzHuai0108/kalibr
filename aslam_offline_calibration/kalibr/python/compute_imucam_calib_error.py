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


def parseImuCameraCalibrationResult(resulttxt):
    statList = []
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
        print("Usage: {} <calibration result folder> <reference yaml> <output csv file in append mode>".format(
            sys.argv[0]))
        sys.exit(1)

    folder = sys.argv[1]
    referenceYaml = sys.argv[2]
    outputCsv = sys.argv[3]

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
    camNr = 0
    T_cam_imu = estimatedChain.getExtrinsicsImuToCam(camNr)
    ref_T_cam_imu = referenceChain.getExtrinsicsImuToCam(camNr)
    deltaT = ref_T_cam_imu.inverse() * T_cam_imu
    translationError = np.linalg.norm(deltaT.t()) * 1000
    rotationVector = sm.quat2AxisAngle(deltaT.q())
    rotationError = abs(math.atan(math.tan(np.linalg.norm(rotationVector))) * 180 / math.pi)

    deltaTime = referenceChain.getTimeshiftCamImu(camNr) - estimatedChain.getTimeshiftCamImu(camNr)
    deltaTime = deltaTime * 1000000
    deltaLineDelay = estimatedChain.getLineDelay(camNr) / 1000.0

    statList = parseImuCameraCalibrationResult(resulttxt)

    existingCsv = os.path.isfile(outputCsv)
    with open(outputCsv, 'a') as stream:
        if not existingCsv:
            stream.write("folder, translation_error(mm), rotation error(deg), time offset error(us), line delay (us), "
                         "reprojection error (mean, median, std), gyro error (mean, median, std), "
                         "accel error (mean, median, std)\n")
        stream.write("{}, {}, {}, {}, {}, {}\n".format(
            folder, translationError, rotationError, deltaTime, deltaLineDelay, ', '.join(map(str, statList))))


if __name__ == '__main__':
    main()
