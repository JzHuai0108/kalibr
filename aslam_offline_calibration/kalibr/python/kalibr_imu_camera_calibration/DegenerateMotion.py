import argparse
import copy
import math
import sys

import numpy as np
import sm


class TranslationX(object):
    def __init__(self, a = 0.85, b = 0.35, f = 0.15, axisIndex = 0):
        self.__a = a
        self.__b = b
        self.__f = f
        self.__axis = axisIndex
        self.__startTime = 20

    def samplePoses(self, referencePose, sampleFrequency=30, sampleDuration=60):
        sampleInterval = 1.0 / sampleFrequency
        timeVector = np.arange(0, sampleDuration, sampleInterval) + self.__startTime
        samplePoses = []
        for t in timeVector:
            v = self.__a * math.sin(2 * math.pi * self.__f * t) + self.__b
            pose = copy.deepcopy(referencePose)
            pose[self.__axis] = v
            samplePoses.append(pose)
        return timeVector, samplePoses


class RotationZ(object):
    def __init__(self, a = 1.5, b = 0.0, f = 0.15, axisIndex = 2):
        self.__a = a
        self.__b = b
        self.__f = f
        self.__axis = axisIndex
        self.__startTime = 20

    def samplePoses(self, referencePose, sampleFrequency=30, sampleDuration=60):
        sampleInterval = 1.0 / sampleFrequency
        timeVector = np.arange(0, sampleDuration, sampleInterval) + self.__startTime
        samplePoses = []
        for t in timeVector:
            v = self.__a * math.sin(2 * math.pi * self.__f * t) + self.__b
            pose = copy.deepcopy(referencePose)
            aa = np.zeros((3, 1))
            aa[self.__axis] = v
            dq = sm.axisAngle2quat(aa)
            # quatInv converts JPL quaternion to Halmilton quaternion (x,y,z,w).
            pose[3:7] = sm.quatInv(np.dot(sm.quatPlus(sm.quatInv(np.array(pose[3:7]))), dq))
            samplePoses.append(pose)
        return timeVector, samplePoses


class TranslationY(TranslationX):
    def __init__(self, a = 0.7, b = 0.25, f = 0.15, axisIndex = 1):
        super(TranslationY, self).__init__(a, b, f, axisIndex)


class TranslationZ(TranslationX):
    def __init__(self, a = 0.85, b = 1.3, f = 0.15, axisIndex = 2):
        super(TranslationZ, self).__init__(a, b, f, axisIndex)


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('axis', type=int, default=0,
                        help='simulate degenerate motion at which axis? 0, 1, 2: translate along x, y, z; '
                             '3, 4, 5: rotate about x, y, z (default: %(default)s)')
    parser.add_argument('--output-file', type=str, default=0, dest='outputFile',
                        help='output pose file, each line time (sec), txyz, qxyzw', required=False)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(2)
    return parser.parse_args()


def main():
    parsed = parseArgs()

    if parsed.axis == 0:
        sampler = TranslationX()
        referencePose = [0.35, 0.35, 1.5, 0.9984003, -0.0399641, -0.0016026, -0.0399641]
    elif parsed.axis == 1:
        sampler = TranslationY()
        referencePose = [0.35, 0.35, 1.5, 0.9984003, -0.0399641, -0.0016026, -0.0399641]
    elif parsed.axis == 2:
        sampler = TranslationZ()
        referencePose = [0.35, 0.35, 1.5, 0.9984003, -0.0399641, -0.0016026, -0.0399641]
    else:
        sampler = RotationZ(axisIndex = parsed.axis - 3)
        referencePose = [0.35, 0.35, 1.3, 0.9984003, -0.0399641, -0.0016026, -0.0399641]

    times, poses = sampler.samplePoses(referencePose)
    with open(parsed.outputFile, 'w') as stream:
        for index, time in enumerate(times):
            stream.write('{:.9f}, {}\n'.format(time, ', '.join(map(str, poses[index]))))


if __name__ == "__main__":
    main()
