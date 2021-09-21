import numpy as np
from . import BSplineIO

def testSerializePoses():
    timeList, smTList = BSplineIO.generateRandomPoses()
    testFile = "testPoses.txt"
    BSplineIO.savePoses(timeList, smTList, testFile)

    newTimeList, newsmTList = BSplineIO.loadPoses(testFile)
    assert len(timeList) == len(newTimeList)
    for index, time in enumerate(timeList):
        assert time == newTimeList[index]
        assert np.allclose(smTList[index].T(), newsmTList[index].T())

    nulltxPoseList = BSplineIO.projectPoses(smTList, 'tz')
    nullrzPoseList = BSplineIO.projectPoses(smTList, 'rz')
    for pose in nulltxPoseList:
        print("{}\n".format(pose.T()))
    for pose in nullrzPoseList:
        print("{}\n".format(pose.T()))

