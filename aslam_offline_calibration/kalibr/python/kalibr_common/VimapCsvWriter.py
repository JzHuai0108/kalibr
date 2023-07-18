"""
Input and output functions for CSV files of a maplab vimap
"""
import os
import numpy as np

from kalibr_common import ConfigReader as cr
import sm


def saveObservations(landmarkObservations, observationCsv):
    with open(observationCsv, 'w') as stream:
        header = ', '.join(["vertex index", "frame index", "keypoint index", "landmark index"])
        stream.write('{}\n'.format(header))
        for landmarkId, observationList in sorted(landmarkObservations.items()):
            for observation in observationList:
                stream.write('{}, {}, {}, {}\n'.format(observation[0], observation[1], observation[2], landmarkId))


def saveTracks(frameKeypointMap, trackCsv):
    with open(trackCsv, 'w') as stream:
        header = ', '.join(
            ["timestamp [ns]", "vertex index", "frame index", "keypoint index", "keypoint measurement 0 [px]",
                "keypoint measurement 1 [px]", "keypoint measurement uncertainty", "keypoint scale",
                "keypoint track id"])
        stream.write('{}\n'.format(header))
        for vertexId, frameKeypointList in frameKeypointMap.items():
            for frameKeypoints in frameKeypointList:
                timeString = acvTimeToNanosecondString(frameKeypoints[0])
                cameraId = frameKeypoints[1]
                for keypoint in frameKeypoints[2]:
                    stream.write('{}, {}, {}, {:d}, {:.5f}, {:.5f}, {}, {}, {}\n'.format(
                        timeString, vertexId, cameraId, keypoint[1], keypoint[2], keypoint[3], keypoint[4],
                        keypoint[5], keypoint[6]))

def findNearestTimeSince(timeList, index, time):
    """
    timeList: acv Time list in increasing order.
    index: start index.
    time: acv Time.
    """
    diff = abs(timeList[index].toSec() - time.toSec())

    while index + 1 < len(timeList):
        nextdiff = abs(timeList[index + 1].toSec() - time.toSec())
        if nextdiff < diff:
            diff = nextdiff
            index = index + 1
        else:
            break
    return index


def acvTimeToNanosecondString(acvTime):
    return "{}{:09d}".format(acvTime.sec, acvTime.nsec)


def writeOpenCVYaml(initialParameters, yamlFile):
    """
    write the params in initialParameters to a yaml in opencv format
    """
    with open(yamlFile, 'w') as stream:
        stream.write("%YAML:1.0\n")
        for key, val in initialParameters.items():
            if key.startswith("cam"):
                stream.write("{}:\n".format(key))
                for key_c, val_c in val.items():
                    if key_c == "T_imu_cam":
                        stream.write("  T_imu_cam: !!opencv-matrix\n")
                        stream.write("    cols: {}\n".format(val_c.shape[1]))
                        stream.write("    rows: {}\n".format(val_c.shape[0]))
                        stream.write("    dt: d\n")
                        stream.write("    data: [")
                        for i in range(val_c.shape[0]):
                            for j in range(val_c.shape[1]):
                                if i == val_c.shape[0] - 1 and j == val_c.shape[1] - 1:
                                    stream.write("{}".format(val_c[i][j]))
                                else:
                                    stream.write("{}, ".format(val_c[i][j]))
                            if i < val_c.shape[0] - 1:
                                stream.write("\n           ")
                            else:
                                stream.write("]\n")
                    else:
                        stream.write("  {}: {}\n".format(key_c, val_c))
            else:
                stream.write("{}: [".format(key))
                for i in range(len(val)):
                    if i == len(val) - 1:
                        stream.write("{}]\n".format(val[i]))
                    else:
                        stream.write("{}, ".format(val[i]))


def saveVimap(cself, outputDir):
    """save extracted image keypoints and IMU data in maplab csv format
    If bag_from_to is provided to the entry program, this will only save the data within that interval.
    If perform_synchronization is provided to the entry program, this will save synced local times for cameras and IMUs.
    """
    vertexCsv = os.path.join(outputDir, "vertices.csv")
    # create a vertex for every image of camera 0.
    numFrames = len(cself.CameraChain.camList[0].targetObservations)
    states = np.zeros((numFrames, 16))
    timeList = []
    vertexIdList = []
    T_c_b = cself.CameraChain.camList[0].T_extrinsic.T()
    for vertexId, obs in enumerate(cself.CameraChain.camList[0].targetObservations):
        T_t_b = np.dot(obs.T_t_c().T(), T_c_b)
        sm_T_w_b = sm.Transformation(T_t_b)
        states[vertexId, 0:3] = sm_T_w_b.t()
        # quatInv converts JPL quaternion to Halmilton quaternion (x,y,z,w).
        states[vertexId, 3:7] = sm.quatInv(sm_T_w_b.q())
        # velocity, gyro bias, and accelerometer bias are initialized to zeros.
        timeList.append(obs.time())
        vertexIdList.append(vertexId)

    with open(vertexCsv, 'w') as stream:
        stream.write('vertex index, timestamp [ns], position x [m], position y [m], position z [m], '
                     'quaternion x, quaternion y, quaternion z, quaternion w, velocity x [m/s], '
                     'velocity y [m/s], velocity z [m/s], acc bias x [m/s^2], acc bias y [m/s^2], '
                     'acc bias z [m/s^2], gyro bias x [rad/s], gyro bias y [rad/s], gyro bias z [rad/s]\n')
        for index, row in enumerate(states):
            msg = ', '.join(map(str, row))
            stream.write("{:d}, {}, {}\n".format(vertexIdList[index], acvTimeToNanosecondString(timeList[index]), msg))

    yamlFile = os.path.join(outputDir, "initial_camchain_imu.yaml")
    initialParameters = cr.ParametersBase(yamlFile, "CameraChainParameters", True)
    for camNr, camera in enumerate(cself.CameraChain.camList):
        T_c_b = cself.CameraChain.getResultTrafoImuToCam(camNr)
        T_b_c = T_c_b.inverse()
        timeShiftPrior = float(camera.timeshiftCamToImuPrior)
        camName = "cam{}".format(camNr)
        initialParameters.data[camName] = dict()
        initialParameters.data[camName]["T_imu_cam"] = T_b_c.T()
        initialParameters.data[camName]["timeshift_cam_imu"] = timeShiftPrior
    initialParameters.data["gravity_in_target"] = cself.CameraChain.getEstimatedGravity().tolist()
    # initialParameters.writeYaml()
    writeOpenCVYaml(initialParameters.data, yamlFile)

    landmarks = cself.CameraChain.camList[0].detector.target().points()

    frameKeypointList = list()
    landmarkObservations = dict()

    for landmarkId in range(len(landmarks)):
        landmarkObservations[landmarkId] = list()

    # assign observations into vertices and cameras.
    for cameraId in range(len(cself.CameraChain.camList)):
        reprojectionSigma = cself.CameraChain.camList[cameraId].cornerUncertainty
        vertexId = 0
        for obs in cself.CameraChain.camList[cameraId].targetObservations:
            # find the corresponding vertex by timestamp
            rawTime = obs.time()
            vertexId = findNearestTimeSince(timeList, vertexId, rawTime)
            assert timeList[vertexId].toSec() - rawTime.toSec() < 0.05
            cornersInImage = obs.getCornersImageFrame()
            landmarkIds = obs.getCornersIdx()

            frameKeypoints = list()
            for keypointIndex in range(cornersInImage.shape[0]):
                landmarkObservations[landmarkIds[keypointIndex]].append((vertexId, cameraId, keypointIndex))
                frameKeypoints.append((landmarkIds[keypointIndex], keypointIndex, cornersInImage[keypointIndex][0],
                cornersInImage[keypointIndex][1], reprojectionSigma, 12, -1))
            frameKeypointList.append((rawTime, vertexId, cameraId, frameKeypoints))

    observationCsv = os.path.join(outputDir, "observations.csv")
    saveObservations(landmarkObservations, observationCsv)

    trackCsv = os.path.join(outputDir, "tracks.csv")
    saveTracks(frameKeypointList, trackCsv)

    landmarkCsv = os.path.join(outputDir, "landmarks.csv")
    with open(landmarkCsv, 'w') as stream:
        header = ', '.join(["landmark index", "landmark position x [m]",
                            "landmark position y [m]", "landmark position z [m]"])
        stream.write('{}\n'.format(header))
        for index, row in enumerate(landmarks):
            stream.write("{}, {}, {}, {}\n".format(index, row[0], row[1], row[2]))
    for index, imu in enumerate(cself.ImuList):
        imuCsv = os.path.join(outputDir, "imu{}.csv".format(index))
        with open(imuCsv, "w") as stream:
            header = ', '.join(["timestamp [ns]", "acc x [m/s^2]", "acc y [m/s^2]", "acc z [m/s^2]",
                                "gyro x [rad/s]", "gyro y [rad/s]", "gyro z [rad/s]"])
            stream.write('{}\n'.format(header))
            for time, omega, alpha in imu.dataset:
                alphaString = ', '.join(map(str, alpha))
                omegaString = ', '.join(map(str, omega))
                stream.write("{}, {}, {}\n".format(acvTimeToNanosecondString(time), alphaString, omegaString))
