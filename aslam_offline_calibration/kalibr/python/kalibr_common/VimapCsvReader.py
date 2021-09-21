import csv
import os

import numpy as np

import aslam_cv as acv
import sm


class FrameObservation(object):
    def __init__(self):
        self._stamp = acv.Time(0.0)
        self._T_t_c = sm.Transformation()
        self.landmarkList = [] # observed landmarks
        self.landmarkIdList = []
        self.observationList = [] # landmark observations in image

    def T_t_c(self):
        """
        :return: transform from camera to target
        """
        return self._T_t_c

    def set_T_t_c(self, T_t_c):
        self._T_t_c = T_t_c

    def time(self):
        return self._stamp

    def setTime(self, time):
        self._stamp = time

    def getCornersImageFrame(self):
        return self.observationList

    def getCornersTargetFrame(self):
        return self.landmarkList

    def getCornersIdx(self):
        return self.landmarkIdList

    def hasSuccessfulObservation(self):
        return len(self.observationList) > 0

    def getCornerReprojection(self, cameraGeometry):
        """

        :param cameraGeometry: eg, DistortedPinholeCameraGeometry
        :return:
        """
        raise NotImplementedError

    @staticmethod
    def landmarkSentinel():
        return np.array([0, 0, 1e8])

    def removeUnusedKeypoints(self):
        keep = []

        for landmark in self.landmarkList:
            if np.allclose(landmark, FrameObservation.landmarkSentinel()):
                keep.append(False)
            else:
                keep.append(True)

        self.landmarkList = [landmark for index, landmark in enumerate(self.landmarkList) if keep[index]]
        self.landmarkIdList = [id for index, id in enumerate(self.landmarkIdList) if keep[index]]
        self.observationList = [observation for index, observation in enumerate(self.observationList) if keep[index]]

    def appendLandmark(self, landmark, landmarkId):
        self.landmarkList.append(landmark)
        self.landmarkIdList.append(landmarkId)

    def appendObservation(self, observation):
        self.observationList.append(observation)

    def __str__(self):
        header = "{:.9f} {} {}\n".format(self._stamp.toSec(), ' '.join(map(str, self._T_t_c.t())),
                                        ' '.join(map(str, sm.quatInv(self._T_t_c.q()))))
        msg = header
        for index, landmark in enumerate(self.landmarkList):
            msg += '{}: {} {}\n'.format(self.landmarkIdList[index], landmark, self.observationList[index])
        return msg


class VimapCsvReaderIterator(object):
  def __init__(self, dataset, indices=None):
    self.dataset = dataset
    if indices is None:
      self.indices = np.arange(dataset.numImages())
    else:
      self.indices = indices
    self.iter = self.indices.__iter__()

  def __iter__(self):
    return self

  def __next__(self):
    idx = next(self.iter)
    return self.dataset.getImage(idx)

  next = __next__  # Python 2

class VimapCsvReader(object):
    def __init__(self, folder, topic, from_to=None, perform_synchronization=False):
        self.folder = folder
        self.bagfile = folder
        self.camera_index = int(topic[topic.find('image_raw') - 2])
        self.topic = topic
        self.from_to = from_to
        self.numCameras = -1
        self.numVertices = -1
        self.tracks = None
        self.observations = None
        self.landmarks = None
        self.vertices = None
        self.loadVimap(folder)
        self.targetObservations = self.composeFeatureAssociations()
        self.indices = np.arange(len(self.targetObservations))
        if from_to:
            self.indices = self.truncateIndicesFromTime(self.indices, from_to)

    def hasFeatureAssociations(self):
        return True

    def composeFeatureAssociations(self):
        targetObservations = []
        for j in range(self.numVertices):
            targetObservations.append(FrameObservation())

        # camera poses for target observations will be initialized in calibrator by PnP.

        for keypoint in self.tracks:
            if keypoint.camera_idx != self.camera_index:
                continue
            assert keypoint.keypoint_idx == len(
                targetObservations[keypoint.vertex_idx].landmarkList)
            assert keypoint.keypoint_idx == len(
                targetObservations[keypoint.vertex_idx].observationList)
            targetObservations[keypoint.vertex_idx].setTime(keypoint.timestamp)
            targetObservations[keypoint.vertex_idx].appendLandmark(
                FrameObservation.landmarkSentinel(), -1)
            targetObservations[keypoint.vertex_idx].appendObservation(keypoint.measurement)

        for observation in self.observations:
            if observation.camera_idx != self.camera_index:
                continue
            targetObservations[observation.vertex_idx].landmarkList[
                observation.keypoint_idx] = self.landmarks[observation.landmark_idx].position
            targetObservations[observation.vertex_idx].landmarkIdList[
                observation.keypoint_idx] = observation.landmark_idx

        for j in range(self.numVertices):
            targetObservations[j].removeUnusedKeypoints()

        return targetObservations

    def truncateIndicesFromTime(self, indices, bag_from_to):
        # get the timestamps
        timestamps = [observation.time().toSec() for observation in self.targetObservations]

        bagstart = min(timestamps)
        baglength = max(timestamps) - bagstart

        # some value checking
        if bag_from_to[0] >= bag_from_to[1]:
            raise RuntimeError("Bag start time must be bigger than end time.".format(bag_from_to[0]))
        if bag_from_to[0] < 0.0:
            sm.logWarn("Bag start time of {0} s is smaller 0".format(bag_from_to[0]))
        if bag_from_to[1] > baglength:
            sm.logWarn("Bag end time of {0} s is bigger than the total length of {1} s".format(
                bag_from_to[1], baglength))

        # find the valid timestamps
        valid_indices = []
        for idx, timestamp in enumerate(timestamps):
            if timestamp >= (bagstart + bag_from_to[0]) and timestamp <= (bagstart + bag_from_to[1]):
                valid_indices.append(idx)
        sm.logWarn(
            "VimapCsvReader: truncated {0} / {1} images.".format(len(indices) - len(valid_indices),
                                                                 len(indices)))
        return valid_indices

    def __iter__(self):
        # Reset the bag reading
        return self.readDataset()

    def readDataset(self):
        return VimapCsvReaderIterator(self, self.indices)

    def numImages(self):
        return len(self.indices)

    def getImage(self, idx):
        return self.targetObservations[idx]

    def getFeatureAssociations(self):
        pnpObservations = []
        for index in self.indices:
            observation = acv.PnPObservation()
            frameObservation = self.targetObservations[index]
            corners = frameObservation.getCornersImageFrame()
            if len(corners) > 7:
                observation.setCornersImageFrame(np.array(corners))
                observation.setCornersTargetFrame(np.array(frameObservation.getCornersTargetFrame()))
                observation.setCornersIdx(np.array(frameObservation.getCornersIdx()))
                observation.setTime(frameObservation.time())
                observation.set_T_t_c(frameObservation.T_t_c())
                pnpObservations.append(observation)
        return pnpObservations

    def loadVimap(self, folder):
        trackCsv = os.path.join(folder, 'tracks.csv')
        self.tracks, self.numCameras, self.numVertices = loadTrackCsv(trackCsv)
        observationCsv = os.path.join(folder, "observations.csv")
        self.observations = loadObservationCsv(observationCsv)
        landmarkCsv = os.path.join(folder, "landmarks.csv")
        self.landmarks = loadLandmarkCsv(landmarkCsv)
        vertexCsv = os.path.join(folder, "vertices.csv")
        self.vertices = loadVertexCsv(vertexCsv)
        assert (self.numVertices == len(self.vertices))


class VimapImuCsvReaderIterator(object):
    def __init__(self, dataset, indices=None):
        self.dataset = dataset
        if indices is None:
            self.indices = np.arange(dataset.numMessages())
        else:
            self.indices = indices
        self.iter = self.indices.__iter__()

    def __iter__(self):
        return self

    def __next__(self):
        idx = next(self.iter)
        return self.dataset.getMessage(idx)

    next = __next__  # Python 2

class VimapImuCsvReader(object):
    def __init__(self, folder, topic, from_to=None, perform_synchronization=False):
        self.folder = folder
        self.imu_index = int(topic[topic.find('imu') + 3])
        self.topic = topic

        imuCsv = os.path.join(folder, "imu.csv")
        self.timestamps, self.accelData, self.gyroData = loadImuData(imuCsv)
        self.indices = np.arange(len(self.timestamps))
        if from_to:
            self.indices = self.truncateIndicesFromTime(self.indices, from_to)

    def __iter__(self):
        # Reset the bag reading
        return self.readDataset()

    def readDataset(self):
        return VimapImuCsvReaderIterator(self, self.indices)

    def truncateIndicesFromTime(self, indices, bag_from_to):
        timestampSecs = [time.toSec() for time in self.timestamps]
        bagstart = min(timestampSecs)
        baglength = max(timestampSecs) - bagstart

        # some value checking
        if bag_from_to[0] >= bag_from_to[1]:
            raise RuntimeError("Bag start time must be bigger than end time.".format(bag_from_to[0]))
        if bag_from_to[0] < 0.0:
            sm.logWarn("Bag start time of {0} s is smaller 0".format(bag_from_to[0]))
        if bag_from_to[1] > baglength:
            sm.logWarn("Bag end time of {0} s is bigger than the total length of {1} s".format(
                bag_from_to[1], baglength))

        # find the valid timestamps
        valid_indices = []
        for idx, timestamp in enumerate(timestampSecs):
            if bagstart + bag_from_to[0] <= timestamp <= bagstart + bag_from_to[1]:
                valid_indices.append(idx)
        sm.logWarn(
            "VimapImuCsvReader: truncated {0} / {1} IMU data.".format(len(indices) - len(valid_indices),
                                                                      len(indices)))
        return valid_indices

    def numMessages(self):
        return len(self.indices)

    def getMessage(self, idx):
        timestamp = self.timestamps[idx]
        omega = self.gyroData[idx]
        alpha = self.accelData[idx]
        return (timestamp, omega, alpha)


class Keypoint(object):
    def __init__(self, timestamp, vertex_idx, camera_idx, keypoint_idx, measurement, track_idx):
        self.timestamp = timestamp
        self.vertex_idx = vertex_idx
        self.camera_idx = camera_idx
        self.keypoint_idx = keypoint_idx
        self.measurement = measurement
        self.track_idx = track_idx

    def __str__(self):
        return "{:.9f} {} {} {} {} {}".format(self.timestamp.toSec(), self.vertex_idx,
                                              self.camera_idx, self.keypoint_idx,
                                              self.measurement, self.track_idx)


class Observation(object):
    def __init__(self, vertex_idx, camera_idx, keypoint_idx, landmark_idx):
        self.vertex_idx = vertex_idx
        self.camera_idx = camera_idx
        self.keypoint_idx = keypoint_idx
        self.landmark_idx = landmark_idx

    def __str__(self):
        return "{} {} {} {}".format(self.vertex_idx, self.camera_idx, self.keypoint_idx, self.landmark_idx)


class Vertex(object):
    def __init__(self, vertex_id, timestamp):
        self.vertex_idx = vertex_id
        self.timestamp = timestamp
        self._T_w_b = sm.Transformation()

    def set_T_w_b(self, qxyzw, pxyz):
        qxyzw[:3] = - qxyzw[:3]  # Hamilton to JPL convention.
        self._T_w_b = sm.Transformation(qxyzw, pxyz)

    def T_w_b(self):
        return self._T_w_b

    def __str__(self):
        return "{} {:.9f} {} {}".format(self.vertex_idx, self.timestamp.toSec(),
                                        ' '.join(map(str, self._T_w_b.t())),
                                        ' '.join(map(str, sm.quatInv(self._T_w_b.q()))))


class Landmark(object):
    def __init__(self, id, point):
        self.id = id
        self.position = point

    def __str__(self):
        return "{} {}".format(self.id, self.position)


def loadImuData(csvFile):
    timestamps = []  # list of acv.Time
    accelData = []  # list of numpy array (3,)
    gyroData = []  # list of numpy array (3,)
    with open(csvFile, 'r') as stream:
        csv_reader = csv.reader(stream , delimiter=',')

        for index, row in enumerate(csv_reader):
            if index == 0:
                continue
            timeString = row[0].strip()
            time = acv.Time(int(timeString[:-9]), int(timeString[-9:]))
            timestamps.append(time)
            accelData.append(np.array([float(row[1]), float(row[2]), float(row[3])]))
            gyroData.append(np.array([float(row[4]), float(row[5]), float(row[6])]))

    # print('Total imu data {}'.format(len(timestamps)))
    # print('Time(s), accel, gyro')
    # print("First row: {:.9f}, {}, {}".format(timestamps[0].toSec(), accelData[0], gyroData[0]))
    # print("Last row: {:.9f}, {}, {}".format(timestamps[-1].toSec(), accelData[-1], gyroData[-1]))
    return timestamps, accelData, gyroData


def loadTrackCsv(trackCsv):
    trackList = []
    maxCameraId = -1
    maxVertexId = -1
    with open(trackCsv, 'r') as stream:
        csv_reader = csv.reader(stream, delimiter=',')
        for index, row in enumerate(csv_reader):
            if index == 0:
                continue
            timeString = row[0].strip()
            time = acv.Time(int(timeString[:-9]), int(timeString[-9:]))
            cameraId = int(row[2])
            vertexId = int(row[1])
            keypoint = Keypoint(time, vertexId, cameraId, int(row[3]),
                                np.array([float(row[4]), float(row[5])]), int(row[8]))
            if vertexId > maxVertexId:
                maxVertexId = vertexId
            if cameraId > maxCameraId:
                maxCameraId = cameraId
            trackList.append(keypoint)
    # print('Total keypoints {}'.format(len(trackList)))
    # print('First keypoint {}'.format(trackList[0]))
    # print('Last keypoint {}'.format(trackList[-1]))
    return trackList, maxCameraId + 1, maxVertexId + 1


def loadObservationCsv(observationCsv):
    observationList = []
    with open(observationCsv, 'r') as stream:
        csv_reader = csv.reader(stream, delimiter=',')
        for index, row in enumerate(csv_reader):
            if index == 0:
                continue
            observation = Observation(int(row[0]), int(row[1]), int(row[2]), int(row[3]))
            observationList.append(observation)
    # print('Total observations {}'.format(len(observationList)))
    # print('First observation {}'.format(observationList[0]))
    # print('Last observation {}'.format(observationList[-1]))
    return observationList


def loadLandmarkCsv(landmarkCsv):
    landmarkList = []
    with open(landmarkCsv, 'r') as stream:
        csv_reader = csv.reader(stream, delimiter=',')
        for index, row in enumerate(csv_reader):
            if index == 0:
                continue
            landmark = Landmark(int(row[0]), np.array([float(row[1]), float(row[2]), float(row[3])]))
            landmarkList.append(landmark)
    # print('Total landmarks {}'.format(len(landmarkList)))
    # print('First landmark {}'.format(landmarkList[0]))
    # print('Last landmark {}'.format(landmarkList[-1]))
    return landmarkList


def loadVertexCsv(vertexCsv):
    vertexList = []
    with open(vertexCsv, 'r') as stream:
        csv_reader = csv.reader(stream, delimiter=',')
        for index, row in enumerate(csv_reader):
            if index == 0:
                continue
            timeString = row[1].strip()
            time = acv.Time(int(timeString[:-9]), int(timeString[-9:]))
            vertex = Vertex(int(row[0]), time)
            pxyz = np.array([float(row[2]), float(row[3]), float(row[4])])
            qxyzw = np.array([float(row[5]), float(row[6]), float(row[7]), float(row[8])])
            vertex.set_T_w_b(qxyzw, pxyz)
            vertexList.append(vertex)
    # print('Total vertices {}'.format(len(vertexList)))
    # print('First vertex {}'.format(vertexList[0]))
    # print('Last vertex {}'.format(vertexList[-1]))
    return vertexList
