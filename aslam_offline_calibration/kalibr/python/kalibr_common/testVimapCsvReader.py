import os
import sm
from . import VimapCsvReader, FrameObservation
from . import VimapCsvWriter

vimapFolder = ""


def testLoadImuData():
    csv = os.path.join(vimapFolder, 'imu.csv')
    if os.path.isfile(csv):
        VimapCsvReader.loadImuData(csv)


def testLoadTrackCsv():
    csv = os.path.join(vimapFolder, 'tracks.csv')
    if os.path.isfile(csv):
        a, b, c = VimapCsvReader.loadTrackCsv(csv)
        print("{} {}".format(b, c))


def testLoadObservationCsv():
    csv = os.path.join(vimapFolder, 'observations.csv')
    if os.path.isfile(csv):
        a = VimapCsvReader.loadObservationCsv(csv)


def testLoadLandmarkCsv():
    csv = os.path.join(vimapFolder, 'landmarks.csv')
    if os.path.isfile(csv):
        a = VimapCsvReader.loadLandmarkCsv(csv)


def testLoadVertexCsv():
    csv = os.path.join(vimapFolder, 'vertices.csv')
    if os.path.isfile(csv):
        a = VimapCsvReader.loadVertexCsv(csv)


def testVimapImuCsvReader():
    if os.path.isdir(vimapFolder):
        dataset = VimapCsvReader.VimapImuCsvReader(vimapFolder, '/imu0', [2, 2.5], False)
        for timestamp, omega, alpha in dataset:
            print("{:.6f} {} {}".format(timestamp.toSec(), omega, alpha))


def testRemoveElements():
    a = [True, False, True, False, False, True, False]
    b = [1, 2, 3, 4, 5, 6, 7]
    c = [1, 3, 6]
    d = [e for i, e in enumerate(b) if a[i]]
    assert c == d


def testFrameObservation():
    print("sentinel {}".format(FrameObservation.landmarkSentinel()))


def testVimapCsvReader():
    if os.path.isdir(vimapFolder):
        dataset = VimapCsvReader.VimapCsvReader(vimapFolder, '/cam0/image_raw', [1, 2], False)
        targetObservations = dataset.getFeatureAssociations()
        print('Total frames {}'.format(dataset.numImages()))
        print('First frame {}'.format(targetObservations[0]))
        print('Last frame {}'.format(targetObservations[-1]))


def testPnPObservation():
    import aslam_cv as acv
    import numpy as np
    l = [np.random.rand(3), np.random.rand(3), np.random.rand(3)]
    o = [np.random.rand(2), np.random.rand(2), np.random.rand(2)]
    i = [1, 2, 3]

    stamp = acv.Time(100.0)
    T_t_c = sm.Transformation(np.array([0.5, 0.5, -0.5, 0.5]), np.array([4, 5, 6]))

    obs = acv.PnPObservation()
    a = obs.getCornersTargetFrame()
    b = obs.getCornersImageFrame()
    c = obs.getCornersIdx()
    assert len(a) == 0
    assert len(b) == 0
    assert len(c) == 0
    assert not obs.hasSuccessfulObservation()

    obs.setCornersTargetFrame(np.array(l))
    obs.setCornersImageFrame(np.array(o))
    obs.setCornersIdx(np.array(i))
    obs.setTime(stamp)
    obs.set_T_t_c(T_t_c)

    a = obs.getCornersTargetFrame()
    b = obs.getCornersImageFrame()
    c = obs.getCornersIdx()
    d = obs.time()
    e = obs.T_t_c()

    assert a.shape[0] == len(l) and a.shape[1] == 3
    assert b.shape[0] == len(o) and b.shape[1] == 2
    assert c.shape[0] == len(i)
    assert abs(d.toSec() - stamp.toSec()) < 1e-8
    residual = e.inverse() * T_t_c
    assert np.allclose(residual.t(), np.array([0, 0, 0]))
    assert np.allclose(residual.q(), np.array([0, 0, 0, 1]))
    assert obs.hasSuccessfulObservation()


def testFindNearestTimeSince():
    import aslam_cv as acv
    timeList = []
    timeList.append(acv.Time(10.0))
    timeList.append(acv.Time(10.1))
    timeList.append(acv.Time(10.2)),
    timeList.append(acv.Time(10.3))
    index = 0
    time = acv.Time(9.0)
    newindex = VimapCsvWriter.findNearestTimeSince(timeList, index, time)
    assert newindex == 0

    time = acv.Time(10.04)
    newindex = VimapCsvWriter.findNearestTimeSince(timeList, newindex, time)
    assert newindex == 0
    
    time = acv.Time(10.06)
    newindex = VimapCsvWriter.findNearestTimeSince(timeList, newindex, time)
    assert newindex == 1

    time = acv.Time(10.3)
    newindex = VimapCsvWriter.findNearestTimeSince(timeList, newindex, time)
    assert newindex == 3

    time = acv.Time(10.4)
    newindex = VimapCsvWriter.findNearestTimeSince(timeList, newindex, time)
    assert newindex == 3
