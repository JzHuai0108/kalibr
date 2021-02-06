import os
import sm
import VimapCsvReader

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
    print("sentinel {}".format(VimapCsvReader.FrameObservation.landmarkSentinel()))

def testVimapCsvReader():
    if os.path.isdir(vimapFolder):
        T_cN_imu = sm.Transformation()
        dataset = VimapCsvReader.VimapCsvReader(vimapFolder, '/cam0/image_raw', T_cN_imu, [1, 2], False)
        targetObservations = dataset.getFeatureAssociations()
        print('Total frames {}'.format(dataset.numImages()))
        print('First frame {}'.format(targetObservations[0]))
        print('Last frame {}'.format(targetObservations[-1]))

