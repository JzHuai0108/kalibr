"""
Simulate monocular camera calibration data given poses, target, and the camera configuration yaml
"""
import sm
import kalibr_common as kc
import aslam_cv as acv
import aslam_cameras_april as acv_april

import argparse
from scipy import io
import numpy as np
from random import gauss

def parseArgs():
    # camera yaml
    # target yaml
    # pose mat
    # noise std
    parser = argparse.ArgumentParser(
        description="input: include camera.yaml, target.yaml, pose.mat and noise std")
    parser.add_argument("--camera_yaml",
                        default='[./data/camera.yaml]',
                        help="camera.yaml")
    parser.add_argument("--target_yaml",
                        default='[./data/target.yaml]',
                        help="target.yaml")
    parser.add_argument("--posemat",
                        default='[./data/pose.mat]',
                        help="pose.mat")
    parser.add_argument('--noise_std',
                        type=float,
                        default=[0.01],
                        help='noise std')
    parser.add_argument("--out_cornerfile",
                        default='[./data/outcorner.mat]',
                        help="outcorner.mat")
    args = parser.parse_args()

    return args

# # References
# # [1] https://github.com/JzHuai0108/kalibr/blob/develop/aslam_offline_calibration/kalibr/python/kalibr_imu_camera_calibration/Simulator.py

def printExtraCameraDetails(camConfig):
    resolution = camConfig.getResolution()
    print('  Camera resolution: {}'.format(resolution))
    imageNoise = camConfig.getImageNoise()
    print('  Image noise std dev: {}'.format(imageNoise))
    lineDelay = camConfig.getLineDelayNanos()
    print("  Line delay: {} ns".format(lineDelay))
    updateRate = camConfig.getUpdateRate()
    print("  Update rate: {} Hz".format(updateRate))

def loadCamera(camera_yaml):
    print("Camera chain from {}".format(camera_yaml))
    chain = kc.CameraChainParameters(camera_yaml)
    T_imu_cam_list = []
    timeOffsetList = []
    camGeometryList = []
    numCameras = chain.numCameras()
    for i in range(numCameras):
        camConfig = chain.getCameraParameters(i)
        camConfig.printDetails()
        # printExtraCameraDetails(camConfig)
        # These parameters are set to default values assuming no IMU is present.
        T_imu_cam_list.append(sm.Transformation())
        timeOffsetList.append(0)
        camera = kc.AslamCamera.fromParameters(camConfig)
        camGeometryList.append(camera.geometry)
    
    return camGeometryList

def loadTarget(target_yaml):
    targetConfig = kc.CalibrationTargetParameters(target_yaml)
    print("Target used in the simulation:")
    targetConfig.printDetails()
    targetObservation = None
    allTargetCorners = None
    # setupCalibrationTarget(targetConfig, showExtraction=False, showReproj=False, imageStepping=False)

    # load the calibration target configuration
    targetParams = targetConfig.getTargetParams()
    targetType = targetConfig.getTargetType()

    if targetType == 'checkerboard':
        options = acv.CheckerboardOptions()
        options.filterQuads = True
        options.normalizeImage = True
        options.useAdaptiveThreshold = True
        options.performFastCheck = False
        options.windowWidth = 5
        options.showExtractionVideo = False
        grid = acv.GridCalibrationTargetCheckerboard(targetParams['targetRows'],
                                                        targetParams['targetCols'],
                                                        targetParams['rowSpacingMeters'],
                                                        targetParams['colSpacingMeters'],
                                                        options)
    elif targetType == 'circlegrid':
        options = acv.CirclegridOptions()
        options.showExtractionVideo = False
        options.useAsymmetricCirclegrid = targetParams['asymmetricGrid']
        grid = acv.GridCalibrationTargetCirclegrid(targetParams['targetRows'],
                                                    targetParams['targetCols'],
                                                    targetParams['spacingMeters'],
                                                    options)
    elif targetType == 'aprilgrid':
        options = acv_april.AprilgridOptions()
        options.showExtractionVideo = False
        options.minTagsForValidObs = int(np.max([targetParams['tagRows'], targetParams['tagCols']]) + 1)

        grid = acv_april.GridCalibrationTargetAprilgrid(targetParams['tagRows'],
                                                        targetParams['tagCols'],
                                                        targetParams['tagSize'],
                                                        targetParams['tagSpacing'],
                                                        options)
    else:
        raise RuntimeError("Unknown calibration target.")

    options = acv.GridDetectorOptions()
    options.imageStepping = False
    options.plotCornerReprojection = False
    options.filterCornerOutliers = True

    targetObservation = acv.GridCalibrationTargetObservation(grid)
    allTargetCorners = targetObservation.getAllCornersTargetFrame()  # nx3
    assert allTargetCorners.shape[0] == targetObservation.getTotalTargetPoint()

    return targetObservation

def loadPoses(posemat):
    data = (io.loadmat(posemat))['corners']
    poses = []
    x = []
    cspond = []
    gused = []
    print(data.shape[1])
    for i in range(data.shape[1]):
        ptemp = data[0, i]['t_T_c'][0, 0][0, :3]
        qtemp = data[0, i]['t_T_c'][0, 0][0, 3:]
        # Note convert the pose to aslam::Transformation() by toSmTransformation
        poses.append(toSmTransformation(qtemp, ptemp))
        x.append(data[0, i]['x'][0, 0])
        cspond.append(data[0, i]['cspond'][0, 0])
        gused.append(data[0, i]['used'][0, 0])

    # Note the cspond index is in matlab format, starts from 1.
    resolution = (io.loadmat(posemat))['imgsize'][0]
    # print(resolution)
    times = (io.loadmat(posemat))['times'][0]

    return poses, resolution, x, cspond, gused, times

def saveMat(corners_mat, imgsize, times, used, marked, out_cornerfile):
    # # save board def mat
    # io.savemat('board.mat',{"boards":{"Rt":Rt,"X":x_board}})s
    io.savemat(out_cornerfile, {"corners": corners_mat, "imgsize": imgsize,
                                "times": times,
                                'used':used, 'marked': marked})

def toSmTransformation(qxyzw, pxyz):
    qxyzw[:3] = - qxyzw[:3]  # Hamilton to JPL convention.
    return sm.Transformation(qxyzw, pxyz)

def noisyValue(x, upperbound, noise):
    if x <= 1 or x >= upperbound - 1:
        noisyx = x
    elif x + noise < 0:
        noisyx = x - noise
    elif x + noise > upperbound:
        noisyx = x - noise
    else:
        noisyx = x + noise
    return noisyx


def main():
    args = parseArgs() # TODO: binliang refer to https://github.com/JzHuai0108/vio_common/blob/master/python/rgbd_bag_to_synced_images.py#L56-L81
    cam = loadCamera(args.camera_yaml) # TODO: refer to https://github.com/JzHuai0108/kalibr/blob/develop/aslam_offline_calibration/kalibr/python/kalibr_imu_camera_calibration/Simulator.py#L93-L118
    targetObservation = loadTarget(args.target_yaml) # TODO: refer to https://github.com/JzHuai0108/kalibr/blob/develop/aslam_offline_calibration/kalibr/python/kalibr_imu_camera_calibration/Simulator.py#L120-L168
    poses, resolution, origx, origcspond, origused, times  = loadPoses(args.posemat) # TODO: load matlab mat file of poses by scipy.io.loadmat

    numLandmarks = targetObservation.getTotalTargetPoint()
    imageWidth = resolution[0]
    imageHeight = resolution[1]

    numFailedProjection = 0
    corners_mat = []
    print(numLandmarks)
    for j in range(len(poses)):
        x = []
        cspond = []
        sm_T_w_c = poses[j]
        for iota in range(numLandmarks):
            validProjection, imagePoint = targetObservation.projectATargetPoint(cam[0], sm_T_w_c,
                                                                                         iota) # 3x1.
            if not validProjection:
                numFailedProjection += 1
                continue
            # add noise # TODO refer to https://github.com/JzHuai0108/kalibr/blob/develop/aslam_offline_calibration/kalibr/python/kalibr_imu_camera_calibration/Simulator.py#L359-L362
            xnoise = gauss(0.0, args.noise_std[0])
            ynoise = gauss(0.0, args.noise_std[0])
            noisyPoint = [noisyValue(imagePoint[0, 0], imageWidth, xnoise),
                            noisyValue(imagePoint[1, 0], imageHeight, ynoise)]
            x.append(noisyPoint)
            spd = origcspond[j][0]
            cspond = origcspond[j]
            m = np.where(spd==iota+1)
            if m[0].shape[0]==0:
                continue
            origPoint = [origx[j][:, m][0][0][0],origx[j][:, m][1][0][0]]
            subPoint = [noisyPoint[0]-origPoint[0], noisyPoint[1]-origPoint[1]]
            if np.linalg.norm(subPoint)>5:
                print('warn: Dist(noisyPoint-origPoint)>5')
                print(noisyPoint)
                print(origPoint)
        
        np_T_tc = np.zeros(7)
        np_T_tc[0:3] = sm_T_w_c.t()
        # quatInv converts JPL quaternion to Hamilton quaternion (x,y,z,w).
        np_T_tc[3:7] = sm.quatInv(sm_T_w_c.q())
        x = list(map(list, zip(*x)))
        corners_mat.append({"x": x, "cspond": cspond, 't_T_c': np_T_tc, 'used' : origused[j]})

    # refer to https://github.com/castacks/tartancalib/blob/main/aslam_offline_calibration/kalibr/python/tartan_calibrate#L586-L593
    saveMat(corners_mat, resolution, times, 32, 32, args.out_cornerfile)

if __name__ == "__main__":
    main()