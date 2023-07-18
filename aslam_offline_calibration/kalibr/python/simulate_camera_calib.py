"""
Simulate 2d observations for monocular camera calibration given poses, target, and the camera configuration yaml.
No rolling shutter effect is added since camera calibration can do away with it by holding the target static.
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
    # camera yaml, target yaml, pose mat, noise std
    parser = argparse.ArgumentParser(
        description="input: include camera.yaml, target.yaml, pose.mat and noise std")
    parser.add_argument("--camera-yaml",
                        default='./data/camera.yaml',
                        help="camera.yaml")
    parser.add_argument("--target-yaml",
                        default='./data/target.yaml',
                        help="target.yaml")
    parser.add_argument("--posemat",
                        default="corners.mat", 
                        help=("corners saved in babelcalib's matlab format by using the tartancalib tool at"
                              "git@github.com:JzHuai0108/tartancalib.git"))
    parser.add_argument('--noise-std',
                        type=float,
                        default=0.01,
                        help='noise std')
    parser.add_argument("--outputmat",
                        default='./data/outcorner.mat',
                        help="outcorner.mat")

    parser.add_argument('--sim-only-used', action='store_true', help='only simulate used poses')
    args = parser.parse_args()

    return args

def printExtraCameraDetails(camConfig):
    resolution = camConfig.getResolution()
    print('  Camera resolution: {}'.format(resolution))
    imageNoise = camConfig.getImageNoise()
    print('  Image noise std dev: {}'.format(imageNoise))
    lineDelay = camConfig.getLineDelayNanos()
    print("  Line delay: {} ns".format(lineDelay))
    updateRate = camConfig.getUpdateRate()
    print("  Update rate: {} Hz".format(updateRate))

def isKB(modelpair):
    return modelpair[0] == "pinhole" and modelpair[1] == "equidistant"

def loadCamera(camera_yaml):
    print("Camera chain from {}".format(camera_yaml))
    chain = kc.CameraChainParameters(camera_yaml)
    camGeometryList = []
    cammodels = []
    numCameras = chain.numCameras()
    for i in range(numCameras):
        camConfig = chain.getCameraParameters(i)
        camConfig.printDetails()
        camera_model, intrinsics = camConfig.getIntrinsics()
        dist_model, dist_coeff = camConfig.getDistortion()
        camera = kc.AslamCamera.fromParameters(camConfig)
        camGeometryList.append(camera.geometry)
        cammodels.append((camera_model, dist_model))
    return camGeometryList, cammodels

def loadTarget(target_yaml):
    targetConfig = kc.CalibrationTargetParameters(target_yaml)
    print("Target used in the simulation:")
    targetConfig.printDetails()
    targetObservation = None
    allTargetCorners = None

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
    '''Note the cspond index is in matlab format, starts from 1.'''
    data = (io.loadmat(posemat))['corners']
    poses = []
    x = []
    cspond = []
    gused = []
    for i in range(data.shape[1]):
        ptemp = data[0, i]['t_T_c'][0, 0][0, :3]
        qtemp = data[0, i]['t_T_c'][0, 0][0, 3:]
        poses.append(toSmTransformation(qtemp, ptemp))
        x.append(data[0, i]['x'][0, 0])
        cspond.append(data[0, i]['cspond'][0, 0])
        gused.append(data[0, i]['used'][0, 0])

    resolution = (io.loadmat(posemat))['imgsize'][0]
    times = (io.loadmat(posemat))['times'][0]

    return poses, resolution, x, cspond, gused, times

def saveMat(corners_mat, imgsize, times, used, marked, outputmat):
    io.savemat(outputmat, {"corners": corners_mat, "imgsize": imgsize,
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
    args = parseArgs() 
    cam, models = loadCamera(args.camera_yaml) 
    targetObservation = loadTarget(args.target_yaml) 
    poses, resolution, origx, origcspond, origused, origtimes = loadPoses(args.posemat)

    numLandmarks = targetObservation.getTotalTargetPoint()
    imageWidth = resolution[0]
    imageHeight = resolution[1]

    numFailedProjection = 0
    corners_mat = []
    times = []
    numusedframes = 0
    numcheckedpoints = 0
    iskb = isKB(models[0])
    for j, pose in enumerate(poses):
        if args.sim_only_used and origused[j] == 0:
            continue
        x = []
        cspond = []
        sm_T_w_c = poses[j]
        for iota in range(numLandmarks):
            validProjection, imagePoint = targetObservation.projectATargetPoint(cam[0], sm_T_w_c, iota, iskb)
            # imagePoint 3x1
            if not validProjection:
                numFailedProjection += 1
                continue
            xnoise = gauss(0.0, args.noise_std)
            ynoise = gauss(0.0, args.noise_std)
            noisyPoint = [noisyValue(imagePoint[0, 0], imageWidth, xnoise),
                            noisyValue(imagePoint[1, 0], imageHeight, ynoise)]
            x.append(noisyPoint)
            spd = origcspond[j][0]
            cspond.append([iota+1, 1])
            found = [i for i, y in enumerate(spd) if y == iota + 1]
            if len(found) == 0:
                continue
            m = int(found[0])
            origPoint = origx[j][:, m]
            subPoint = [noisyPoint[0]-origPoint[0], noisyPoint[1]-origPoint[1]]
            numcheckedpoints += 1
            d = np.linalg.norm(subPoint)
            if d > 5:
                print('Warn: Dist(noisyPoint({}) - origPoint({})) {} > 5'.format(noisyPoint, origPoint, d))

        np_T_tc = np.zeros(7)
        np_T_tc[0:3] = sm_T_w_c.t()
        # quatInv converts JPL quaternion to Hamilton quaternion (x,y,z,w).
        np_T_tc[3:7] = sm.quatInv(sm_T_w_c.q())
        x = np.array(x).transpose()
        cspond = np.array(cspond).transpose()
        corners_mat.append({"x": x, "cspond": cspond, 't_T_c': np_T_tc, 'used' : 1})
        times.append(origtimes[j])
        numusedframes += 1
    print("Simulated {} frames, checked {} points.".format(numusedframes, numcheckedpoints))
    saveMat(corners_mat, resolution, times, numusedframes, numusedframes, args.outputmat)


if __name__ == "__main__":
    main()