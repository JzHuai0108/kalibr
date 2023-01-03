"""
Simulate monocular camera calibration data given poses, target, and the camera configuration yaml
"""
import sm

def parseArgs():
    # camera yaml
    # target yaml
    # pose mat
    # noise std

    pass

# References
# [1] https://github.com/JzHuai0108/kalibr/blob/develop/aslam_offline_calibration/kalibr/python/kalibr_imu_camera_calibration/Simulator.py


def toSmTransformation(qxyzw, pxyz):
    qxyzw[:3] = - qxyzw[:3]  # Hamilton to JPL convention.
    return sm.Transformation(qxyzw, pxyz)


def main():
    args = parseArgs() # TODO: binliang refer to https://github.com/JzHuai0108/vio_common/blob/master/python/rgbd_bag_to_synced_images.py#L56-L81
    cam = loadCamera(args.camera_yaml) # TODO: refer to https://github.com/JzHuai0108/kalibr/blob/develop/aslam_offline_calibration/kalibr/python/kalibr_imu_camera_calibration/Simulator.py#L93-L118
    target = loadTarget(args.target_yaml) # TODO: refer to https://github.com/JzHuai0108/kalibr/blob/develop/aslam_offline_calibration/kalibr/python/kalibr_imu_camera_calibration/Simulator.py#L120-L168
    poses = loadPoses(args.posemat) # TODO: load matlab mat file of poses by scipy.io.loadmat
    # Note convert the pose to aslam::Transformation() by toSmTransformation
    # Note the cspond index is in matlab format, starts from 1.

    numLandmarks = targetObservation.getTotalTargetPoint()
    imageWidth = resolution[0]
    imageHeight = resolution[1]

    numFailedProjection = 0
    x = []
    cspond = []
    for sm_T_w_c in poses:
        for iota in range(numLandmarks):
            validProjection, imagePoint = targetObservation.projectATargetPoint(cam, sm_T_w_c,
                                                                                         iota) # 3x1.

            if not validProjection:
                numFailedProjection += 1
                continue
            # add noise # TODO refer to https://github.com/JzHuai0108/kalibr/blob/develop/aslam_offline_calibration/kalibr/python/kalibr_imu_camera_calibration/Simulator.py#L359-L362
            x.append(imagePoint)

    # check the project points are close to the values in the mat file, x

    saveMat(boardmat) # refer to https://github.com/castacks/tartancalib/blob/main/aslam_offline_calibration/kalibr/python/tartan_calibrate#L586-L593
    saveMat(cornersmat)


if __name__ == "__main__":
    main()