# Install Kalibr from source on Ubuntu 16.04 with ROS Kinetic

## Install dependencies as instructed on the Wiki of Kalibr on github
```
sudo apt-get install python-setuptools python-rosinstall ipython libeigen3-dev libboost-all-dev doxygen libopencv-dev ros-kinetic-vision-opencv ros-kinetic-image-transport-plugins ros-kinetic-cmake-modules python-software-properties software-properties-common libpoco-dev python-matplotlib python-scipy python-git python-pip ipython libtbb-dev libblas-dev liblapack-dev python-catkin-tools libv4l-dev 
```
## Install the additional python packages
```
sudo apt-get install python-wxversion python-wxtools
```
## Install python-igraph via
```
sudo pip install python-igraph --upgrade
```
In installing python-igraph, there may arise an error reads as in the below reference
reference: https://stackoverflow.com/questions/37495375/python-pip-install-throws-typeerror-unsupported-operand-types-for-retry
Even using the solution suggested in the reference page, the error in installing python-igraph persists.

Don't fuss about it. Simply ignore it and go ahead to build kalibr with catkin_make. The building procedure succeeds in my case.

# Simulate monocular inertial data

1. Record an imu camera calibration dataset and generate a B-spline model

The dataset can be the imu camera calibration sample data provided by Kalibr

```
kalibr_calibrate_imu_camera --cam camchain.yaml --target april_6x6.yaml --imu imu_adis16448.yaml --bag dynamic.bag \
  --bag-from-to 5 45 --dont-show-report

```
The output results will include the B-spline model, knotCoeffT, and ref_pose, ref_state, ref_imu_meas generated from the model. More details refer to saveBSpline in python/kalibr_imu_camera_calibration/IccUtil.py.

2. Based on the B-spline model, simulate rolling shutter camera measurements.
Templates for the camera-IMU system, and the IMU is at the config_templates folder.
```
kalibr_simulate_imu_camera $output_dir/bspline_pose.txt --cam $script_dir/camchain_template.yaml --imu $script_dir/imu_template.yaml \
  --target $data_dir/april_6x6.yaml --output_dir $output_dir
```

# A crash course on design variables of calibration by using B splines.

## kalibr\_calibrate\_imu_camera

Design variables

```
poseDv: asp.BSplinePoseDesignVariable
gravityDv: aopt.EuclideanPointDv or aopt.EuclideanDirection
imu Dvs: 
    gyroBiasDv: asp.EuclideanBSplineDesignVariable
    accelBiasDv: asp.EuclideanBSplineDesignVariable
    q_i_b_Dv: aopt.RotationQuaternionDv (can be inactive)
    r_b_Dv: aopt.EuclideanPointDv (can be inactive)
    potential imu Dvs for scaledMisalignment: q_gyro_i_Dv, M_accel_Dv, M_gyro_Dv, M_accel_gryo_Dv
    potential imu Dvs for size effect: rx_i_Dv, ry_i_Dv, rz_i_Dv, Ix_Dv, Iy_Dv, Iz_Dv
camera Dvs: 
    T_c_b_Dv: aopt.TransformationDv (T_cNplus1_cN)
    cameraTimetoImuTimeDv: aopt.Scalar
```

Access content 
```
gravityDv: toEuclidean()
accelBiasDv/gyroBiasDv: spline().eval(t) or evalD(t, 0)
q_i_b_Dv: toRotationMatrix() 
r_b_Dv: toEuclidean()
poseSplineDv: sm.Transformation(T_w_b.toTransformationMatrix()) where T_w_b=transformationAtTime(timeExpression, 0.0, 0.0)
```

Error terms
```
CameraChainErrorTerms: error_t(frame, pidx, p) where error_t = self.camera.reprojectionErrorType + setMEstimatorPolicy

The realizations of reprojectionErrorType derive from the SimpleReprojectionError C++ class which is exported to python in exportReprojectionError().
The different error types are grouped into a variety of camera models in terms of python classes defined in aslam_cv/aslam_cv_backend_python/python/aslam_cv_backend/__init__.py.
These errors are independent of camera parameters, thus simple.

AccelerometerErrorTerms: ket.EuclideanError + setMEstimatorPolicy
GyroscopeErrorTerms: ket.EuclideanError + setMEstimatorPolicy
Accel and gyro BiasMotionTerms: BSplineEuclideanMotionError
PoseMotionTerms: MarginalizationPriorErrorTerm (by default inactive)
```

## kalibr_calibrate_rs_cameras
This calibration procedure supports only one camera.

Design variables
```
landmark_w_dv: aopt.HomogeneousPointDv (by default inactive)
__poseSpline_dv: asp.BSplinePoseDesignVariable
__camera_dv:
    projection: DesignVariableAdapter<projection_t>
    distortion: DesignVariableAdapter<distortion_t>
    shutter: DesignVariableAdapter<shutter_t>
```

Error terms
```
For rolling shutter models, reprojectionErrorAdaptiveCovariance. These error types derive from the CovarianceReprojectionError C++ class, which is exported to python by exportCovarianceReprojectionError. 
The different error types are grouped into a variety of camera models in terms of python classes defined in aslam_cv/aslam_cv_backend_python/python/aslam_cv_backend/__init__.py.
Reprojection errors with adaptive covariance is developed solely for rolling shutter cameras as discussed in
3.5. Error Term Standardisation of Oth et. al. CVPR Rolling shutter camera calibration.
Because of the error standardisation, these reprojection errors depend on not only the camera parameters, 
but also the pose B splines.

For global shutter models, reprojectionError. These error types derive from the ReprojectionError C++ class which is exported to python by exportReprojectionError. These errors depend on the camera parameters which may be optimized in the kalibr_calibrate_rs_cameras procedure.

regularizer: asp.BSplineMotionError of aslam_nonparametric_estimation/aslam_splines/include/aslam/backend.
to disambiguate, there is another BSplineMotionError in aslam_cv/aslam_cv_error_terms/include/aslam/backend.
The two implementations are more or less the same, the one in aslam_splines looks newer.
```

