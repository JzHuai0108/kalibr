# Appendix
This page details how to use features in this extended Kalibr version, in the order of practical value.
* a. calibrating the spatiotemporal parameters of a rolling shutter (RS) camera - IMU system, 
* b. simulating data of a rolling shutter camera - IMU system,
* c. calibrating a rolling shutter camera with optional IMU data,


## Install Kalibr
Follow instructions on [Kalibr installation](https://github.com/ethz-asl/kalibr/wiki/installation).

## Calibrate a (multiple) RS camera - IMU system
This runs very similarly to the default global shutter (GS) [camera - IMU calibration](https://github.com/ethz-asl/kalibr/wiki/camera-imu-calibration).
To calibrate a RS camera - IMU system, only two additional parameters are needed.
* add parameter *line_delay_nanoseconds* in the camera configuration yaml with an nonzero value, 
see [a template](./aslam_offline_calibration/kalibr/config_templates/camchain_template.yaml) for example.
* pass *--estimate-line-delay* to the *kalibr_calibrate_imu_camera* command.

You may try RS camera - IMU calibration with the 
[sample camera - IMU calibration dataset](https://github.com/ethz-asl/kalibr/wiki/downloads).
The final calibration result will include the refined line_delay_nanoseconds with a value very close to 0.

Unsurprisingly, aside from the RS effect, it is also possible to calibrate the IMU intrinsic parameters as in the original 
[GS camera - IMU calibration](https://github.com/ethz-asl/kalibr/wiki/Multi-IMU-and-IMU-intrinsic-calibration).

## Simulate RS camera - IMU data from real camera - IMU data

### Record an camera - IMU calibration dataset

The [sample camera - IMU dataset](https://github.com/ethz-asl/kalibr/wiki/downloads) also works.

### Run camera - IMU calibration and save the resulting B-spline models

To enable saving the B-spline model, pass the argument *--save-splines* to the *kalibr_calibrate_imu_camera* routine.
This will save the B spline models for the pose trajectory, gyro biases, and accelerometer biases, among others.

### Prepare camera and IMU configuration yamls for simulation

Refer to the [camera and IMU configuration yaml templates](./aslam_offline_calibration/kalibr/config_templates/) for examples.

### Based on the B-spline motion models, simulate RS camera and IMU data

Suppose output_dir is where the B-spline models are saved, simulate the RS camera - IMU with the below command.
```
kalibr_simulate_imu_camera $output_dir/bspline_pose.txt --cam camchain_imucam.yaml --imu imu.yaml \
  --target april_6x6.yaml --output_dir $output_dir
```

Note that currently the simulation only supports simulating for one camera.

The output files is in the [maplab csv dataset format](https://github.com/ethz-asl/maplab/wiki/CSV-Dataset-Format).

## Calibrate a RS camera with optional IMU data (very experimental)
The original RS camera calibration routine calibrates the RS effect and optional camera intrinsic parameters with only camera data.
Intuitively, its accuracy and stability can be boosted with the IMU data.
While inheriting the original functionality of kalibr_calibrate_rs_cameras,
the present implementation can take additional IMU data for constraints.
For simplicity, the calibrated IMU model is used where the IMU data are modeled with true values, biases, and noises, 
but not scale and misalignment.

To allow use of additional IMU data, 
* make sure the rosbag dataset has the IMU data,
* pass the IMU configuration yaml via *--imu* argument to *kalibr_calibrate_rs_cameras*.

## About covariance recovery
Theoretically, it is possible to recover the covariances for estimated parameters by using marginalization techniques.
However, the computation often takes too long so the covariance recovery functions in kalibr_calibrate_imu_camera and
kalibr_calibrate_rs_cameras are literally useless. 
The cause may be the strong coupling between adjacent control points in B-splines.

## A crash course on design variables in calibration with B-splines.

### kalibr\_calibrate\_imu_camera

* Design variables

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

* Access values
```
gravityDv: toEuclidean()
accelBiasDv/gyroBiasDv: spline().eval(t) or evalD(t, 0)
q_i_b_Dv: toRotationMatrix()
r_b_Dv: toEuclidean()
poseSplineDv: sm.Transformation(T_w_b.toTransformationMatrix()) where T_w_b=transformationAtTime(timeExpression, 0.0, 0.0)
```

* Error terms
```
CameraChainErrorTerms: error_t(frame, pidx, p) where error_t = self.camera.reprojectionErrorType + setMEstimatorPolicy

The realizations of reprojectionErrorType derive from the SimpleReprojectionError C++ class which is exported to python in exportReprojectionError().
The different error types are grouped into a variety of camera models in terms of python classes defined in 
aslam_cv/aslam_cv_backend_python/python/aslam_cv_backend/__init__.py.
These errors are independent of camera parameters.

AccelerometerErrorTerms: ket.EuclideanError + setMEstimatorPolicy
GyroscopeErrorTerms: ket.EuclideanError + setMEstimatorPolicy
Accel and gyro BiasMotionTerms: BSplineEuclideanMotionError
PoseMotionTerms: MarginalizationPriorErrorTerm (by default inactive)
```

### kalibr_calibrate_rs_cameras
This calibration procedure supports only one camera.

* Design variables
```
landmark_w_dv: aopt.HomogeneousPointDv (by default inactive)
__poseSpline_dv: asp.BSplinePoseDesignVariable
__camera_dv: The camera design variables are created by cameraModel.designVariable(self.geometry). 
The camera design variables are exported to python by exportCameraDesignVariables.
    projection: DesignVariableAdapter<projection_t>
    distortion: DesignVariableAdapter<distortion_t>
    shutter: DesignVariableAdapter<shutter_t>
```

* Error terms
```
For rolling shutter models, reprojectionErrorAdaptiveCovariance. These error types derive from the CovarianceReprojectionError C++ class, which is exported to python by exportCovarianceReprojectionError.
The different error types are grouped into a variety of camera models in terms of python classes defined in aslam_cv/aslam_cv_backend_python/python/aslam_cv_backend/__init__.py.
Reprojection errors with adaptive covariance is developed solely for rolling shutter cameras as discussed in
Section 3.5 Error Term Standardisation, Oth et. al., Rolling shutter camera calibration.
Because of the error standardisation, these reprojection errors depend on not only the camera parameters, 
but also the pose B splines.

For global shutter models, reprojectionError. These error types derive from the ReprojectionError C++ class which
is exported to python by exportReprojectionError. 
These errors depend on the camera parameters which may be optimized in the kalibr_calibrate_rs_cameras procedure.

regularizer: asp.BSplineMotionError of aslam_nonparametric_estimation/aslam_splines/include/aslam/backend.
```
