![Kalibr](https://raw.githubusercontent.com/wiki/ethz-asl/kalibr/images/kalibr_small.png)

- [Introduction](#introduction)
- [Installation](#installation)
  - [Ubuntu 18.04 + ROS1 melodic](#ubuntu-1804--ros1-melodic)
  - [Ubuntu 20.04 + ROS1 noetic](#ubuntu-2004--ros1-noetic)
- [RS camera-IMU calibration](#rs-camera-imu-calibration)
- [Simulate RS camera-IMU data from real data](#simulate-rs-camera-imu-data-from-real-data)
  - [1. Find or record a camera-IMU calibration dataset](#1-find-or-record-a-camera-imu-calibration-dataset)
  - [2. Run camera-IMU calibration and save the trajectory](#2-run-camera-imu-calibration-and-save-the-trajectory)
  - [3. Prepare camera and IMU configuration yamls for simulation](#3-prepare-camera-and-imu-configuration-yamls-for-simulation)
  - [4. Simulate RS camera and IMU data](#4-simulate-rs-camera-and-imu-data)
- [IMU noise identification](#imu-noise-identification)
- [Static frame detection](#static-frame-detection)
- [RS camera calibration with IMU data (experimental)](#rs-camera-calibration-with-imu-data-experimental)
- [Docker](#docker)
- [A crash course on calibration with B-splines](#a-crash-course-on-calibration-with-b-splines)
  - [kalibr_calibrate_imu_camera](#kalibr_calibrate_imu_camera)
    - [Design variables](#design-variables)
    - [Get values of design variables](#get-values-of-design-variables)
    - [Error terms](#error-terms)
  - [kalibr_calibrate_rs_cameras](#kalibr_calibrate_rs_cameras)
    - [Design variables](#design-variables-1)
    - [Error terms](#error-terms-1)
- [Citing](#citing)
- [Compiled rolling shutter parameters for consumer products](#compiled-rolling-shutter-parameters-for-consumer-products)
- [TODOs](#todos)

## Introduction
Kalibr is a popular and excellent calibration toolbox based on continuous-time B-splines
for calibrating cameras and camera-IMU systems.

This fork extends the [original Kalibr](https://github.com/ethz-asl/kalibr) developed by ethz-asl 
to support rolling shutter (RS) camera-IMU calibration, IMU noise identification, and 
static frame detection.
Moreover, this fork supports Ubuntu 20.04 by assimilating the [ori-drs fork](https://github.com/ori-drs/kalibr).

Summarily, this package solves the following problems:
1. **Rolling shutter camera-IMU calibration**:
   extrinsic and temporal calibration of a rolling shutter camera relative to an rigidly attached IMU
2. **Simulation from real data for a mono RS camera-IMU system**

3. **IMU noise identification within camera-IMU calibration**:
   identify the IMU noise parameters for camera-IMU calibration by using the same calibration data 
   so as to avoid the long Allan variance analysis

4. **Static frame detection**:
   detect static frames by checking the optic flow of target corners to 
   simplify the data collection procedure for rolling shutter camera calibration

in addition to the below calibration problem solved by the original Kalibr package:
>1. **Multiple camera calibration**: 
>    intrinsic and extrinsic calibration of a camera-systems with non-globally shared overlapping fields of view
>2. **Visual-inertial calibration calibration (camera-IMU)**:
>    spatial and temporal calibration of an IMU w.r.t a camera-system
>3. **Rolling shutter camera calibration**:
>    full intrinsic calibration (projection, distortion and shutter parameters) of rolling shutter cameras
The [wiki pages](https://github.com/ethz-asl/kalibr/wiki) of the original Kalibr provide clear and detailed explanations on how to 
perform the latter three tasks.
Here we focus on the former three tasks.

## Installation

### Ubuntu 18.04 + ROS1 melodic
For Ubuntu <=18.04 + ROS1 <= melodic, follow instructions at [here](https://github.com/ethz-asl/kalibr/wiki/installation).
In addition, install suitesparse by
```
sudo apt-get install libsuitesparse-dev
```
because this fork uses the system wide suitesparse whereas the original Kalibr builds suitesparse in the Kalibr workspace.
This installation procedure can be greatly simplified by using the provided Dockerfiles, see the [Docker](#docker) section.

### Ubuntu 20.04 + ROS1 noetic

```
sudo apt update
sudo apt-get install python3-setuptools python3-rosinstall ipython3 libeigen3-dev libboost-all-dev doxygen libopencv-dev \
ros-noetic-vision-opencv ros-noetic-image-transport-plugins ros-noetic-cmake-modules python3-software-properties \
software-properties-common libpoco-dev python3-matplotlib python3-scipy python3-git python3-pip libtbb-dev libblas-dev \
liblapack-dev libv4l-dev python3-catkin-tools python3-igraph libsuitesparse-dev

pip3 install wxPython
```
If you encounter errors like "E: Unable to locate package python3-catkin-tools",
then setup the sources.list and keys as instructed [here](http://wiki.ros.org/Installation/Ubuntu).

```
mkdir ~/kalibr_ws/src
cd ~/kalibr_ws/src
git clone --recursive https://github.com/JzHuai0108/kalibr

cd ~/kalibr_ws
source /opt/ros/noetic/setup.bash
catkin init
catkin config --extend /opt/ros/noetic
catkin config --merge-devel
catkin config --cmake-args -DCMAKE_BUILD_TYPE=Release

catkin build -DCMAKE_BUILD_TYPE=Release -j4
```

## RS camera-IMU calibration
To calibrate a RS camera-IMU system, only two additional parameters are needed compared to the default global shutter (GS) [camera-IMU calibration](https://github.com/ethz-asl/kalibr/wiki/camera-imu-calibration).
* add parameter *line_delay_nanoseconds* in the camera configuration yaml with an nonzero value, 
see [a template](./aslam_offline_calibration/kalibr/config_templates/camchain_template.yaml) for example.
* pass *--estimate-line-delay* to the *kalibr_calibrate_imu_camera* command.

Let's try out RS camera-IMU calibration with the 
[kalibr dynamic sample data](https://drive.google.com/file/d/0B0T1sizOvRsUcGpTWUNTRC14RzA/edit?usp=sharing) 
which was collected by a GS camera.
```
source ~/kalibr_ws/devel/setup.bash
cd /path/to/kalibr_dynamic_sample

echo "Calibrate imu-camera system with RS model and estimate-line-delay enabled..."
sed -i "/line_delay_nanoseconds/c\  line_delay_nanoseconds: 5000" camchain.yaml
mkdir -p rs_est
cd rs_est
kalibr_calibrate_imu_camera --target april_6x6.yaml --cam camchain.yaml --imu imu_adis16448.yaml \
 --bag dynamic.bag --bag-from-to 5 45 --estimate-line-delay --dont-show-report
```
In the final report, you should find the estimated line delay is less than 200 nanoseconds, 
confirming that the dataset was collected by a GS camera.

## Simulate RS camera-IMU data from real data

### 1. Find or record a camera-IMU calibration dataset
For instance, you may use the [sample camera-IMU dataset](https://github.com/ethz-asl/kalibr/wiki/downloads).

### 2. Run camera-IMU calibration and save the trajectory
To enable saving the trajectory represented by B-spline models, 
pass the argument *--save-splines* to the *kalibr_calibrate_imu_camera* routine.
This will save the B spline models for the pose trajectory, gyro biases, and accelerometer biases, among others.

### 3. Prepare camera and IMU configuration yamls for simulation
Refer to the [camera and IMU configuration yaml templates](./aslam_offline_calibration/kalibr/config_templates/) for examples.

### 4. Simulate RS camera and IMU data
Suppose $output_dir is where the B-spline models are saved, simulate the RS camera-IMU data with.
```
kalibr_simulate_imu_camera $output_dir/bspline_pose.txt --cam camchain_imucam.yaml --imu imu.yaml \
  --target april_6x6.yaml --output_dir $output_dir
```

Note that currently the simulation only supports simulating for one camera.

The output files are in the [maplab csv dataset format](https://github.com/ethz-asl/maplab/wiki/CSV-Dataset-Format).

## IMU noise identification
Refer to [test_on_dynamic_sample](ci/test_on_dynamic_sample.sh) for shell scripts on IMU noise identification. More instructions will be added soon.

## Static frame detection
Refer to [extract_static_frames](aslam_offline_calibration/kalibr/python/extract_static_frames) for more info. More instructions will be added soon.

## RS camera calibration with IMU data (experimental)
The original RS camera calibration routine calibrates the RS effect and optional camera intrinsic parameters with only camera data.
Intuitively, its accuracy and stability can be boosted with the IMU data.
While inheriting the original functionality of kalibr_calibrate_rs_cameras,
the present implementation can take additional IMU data for constraints.
For simplicity, the calibrated IMU model is used where the IMU data are modeled with true values, biases, and noises, 
excluding scale factor errors and misalignment.

To allow use of additional IMU data, 
* make sure that the rosbag dataset has the IMU data,
* pass the IMU configuration yaml via *--imu* argument to *kalibr_calibrate_rs_cameras*.

## Docker
On Ubuntu, we can use kalibr from a docker container without installing its dependencies system-wide (e.g., ROS1) on the host computer.
Of course, this requires the [docker engine](https://docs.docker.com/engine/install/ubuntu/) installed on the host computer.

To start a docker container where the kalibr can run, we need a docker image which is its blueprint.
The docker image can be created by building the provided dockerfile.
The below commands assume that the ROS1 melodic docker image will be created and used, 
replace melodic with noetic if a ROS1 noetic docker image is desired.
```
cd kalibr/docker/melodic
chmod +x build.sh
./build.sh
```
To confirm that the docker image is created successfully, list the existing docker images by
```
docker images
```
Then we can start a docker container from the image and run Kalibr.

```
cd kalibr/docker
chmod +x run.sh
./run.sh <folder-of-rosbag> melodic 1

```
Note that folder-of-rosbag is the folder on the host computer containing the data bag for calibration.
The run.sh step mounts folder-of-rosbag to /root/data in the docker container, and 
opens an interactive shell session with the Kalibr workspace loaded in the container.
Then in the interactive shell, you may run the calibration commands like
```
kalibr_calibrate_imu_camera --target april_6x6.yaml --cam camchain.yaml --imu imu_adis16448.yaml \
 --bag dynamic.bag --bag-from-to 5 45 --estimate-line-delay --dont-show-report
```
For some tasks, you may want to open extra shell sessions for the container by 
```
docker exec -it <CONTAINER ID> /bin/bash
# CONTAINER ID can be found by 
docker ps
```

## A crash course on calibration with B-splines
For those interested in developing with kalibr, the building blocks for kalibr_calibrate_imu_camera and kalibr_calibrate_rs_cameras
are described below.

### kalibr_calibrate_imu_camera

#### Design variables
poseDv: asp.BSplinePoseDesignVariable
gravityDv: aopt.EuclideanPointDv or aopt.EuclideanDirection
imu Dvs: 
    gyroBiasDv: asp.EuclideanBSplineDesignVariable
    accelBiasDv: asp.EuclideanBSplineDesignVariable
    q_i_b_Dv: aopt.RotationQuaternionDv (can be inactive)
    r_b_Dv: aopt.EuclideanPointDv (can be inactive)
    optional imu Dvs for scaledMisalignment: q_gyro_i_Dv, M_accel_Dv, M_gyro_Dv, M_accel_gryo_Dv
    optional imu Dvs for size effect: rx_i_Dv, ry_i_Dv, rz_i_Dv, Ix_Dv, Iy_Dv, Iz_Dv
camera Dvs: 
    T_c_b_Dv: aopt.TransformationDv (T_cNplus1_cN)
    cameraTimetoImuTimeDv: aopt.Scalar

#### Get values of design variables

gravityDv: toEuclidean()
accelBiasDv/gyroBiasDv: spline().eval(t) or evalD(t, 0)
q_i_b_Dv: toRotationMatrix()
r_b_Dv: toEuclidean()
poseSplineDv: sm.Transformation(T_w_b.toTransformationMatrix()) where T_w_b=transformationAtTime(timeExpression, 0.0, 0.0)

#### Error terms

CameraChainErrorTerms: error_t(frame, pidx, p) where error_t = self.camera.reprojectionErrorType + setMEstimatorPolicy

The variants of reprojectionErrorType derive from the SimpleReprojectionError C++ class which is exported to python in exportReprojectionError().

The different error types are grouped into a variety of camera models in terms of python classes defined in 
aslam_cv/aslam_cv_backend_python/python/aslam_cv_backend/\__init\__.py.

AccelerometerErrorTerms: ket.EuclideanError + setMEstimatorPolicy
GyroscopeErrorTerms: ket.EuclideanError + setMEstimatorPolicy
Accel and gyro BiasMotionTerms: BSplineEuclideanMotionError
PoseMotionTerms: MarginalizationPriorErrorTerm (by default inactive)


### kalibr_calibrate_rs_cameras
This calibration procedure currently (Dec 2021) supports only one camera.

#### Design variables

landmark_w_dv: aopt.HomogeneousPointDv (by default inactive)
__poseSpline_dv: asp.BSplinePoseDesignVariable
__camera_dv: The camera design variables are created by cameraModel.designVariable(self.geometry). 
The camera design variables are exported to python by exportCameraDesignVariables.
    projection: DesignVariableAdapter<projection_t>
    distortion: DesignVariableAdapter<distortion_t>
    shutter: DesignVariableAdapter<shutter_t>


#### Error terms

reprojectionErrorAdaptiveCovariance for RS models.

These error types derive from the CovarianceReprojectionError C++ class, which is exported to python by exportCovarianceReprojectionError.

The different error types are grouped into a variety of camera models in terms of python classes defined in aslam_cv/aslam_cv_backend_python/python/aslam_cv_backend/__init__.py.

Reprojection errors with adaptive covariance is developed solely for rolling shutter cameras as discussed in
Section 3.5 Error Term Standardisation, Oth et. al., Rolling shutter camera calibration.
Because of the error standardisation, these reprojection errors depend on not only the camera parameters, 
but also the pose B splines.

reprojectionError for GS models

These error types derive from the ReprojectionError C++ class which
is exported to python by exportReprojectionError. 
These errors depend on the camera parameters which may be optimized in the kalibr_calibrate_rs_cameras procedure.

regularizer: asp.BSplineMotionError of aslam_nonparametric_estimation/aslam_splines/include/aslam/backend.

## Citing
More information about RS camera-IMU calibration and simulation can be found
at the [report](https://arxiv.org/abs/2108.07200).
If you find the extension useful, please consider citing it.

```
@techreport{huaiContinuoustime2021,
  title = {Continuous-time spatiotemporal calibration of a rolling shutter camera-{{IMU}} system},
  author = {Huai, Jianzhu and Zhuang, Yuan and Yuan, Qicheng and Lin, Yukai},
  year = {2021},
  month = aug,
  eprint = {2108.07200},
  eprinttype = {arxiv},
  pages = {11},
  address = {{Wuhan, Hubei, China}},
  institution = {{Wuhan University}},
  url = {http://arxiv.org/abs/2108.07200},
  urldate = {2021-08-22},
  archiveprefix = {arXiv}
}
```

## Compiled rolling shutter parameters for consumer products
We have gathered rolling shutter parameters for consumer cameras from a variety of sources, and have calibrated such parameters for consumer cameras by a LED panel.
These parameters are provided [here](doc/rolling-shutter-table.md) for reference.


## TODOs
* Support calibrating cameras of FOV >= 180 degrees with the equidistant model.
This can be tested with the calibration sequences of the [TUM VI benchmark](https://vision.in.tum.de/data/datasets/visual-inertial-dataset)
which uses a camera with 192 degree FOV.

* Support covariance recovery.
Theoretically, it is possible to recover the covariances for estimated parameters by using Schur complement techniques.
However, the computation often takes too long so the covariance recovery functions in kalibr_calibrate_imu_camera and
kalibr_calibrate_rs_cameras.
The cause may be the strong coupling between adjacent control points in B-splines.