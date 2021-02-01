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

2. Based on the B-spline model, simulate rolling shutter camera measurements

```
kalibr_simulate_imu_camera $output_dir/bspline_pose.txt --cam $script_dir/camchain_template.yaml --imu $script_dir/imu_template.yaml \
  --target $data_dir/april_6x6.yaml --output_dir $output_dir
```

# A crash course on design variables
Design variable groups used in kalibr\_calibrate\_imu_camera
```
poseDv(spline)
gravityDv(active)
imu Dvs: gyroBiasDv(spline), accelBiasDv(spline), q_i_b_Dv(inactive), r_b_Dv(inactive), 
potential imu Dvs for scaledMisalignment: q_gyro_i_Dv, M_accel_Dv, M_gyro_Dv, M_accel_gryo_Dv
potential imu Dvs for size effect: rx_i_Dv, ry_i_Dv, rz_i_Dv, Ix_Dv, Iy_Dv, Iz_Dv
camera Dvs: T_c_b_Dv, cameraTimetoImuTimeDv
```

```
design variable     gravityDv      accelBiasDv/gyroBiasDv  q_i_b_Dv           r_b_Dv        poseSplineDv
type                trivial        spline()                trivial            trivial       T_w_b=transformationAtTime(timeExpression, 0.0, 0.0)
access its content  toEuclidean()  spline().evalD(t,0)     toRotationMatrix() toEuclidean() sm.Transformation(T_w_b.toTransformationMatrix())
```

