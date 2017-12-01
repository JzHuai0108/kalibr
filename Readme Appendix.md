Install Kalibr from source on Ubuntu 16.04 with ROS Kinetic

# Install dependencies as instructed on the Wiki of Kalibr on github
```
sudo apt-get install python-setuptools python-rosinstall ipython libeigen3-dev libboost-all-dev doxygen libopencv-dev ros-kinetic-vision-opencv ros-kinetic-image-transport-plugins ros-kinetic-cmake-modules python-software-properties software-properties-common libpoco-dev python-matplotlib python-scipy python-git python-pip ipython libtbb-dev libblas-dev liblapack-dev python-catkin-tools libv4l-dev 
```
# Install the additional python packages
```
sudo apt-get install python-wxversion python-wxtools
```
# Install python-igraph via
```
sudo pip install python-igraph --upgrade
```
In installing python-igraph, there may arise an error reads as in the below reference
reference: https://stackoverflow.com/questions/37495375/python-pip-install-throws-typeerror-unsupported-operand-types-for-retry
Even using the solution suggested in the reference page, the error in installing python-igraph persists.

Don't fuss about it. Simply ignore it and go ahead to build kalibr with catkin_make. The building procedure succeeds in my case.

# Simulate monocular inertial data

## Download the imu camera calibration dataset, and generate a B spline model
```
kalibr_calibrate_imu_camera --cam camchain.yaml --target april_6x6.yaml --imu imu_adis16448.yaml --bag dynamic.bag
```

## Simulate rolling shutter camera measurements
```
kalibr_simulate_imu_camera "knotCoeffT.txt" april_6x6.yaml
```


