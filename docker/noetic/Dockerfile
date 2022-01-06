# Install ROS 
FROM osrf/ros:noetic-desktop-full

#(Ubuntu 16.04: kinetic, Ubuntu 14.04: indigo, Ubuntu 18.04: melodic, Ubuntu 20.04: noetic)
ENV ROS_VERSION=noetic

ENV DEBIAN_FRONTEND=noninteractive

# Install kalib
RUN apt-get update && apt-get install -y \ 
        python3-setuptools \
        python3-rosinstall \
        ipython3 \
        libeigen3-dev \
        libboost-all-dev \
        doxygen \
        libopencv-dev \
        ros-noetic-vision-opencv \
        ros-noetic-image-transport-plugins \
        ros-noetic-cmake-modules \
        python3-software-properties \
        software-properties-common \
        libpoco-dev \
        python3-matplotlib \
        python3-scipy \
        python3-git \
        python3-pip \
        libtbb-dev \
        libblas-dev \
        liblapack-dev \
        libv4l-dev \
        python3-catkin-tools \
        python3-igraph \
        libsuitesparse-dev \
        python3-tk

RUN python3 -m pip install --upgrade pip; python3 -m pip install wxPython

ENV KALIBR_WORKSPACE /kalibr_workspace

RUN 	mkdir -p $KALIBR_WORKSPACE/src &&\
	cd $KALIBR_WORKSPACE &&\
	catkin init &&\
	catkin config --extend /opt/ros/$ROS_VERSION &&\
	catkin config --cmake-args -DCMAKE_BUILD_TYPE=Release

RUN 	cd $KALIBR_WORKSPACE/src &&\
	git clone https://github.com/JzHuai0108/kalibr

RUN	cd $KALIBR_WORKSPACE &&\
	catkin build -DCMAKE_BUILD_TYPE=Release -j4

COPY ros_entrypoint.sh /ros_entrypoint.sh
