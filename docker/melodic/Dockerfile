# Install ROS 
FROM osrf/ros:melodic-desktop-bionic

#(Ubuntu 16.04: xenial, Ubuntu 14.04: trusty, Ubuntu 18.04: bionic)
ENV UBUNTU_VERSION=bionic
#(Ubuntu 16.04: kinetic, Ubuntu 14.04: indigo, Ubuntu 18.04: melodic)
ENV ROS_VERSION=melodic

ENV DEBIAN_FRONTEND=noninteractive

# Install kalib
# RUN apt-get update -y && apt-get install -y gnupg2
# RUN apt-key adv --recv-keys --keyserver keyserver.ubuntu.com F42ED6FBAB17C654
RUN apt-get update && apt-get install -y \ 
        python-setuptools \
        python-rosinstall \
        ipython \ 
        libeigen3-dev \
        libboost-all-dev \
        doxygen \
        libopencv-dev \
        ros-melodic-vision-opencv \
        ros-melodic-image-transport-plugins \
        ros-melodic-cmake-modules \
        # python-software-properties \
        software-properties-common \
        libpoco-dev \
        python-matplotlib \
        python-scipy \
        python-git \
        python-pip \
        python-igraph \ 
        ipython \
        libtbb-dev \
        libblas-dev \
        liblapack-dev \
        python-catkin-tools \
        libv4l-dev \
        wget \
        autoconf automake \
	libsuitesparse-dev  \
	python-tk

ENV KALIBR_WORKSPACE /kalibr_workspace

RUN 	mkdir -p $KALIBR_WORKSPACE/src &&\
	cd $KALIBR_WORKSPACE &&\
	catkin init &&\
	catkin config --extend /opt/ros/melodic &&\
	catkin config --cmake-args -DCMAKE_BUILD_TYPE=Release

RUN 	cd $KALIBR_WORKSPACE/src &&\
	git clone https://github.com/JzHuai0108/kalibr

RUN	cd $KALIBR_WORKSPACE &&\
	catkin build -DCMAKE_BUILD_TYPE=Release -j4

# https://answers.ros.org/question/356401/is-ros-already-sourced-if-i-do-it-from-a-dockerfile/
COPY ros_entrypoint.sh /ros_entrypoint.sh