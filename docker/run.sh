#!/bin/bash
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <folder-of-rosbag> [ros_distro] [display]"
  echo "ros_distro can be melodic (default) or noetic"
  echo "display can be 0 or 1 (default)"
  exit 1
fi

data_dir=$1
if [ ! -d "$data_dir" ]; then
  echo "Data directory does not exist: $data_dir."
  exit 2
fi

dros_distro=melodic
if [ "$2" = "noetic" ]; then
  dros_distro=noetic
fi

# Explanations of arguments to the below commands are given in 
# https://stackoverflow.com/questions/43015536/xhost-command-for-docker-gui-apps-eclipse
if [ "$3" = "0" ]; then
  docker run -it -v $data_dir:/root/data kalibr:$dros_distro /bin/bash -c "cd /root/data; /bin/bash"
else
  xhost +local:root;
  docker run -it -e DISPLAY -e QT_X11_NO_MITSHM=1 -v /tmp/.X11-unix:/tmp/.X11-unix:rw -v $data_dir:/root/data kalibr:$dros_distro /bin/bash -c "cd /root/data; /bin/bash"
  xhost -local:root;
fi
