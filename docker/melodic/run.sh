#!/bin/bash

# run
data_dir=$1
if [ ! -d "$data_dir" ]; then
  echo "data directory does not exist: $data_dir"
  exit 1
fi

xhost +
if [ "$2" = "true" ];then
  docker run -it -e DISPLAY -e QT_X11_NO_MITSHM=1 -v /tmp/.X11-unix:/tmp/.X11-unix:rw -v $data_dir:/root/data melodic-kalib-docker /bin/bash -c "cd /root/data; /bin/bash"
else
  docker run -it -v /tmp/.X11-unix:/tmp/.X11-unix:rw -v $data_dir:/root/data melodic-kalib-docker /bin/bash -c "cd /root/data; /bin/bash"
fi


