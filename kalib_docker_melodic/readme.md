BUILD docker:
1. Give the script execution permissions by running: 
    "chmod +x build.sh"
2. Run to build kalib-docker:
     "./build.sh"

RUN docker:
1. Give the script execution permissions by running: 
    "chmod +x run.sh".

2. Run using : 
    "./run.sh <path-to-data-dir>", 
    where <path-to-data-dir> is the path to the directory on your computer which contains your calibration bag file.

3. Step 1-2  will open an interactive shell session inside the docker container. The data directory will be mounted inside the container at the path:
    "/root/data". 
 The shell has the Kalibr workspace loaded. This means you can run your favorite Kalibr commands such as: 
    "kalibr_calibrate_cameras --bag bag.bag --topics /cam_node/left_raw /cam_node/right_raw --models pinhole-radtan pinhole-radtan --target target.yaml"