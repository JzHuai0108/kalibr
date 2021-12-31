# Test the RS camera-IMU calibration with the dynamic sample provided by kalibr.
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/kalibr_ws /path/to/kalibr/sample/dynamic"
    echo "To begin the test, add the line_delay_nanoseconds field for the two cameras in dynamic/camchain.yaml"
    echo " and comment them out like below."
    echo "#  line_delay_nanoseconds: 0"
    echo "..."
    echo "#  line_delay_nanoseconds: 0"
    exit 1
fi

kalibr_ws=$1
datadir=$2
outputdir=$datadir

cd $kalibr_ws
catkin build -DCMAKE_BUILD_TYPE=Release -j4
source devel/setup.bash

test_camera_IMU_calibration() {
cd $outputdir

echo "Calibrate imu-camera system with GS model without specifying line_delay..."
mkdir -p gs_original
cd gs_original
sed -i "/line_delay_nanoseconds/c\#  line_delay_nanoseconds: 0" $datadir/camchain.yaml
kalibr_calibrate_imu_camera --target $datadir/april_6x6.yaml --cam $datadir/camchain.yaml --imu $datadir/imu_adis16448.yaml \
  --bag $datadir/dynamic.bag --bag-from-to 5 45 --dont-show-report --export-poses --saveVimap
cd ..

echo "Calibrate imu-camera system with GS model with specified line_delay 0..."
sed -i "/line_delay_nanoseconds/c\  line_delay_nanoseconds: 0" $datadir/camchain.yaml
mkdir -p gs
cd gs
kalibr_calibrate_imu_camera --target $datadir/april_6x6.yaml --cam $datadir/camchain.yaml --imu $datadir/imu_adis16448.yaml \
  --bag $datadir/dynamic.bag --bag-from-to 5 45 --dont-show-report
cd ..

echo "Calibrate imu-camera system with RS model..."
sed -i "/line_delay_nanoseconds/c\  line_delay_nanoseconds: 5000" $datadir/camchain.yaml
mkdir -p rs
cd rs
kalibr_calibrate_imu_camera --target $datadir/april_6x6.yaml --cam $datadir/camchain.yaml --imu $datadir/imu_adis16448.yaml \
  --bag $datadir/dynamic.bag --bag-from-to 5 45 --dont-show-report
cd ..

echo "Calibrate imu-camera system with RS model and estimate-line-delay enabled..."
sed -i "/line_delay_nanoseconds/c\  line_delay_nanoseconds: 5000" $datadir/camchain.yaml
mkdir -p rs_est
cd rs_est
kalibr_calibrate_imu_camera --target $datadir/april_6x6.yaml --cam $datadir/camchain.yaml --imu $datadir/imu_adis16448.yaml \
 --bag $datadir/dynamic.bag --bag-from-to 5 45 --estimate-line-delay --dont-show-report
cd ..

sed -i "/line_delay_nanoseconds/c\#  line_delay_nanoseconds: 0" $datadir/camchain.yaml
}

test_simulation() {
camimusim_template=$kalibr_ws/src/kalibr/aslam_offline_calibration/kalibr/config_templates/camchain_imucam_template.yaml
imu_template=$kalibr_ws/src/kalibr/aslam_offline_calibration/kalibr/config_templates/imu_template.yaml
targetyaml=$datadir/april_6x6.yaml

cd $outputdir

echo "Generate the pose and bias splines..."
mkdir -p gs_original
cd gs_original
kalibr_calibrate_imu_camera --target $targetyaml --cam $datadir/camchain.yaml --imu $datadir/imu_adis16448.yaml \
  --bag $datadir/dynamic.bag --bag-from-to 5 45 --dont-show-report --save-splines
cd ..

echo "Simulate camera and IMU with pose and bias splines..."
mkdir -p sim
cd sim
delay=41250
cp $camimusim_template camchain_imucam.yaml
sed -i "/line_delay_nanoseconds/c\  line_delay_nanoseconds: $delay" camchain_imucam.yaml
kalibr_simulate_imu_camera $outputdir/gs_original/bspline_pose.txt --cam camchain_imucam.yaml \
    --imu $imu_template --target $targetyaml --output_dir $outputdir/sim --dont-show-report
cd ..

echo "Calibrate imu-camera system with GS model on simulated data..."
mkdir -p gs_sim
cd gs_sim
cp $outputdir/sim/camchain.yaml camchain.yaml
sed -i "/line_delay_nanoseconds/c\  line_delay_nanoseconds: 0" camchain.yaml
kalibr_calibrate_imu_camera --cam camchain.yaml --imu $imu_template \
    --target $targetyaml --bag $outputdir/sim --bag-from-to 0 1000 \
    --timeoffset-padding 0.05 --timeoffset-pattern 0.08 --dont-show-report
cd ..

echo "Calibrate imu-camera system with RS model on simulated data..."
mkdir -p rs_sim
cd rs_sim
cp $outputdir/sim/camchain.yaml camchain.yaml
sed -i "/line_delay_nanoseconds/c\  line_delay_nanoseconds: 66000" camchain.yaml
kalibr_calibrate_imu_camera --cam camchain.yaml --imu $imu_template \
    --target $targetyaml --bag $outputdir/sim --bag-from-to 0 1000 --estimate-line-delay \
    --timeoffset-padding 0.05 --timeoffset-pattern 0.08 --dont-show-report # --recover-covariance
cd ..
}

test_noise_identification() {
echo "Calibrate imu-camera system with noise identification..."
constbias=1
inputdir=$datadir/dynamic.bag
targetyaml=$datadir/april_6x6.yaml

biasarg=""
if [ "$constbias" -ne "0" ]; then
   biasarg=" --constant_bias "
fi

echo "Imu noise identification on $inputdir, outputdir: $outputdir/noiseid, const bias ? $constbias."
mkdir -p $outputdir/noiseid
cd $outputdir/noiseid
cp $datadir/camchain.yaml camchain.yaml
sed -i "/line_delay_nanoseconds/c\  line_delay_nanoseconds: 3000" camchain.yaml
cmd="kalibr_imu_noise_identification --target $targetyaml --cam camchain.yaml \
      --imu $datadir/imu_adis16448.yaml --imu-models="calibrated" --bag $inputdir --bag-from-to 5 48 \
      --estimate-line-delay --timeoffset-padding 0.15 --timeoffset-pattern 0.08 \
      $biasarg --max_noise_opt_iters=3 --dont-show-report"

echo $cmd
$cmd
}

test_camera_IMU_calibration

test_simulation

test_noise_identification
