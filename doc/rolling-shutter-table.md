## Rolling shutter parameters for consumer cameras

### Table 1 The readout times estimated by using a LED panel, a function generator and an oscilloscope, copied from [1].

Note that frame readout time = image height * line delay.

| Camera      | Resolution  | Frame rate | Mean frame readout time | Frame readout time standard deviation | #Estimates |
| - | - | - | - | - | - |
| name and identifier | px x px | Hz | millisecond | millisecond | 1 |
|iPhone 3GS | 640 x 480 | 30 | 30.84 | 0.36 | 7 |
|iPhone 4(#1) | 1280 × 720 | 30 | 31.98 | 0.25 | 6 |
|iPhone 4(#2) | 1280 × 720 | 30 | 32.37 | 0.27  | 4 |
|iPhone 4(#1,front) | 640 × 480 | 30 | 29.99 | 0.16  | 5 |
|iTouch 4 | 1280 × 720 | 30 | 30.07 | 0.25 | 5 |
|iTouch 4 (front) | 640 × 480 | 30 | 31.39 | 0.07 | 5 |
|HTC Desire | 640 × 480 | < 30 | 57.84∗ | 0.051 | 3 |
|HTC Desire | 1280 × 720 | < 30 | 43.834∗ | 0.034 | 4 |
|SE W890i | 320 × 240 | 14.706 | 60.78 | 0.16 | 4 |
|SE Xperia X10 | 640 × 480 | < 30 | 28.386∗ | 0.024 | 4 |
|SE Xperia X10 | 800 × 480 | < 30 | 27.117∗ | 0.040 | 5 |

\* reports one of several readout times mostly likely used by the camera in that setting.

### Table 2 The readout times estimated by panning the camera while viewing a vertical line, copied from [2].

For the full list, please refer to [2].
Note that rolling shutter amount = frame readout time * frame rate.

| Camera      | Resolution  | Frame rate | Rolling shutter amount | Rolling shutter amount std. dev. | Mean frame readout time* | Frame readout time std. dev. |
| - | - | - | - | - | - | - |
| name and identifier | px x px | Hz | 0.01 | 0.01 | millisecond | millisecond |
GoPro Hero3+ Black Edition, wide FOV | 1920x1080 | 60 |  88.89 | |14.82 | |
GoPro Hero3+ Black Edition, narrow FOV | 1920x1080 | 60 | 86.05| |14.34 | |
iPhone 4|1280x720|30|97|2|32.33|0.67
IPhone 5|1920x1080|30|82.1|0.1|27.37|0.03
IPhone 5|1920x1080|24|65.8|0.1|27.42|0.04
IPhone 5, stabilization off, 3rd party app|1920x1080|30|89.6|0.1|29.87|0.03
IPhone 5, 3rd party app|720x480|30|67.8|0.2|22.60|0.07
Samsung Galaxy S3|1920x1080|30|90|1|30.00|0.33

### Table 3 The readout times estimated by moving the camera by hand in front of a checkboard, copied from [3].

| Camera      | Resolution  | Frame rate | Mean line delay | Line delay std. dev. | Mean frame readout time* | Frame readout time std. dev. | #Estimates |
| - | - | - | - | - | - | - | - |
| name and identifier | px x px | Hz | 0.01 | 0.01 | millisecond | millisecond | 1 |
iPhone3GS|640x480|30|64.41|0.11|30.92|0.05|5
iPhone4S|1920x1080|30|24.12|0.52|26.05|0.56|12
Galaxy S3|1920x1080|30|30.25|0.68|32.67|0.73|7

### Table 4 The readout times for drone cameras, copied from [4].

|Drone|Camera|Resolution (px x px)|Shutter type|Sensor|Lens|Horizontal / vertical field of view (deg)|Estimated readout time (ms)
|-|-|-|-|-|-|-|-
|DJI Phantom 2 Vision+|FC200|4384x3288|rolling|1/2.3” CMOS|fisheye|110/80|74
DJI Inspire 1|FC300X|4000x3000|rolling|1/2.3” CMOS|perspective|85/70|30
3DR Solo|GoPro 4 Black|4000x3000|rolling|1/2.3” CMOS|fisheye|125/95|30
senseFly eBee|Canon S110|4000x3000|global|1/1.7” CMOS|perspective|71/56|0


### Table 5 The readout time estimated by a LED panel [5].
| Camera      | Resolution  | Frame rate | Mean line delay | Line delay std. dev. | Nominal frame readout time | Mean frame readout time* | Frame readout time std. dev. | #Estimates |
| - | - | - | - | - | - | - | - | - |
| name and identifier | px x px | Hz | microsecond | microsecond | millisecond | millisecond | millisecond | 1 |
Honor V10|1280x720| 30 | 28.72|0.18|20.733|20.68|0.13|121
iPhone 6S|1280x720| 30 | 20.85|0.49|x|15.01|0.35|122
Phab2 Pro|1280x720| 30 | 43.6|0.41|30.255|31.39|0.30|85
Zenfone AR|1280x720| 26 | 10.3|0.43|5.469|7.42|0.31|100


## References
[1] E. Ringaby and P. Forssén, "Efficient Video Rectification and Stabilisation for Cell-Phones," International Journal of Computer Vision, 96, 335-352, 2011. https://www.diva-portal.org/smash/get/diva2:505943/FULLTEXT01.pdf.

[2] G. Thalin, Deshaker video stabilizer plugin v2.5 for VirtualDub, 2010. http://www.guthspot.se/video/deshaker.htm

[3] L. Oth, P. Furgale, L. Kneip and R. Siegwart, "Rolling Shutter Camera Calibration," 2013 IEEE Conference on Computer Vision and Pattern Recognition, 2013, pp. 1360-1367, doi: 10.1109/CVPR.2013.179.

[4] J. Vautherin, S. Rutishauser, K. Schneider-Zapp, H. F. Choi, V. Chovancova, A. Glass, and C. Strecha, "Photogrammetric accuracy and modeling of rolling shutter cameras," ISPRS Annals of Photogrammetry, Remote Sensing & Spatial Information Sciences, 3(3), 2016.

[5] J. Huai, Y. Zhuang, B. Wang, C. Zhang, Y. Shao, J. Tang, and A. Yilmaz, "Automated rolling shutter calibration with an LED panel," Optics Letters 48, 847-850, 2023.
