import csv
import sys

import numpy as np
import sew


def loadRawImuData(csvFile):
    timestamps = []  # list of acv.Time
    gyroData = []  # list of numpy array (3,)
    accelData = []  # list of numpy array (3,)
    with open(csvFile, 'r') as stream:
        csv_reader = csv.reader(stream , delimiter=',')

        for index, row in enumerate(csv_reader):
            if index == 0:
                continue
            timeString = row[0].strip()
            time = int(timeString[:-9]) + int(timeString[-9:]) * 1e-9
            timestamps.append(time)
            gyroData.append(np.array([float(row[4]), float(row[5]), float(row[6])]))
            accelData.append(np.array([float(row[1]), float(row[2]), float(row[3])]))
    return timestamps, gyroData, accelData


def main():
    if len(sys.argv) < 2:
        print('Usage:{} imu.csv'.format(sys.argv[0]))
        sys.exit(1)

    imucsv = sys.argv[1]
    timestamps, gyroData, accelData = loadRawImuData(imucsv)
    imu_gyro = np.transpose(gyroData)
    imu_acc = np.transpose(accelData)
    imu_t = timestamps

    gyroNoise, accNoise, so3_dt, r3_dt, dt = sew.identifyImuNoiseAndKnotSpacing(imu_t, imu_gyro, imu_acc)

    gyro_weight = np.sqrt(dt) / gyroNoise
    acc_weight = np.sqrt(dt) / accNoise

    print(
        'SO(3): knot spacing {:.1f} ms ({:.1f} Hz), '
        'gyroscope weight {:.2f}'.format(1000 * so3_dt, 1 / so3_dt, gyro_weight))
    print('   R3: knot spacing {:.1f} ms ({:.1f} Hz), '
          'accelerometer weight {:.2f}'.format(1000 * r3_dt, 1 / r3_dt,
                                               acc_weight))
    print('Predicted gyro noise {} acc noise {}'.format(gyroNoise, accNoise))


if __name__ == "__main__":
    main()
