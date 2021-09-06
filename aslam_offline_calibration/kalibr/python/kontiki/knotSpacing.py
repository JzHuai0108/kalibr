
import numpy as np
from sew import knot_spacing_and_variance


def identifyImuNoiseAndKnotSpacing(imu_t, imu_gyro, imu_acc):
    """

    :param imu_t: list or 1D np array, time in seconds of double type
    :param imu_gyro: 3xN array
    :param imu_acc: 3xN array
    :return: predicted gyro noise, predicted acc noise,
        optimal knot placing per gyro data,
        optimal knot placing per accel data.
    """

    min_dt = 0.01  # Don't go lower than this
    q_gyro = 0.99
    q_acc = 0.99
    verbose=False

    sew_opts = dict(min_dt=min_dt, verbose=verbose)
    def removeBias(data):
        mean = data.mean(axis=1)
        data = data - mean[:, np.newaxis]
        return data
    # The direct current / mean value does not affect predicted noises and knot spacing.
    # imu_gyro = removeBias(imu_gyro)
    # imu_acc = removeBias(imu_acc)
    # print("mean {}".format(imu_gyro.mean(axis=1)))

    so3_dt, so3_var = knot_spacing_and_variance(imu_gyro, imu_t, q_gyro, **sew_opts)
    r3_dt, r3_var = knot_spacing_and_variance(imu_acc, imu_t, q_acc, **sew_opts)

    gyroNoiseDiscrete = np.sqrt(so3_var)
    accNoiseDiscrete = np.sqrt(r3_var)
    dt = np.mean(np.diff(imu_t))

    gyroNoise = gyroNoiseDiscrete * np.sqrt(dt)
    accNoise = accNoiseDiscrete * np.sqrt(dt)

    return gyroNoise, accNoise, so3_dt, r3_dt, dt
