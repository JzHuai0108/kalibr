import numpy as np
import kalibr_common as kc
from . import Simulator


def testAddNoise():

    imuParameters = kc.ImuParameters(None, True)
    imuParameters.data["initial_accelerometer_bias"] = np.array([0.2, 0.2, 0.2])
    imuParameters.data["initial_gyro_bias"] = np.array([0.1, 0.1, 0.1])

    imuParameters.data["accelerometer_noise_density"] = 5e-2
    imuParameters.data["accelerometer_random_walk"] = 1e-3
    imuParameters.data["gyroscope_noise_density"] = 1e-2
    imuParameters.data["gyroscope_random_walk"] = 1e-4

    # imuParameters.data["accelerometer_noise_density"] = 1e-2
    # imuParameters.data["accelerometer_random_walk"] = 2e-4
    # imuParameters.data["gyroscope_noise_density"] = 5e-3
    # imuParameters.data["gyroscope_random_walk"] = 4e-6
    rate = 200
    imuParameters.data["update_rate"] = rate

    imuMeasurements = np.random.rand(rate * 50, 6)

    noisyImuMeasurements, trueBiases = Simulator.addNoiseToImuReadings(imuMeasurements, imuParameters)

    noises = noisyImuMeasurements - trueBiases - imuMeasurements

    accGyroNoiseDiscrete = np.std(noises, 0, ddof=1)

    accGyroWalkDiscrete =np.std(np.diff(trueBiases, axis=0), 0, ddof=1)
    rootf = np.sqrt(rate)
    rootdt = 1.0 / rootf

    _, expectedAccWalk, expectedAccNoise = imuParameters.getAccelerometerStatistics()
    _, expectedGyroWalk, expectedGyroNoise = imuParameters.getGyroStatistics()

    accGyroNoise = accGyroNoiseDiscrete * rootdt
    actualAccNoise = np.mean(accGyroNoise[:3])
    actualGyroNoise = np.mean(accGyroNoise[3:])

    accGyroWalk = accGyroWalkDiscrete * rootf
    actualAccWalk = np.mean(accGyroWalk[:3])
    actualGyroWalk = np.mean(accGyroWalk[3:])

    assert abs(actualAccWalk - expectedAccWalk) < expectedAccWalk * 0.01
    assert abs(actualGyroWalk - expectedGyroWalk) < expectedGyroWalk * 0.01

    assert abs(actualAccNoise - expectedAccNoise) < expectedAccNoise * 0.01
    assert abs(actualGyroNoise - expectedGyroNoise) < expectedGyroNoise * 0.01

    # print('Expected noise density {} {} random walk {} {}'.format(
    #     expectedAccNoise, expectedGyroNoise, expectedAccWalk, expectedGyroWalk))
    print('Actual noise density {} {} random walk {} {}'.format(
        actualAccNoise, actualGyroNoise, actualAccWalk, actualGyroWalk))

    biasedNoises = noisyImuMeasurements - imuMeasurements
    accGyroNoiseDiscrete = np.std(biasedNoises, 0, ddof=1)
    accGyroNoise = accGyroNoiseDiscrete * rootdt

    print('Actual noise densities without discounting biases {} {}'.format(
        np.mean(accGyroNoise[:3]), np.mean(accGyroNoise[3:])))
