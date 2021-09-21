import sm
import numpy as np


def testRotationVector():
    """
    This test depends on several lines of import_rotational_kinematics_python() in 
    kalibr/Schweizer-Messer/sm_python/src/export_rotational_kinematics.cpp, 
    which were commented out to disable warnings.
    """
    rvi = sm.RotationVectorImpl()

    rotations = []
    qs = []
    for i in range(1000):
        q = np.random.rand(4)
        q = q / np.linalg.norm(q)
        C = sm.quat2r(q)
        qs.append(q)
        rotations.append(C)

    for i, C in enumerate(rotations):
        a = rvi.rotationMatrixToParametersOriginal(C)
        b = rvi.rotationMatrixToParametersClassic(C)
        d = sm.quat2AxisAngle(qs[i])
        if not np.allclose(a, -b, 1e-6):
            print("\nC: {}\n    Rotation vector: kalibr impl {} conventional impl {}".format(C, a, b))
        assert np.allclose(a, d)

    rotations = [np.array([[1.0000000, 0.0000000, 0.0000000],
                           [0.0000000, -1.0000000, 0.0000000],
                           [0.0000000, 0.0000000, -1.0000000]]),
                 np.array([[0.0000000, -1.0000000, 0.0000000],
                           [0.0000000, 0.0000000, -1.0000000],
                           [1.0000000, 0.0000000, 0.0000000]]),
                 np.array([[0.0000000, -1.0000000, 0.0000000],
                           [-1.0000000, 0.0000000, 0.0000000],
                           [0.0000000, 0.0000000, -1.0000000]]),
                 np.array([[-0.7827305, 0.6175095, 0.0775572],
                           [-0.2027229, -0.1351517, -0.9698647],
                           [-0.5884186, -0.7748652, 0.2309706]])]

    for i, C in enumerate(rotations):
        a = rvi.rotationMatrixToParametersOriginal(C)
        b = rvi.rotationMatrixToParametersClassic(C)
        if not np.allclose(a, -b):
            print('\nC: {}\n    Rotation vector: kalibr impl {} conventional impl {}'.format(C, a, b))
