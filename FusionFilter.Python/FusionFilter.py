import math

class FusionFilter(object):

    radiansToDegrees = 57.29578
    degreesToRadians = 0.0174533

    # <summary>
    # Gets or sets the sample period.
    #
    samplePeriod = None

    # <summary>
    # Gets or sets the Quaternion output.
    #
    quaternion = None

    roll = None

    pitch = None

    yaw = None

    def __init__(self, samplePeriod):
        self.samplePeriod = samplePeriod
        self.quaternion = [1.0, 0.0, 0.0, 0.0]

    def get_fusion(self):
        roll, pitch, yaw = self.get_fusion_in_radians()
        roll *= FusionFilter.radiansToDegrees
        pitch *= FusionFilter.radiansToDegrees
        yaw = yaw * FusionFilter.radiansToDegrees + 180.0
        
        return roll, pitch, yaw

    def get_fusion_in_radians(self):
        return self.roll, self.pitch, self.yaw 

    def _update_angles(self):
        q0, q1, q2, q3 = self.quaternion
        self.roll = math.atan2(q0 * q1 + q2 * q3, 0.5 - q1 * q1 - q2 * q2)
        self.pitch = math.asin(-2.0 * (q1 * q3 - q0 * q2))
        self.yaw = math.atan2(q1 * q2 + q0 * q3, 0.5 - q2 * q2 - q3 * q3)