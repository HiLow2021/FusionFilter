import math
import FusionFilter

# <summary>
# MadgwickAHRS class. Implementation of Madgwick's IMU and AHRS algorithms.
# 
# <remarks>
# See: http:#www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
# 
class MadgwickAHRS(FusionFilter.FusionFilter):

    # <summary>
    # Gets or sets the algorithm gain beta.
    # 
    beta = None
    
    # <summary>
    # Initializes a new instance of the class.
    # 
    # <param name="samplePeriod">
    # Sample period.
    # 
    # <param name="beta">
    # Algorithm gain beta.
    # 
    def __init__(self, samplePeriod,  beta = 0.1):
        super(MadgwickAHRS, self).__init__(samplePeriod)
        self.beta = beta
    
    
    # <summary>
    # Algorithm AHRS update method. Requires only gyroscope and accelerometer data.
    # 
    # <param name="gx">
    # Gyroscope x axis measurement in radians/s.
    # 
    # <param name="gy">
    # Gyroscope y axis measurement in radians/s.
    # 
    # <param name="gz">
    # Gyroscope z axis measurement in radians/s.
    # 
    # <param name="ax">
    # Accelerometer x axis measurement in any calibrated units.
    # 
    # <param name="ay">
    # Accelerometer y axis measurement in any calibrated units.
    # 
    # <param name="az">
    # Accelerometer z axis measurement in any calibrated units.
    # 
    # <param name="mx">
    # Magnetometer x axis measurement in any calibrated units.
    # 
    # <param name="my">
    # Magnetometer y axis measurement in any calibrated units.
    # 
    # <param name="mz">
    # Magnetometer z axis measurement in any calibrated units.
    # 
    # <remarks>
    # Optimised for minimal arithmetic.
    # Total +-: 160
    # Total *: 172
    # Total /: 5
    # Total sqrt: 5
    #  
    def update(self, gx, gy, gz, ax, ay, az, mx, my, mz):
    
        q0, q1, q2, q3 = self.quaternion   # short name local variable for readability        
        recipNorm = None
        s0, s1, s2, s3 = None, None, None, None
        qDot1, qDot2, qDot3, qDot4 = None, None, None, None
        hx, hy = None, None
        _2q0mx, _2q0my, _2q0mz, _2q1mx = None, None, None, None 
        _2bx, _2bz, _4bx, _4bz = None, None, None, None
        _2q0, _2q1, _2q2, _2q3 = None, None, None, None
        _2q0q2, _2q2q3 = None, None
        q0q0, q0q1, q0q2, q0q3 = None, None, None, None
        q1q1, q1q2, q1q3, q2q2, q2q3, q3q3 = None, None, None, None, None, None

        # Use IMU algorithm if magnetometer measurement invalid (as NaN in magnetometer normalisation)
        if (mx == 0.0) and (my == 0.0) and (mz == 0.0):
        
            self.update_without_mag(gx, gy, gz, ax, ay, az)
            return
        

        # Rate of change of quaternion from gyroscope
        qDot1 = 0.5 * (-q1 * gx - q2 * gy - q3 * gz)
        qDot2 = 0.5 * (q0 * gx + q2 * gz - q3 * gy)
        qDot3 = 0.5 * (q0 * gy - q1 * gz + q3 * gx)
        qDot4 = 0.5 * (q0 * gz + q1 * gy - q2 * gx)

        # Compute feedback only if accelerometer measurement valid (as NaN in accelerometer normalisation)
        if not ((ax == 0.0) and (ay == 0.0) and (az == 0.0)):
        
            # Normalise accelerometer measurement
            recipNorm = 1.0 / math.sqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            # Normalise magnetometer measurement
            recipNorm = 1.0 / math.sqrt(mx * mx + my * my + mz * mz)
            mx *= recipNorm
            my *= recipNorm
            mz *= recipNorm

            # Auxiliary variables to a repeated arithmetic
            _2q0mx = 2.0 * q0 * mx
            _2q0my = 2.0 * q0 * my
            _2q0mz = 2.0 * q0 * mz
            _2q1mx = 2.0 * q1 * mx
            _2q0 = 2.0 * q0
            _2q1 = 2.0 * q1
            _2q2 = 2.0 * q2
            _2q3 = 2.0 * q3
            _2q0q2 = 2.0 * q0 * q2
            _2q2q3 = 2.0 * q2 * q3
            q0q0 = q0 * q0
            q0q1 = q0 * q1
            q0q2 = q0 * q2
            q0q3 = q0 * q3
            q1q1 = q1 * q1
            q1q2 = q1 * q2
            q1q3 = q1 * q3
            q2q2 = q2 * q2
            q2q3 = q2 * q3
            q3q3 = q3 * q3

            # Reference direction of Earth's magnetic field
            hx = mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 + _2q1 * my * q2 + _2q1 * mz * q3 - mx * q2q2 - mx * q3q3
            hy = _2q0mx * q3 + my * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my * q1q1 + my * q2q2 + _2q2 * mz * q3 - my * q3q3
            _2bx = math.sqrt(hx * hx + hy * hy)
            _2bz = -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 + _2q1mx * q3 - mz * q1q1 + _2q2 * my * q3 - mz * q2q2 + mz * q3q3
            _4bx = 2.0 * _2bx
            _4bz = 2.0 * _2bz

            # Gradient decent algorithm corrective step
            s0 = -_2q2 * (2.0 * q1q3 - _2q0q2 - ax) + _2q1 * (2.0 * q0q1 + _2q2q3 - ay) - _2bz * q2 * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            s1 = _2q3 * (2.0 * q1q3 - _2q0q2 - ax) + _2q0 * (2.0 * q0q1 + _2q2q3 - ay) - 4.0 * q1 * (1 - 2.0 * q1q1 - 2.0 * q2q2 - az) + _2bz * q3 * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            s2 = -_2q0 * (2.0 * q1q3 - _2q0q2 - ax) + _2q3 * (2.0 * q0q1 + _2q2q3 - ay) - 4.0 * q2 * (1 - 2.0 * q1q1 - 2.0 * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            s3 = _2q1 * (2.0 * q1q3 - _2q0q2 - ax) + _2q2 * (2.0 * q0q1 + _2q2q3 - ay) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            recipNorm = 1.0 / math.sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3) # normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            # Apply feedback step
            qDot1 -= self.beta * s0
            qDot2 -= self.beta * s1
            qDot3 -= self.beta * s2
            qDot4 -= self.beta * s3
        

        # Integrate rate of change of quaternion to yield quaternion
        q0 += qDot1 * self.samplePeriod
        q1 += qDot2 * self.samplePeriod
        q2 += qDot3 * self.samplePeriod
        q3 += qDot4 * self.samplePeriod

        # Normalise quaternion
        recipNorm = 1.0 / math.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
        self.quaternion[0] = q0 * recipNorm
        self.quaternion[1] = q1 * recipNorm
        self.quaternion[2] = q2 * recipNorm
        self.quaternion[3] = q3 * recipNorm
        self._update_angles()
    
    
    # <summary>
    # Algorithm IMU update method. Requires only gyroscope and accelerometer data.
    # 
    # <param name="gx">
    # Gyroscope x axis measurement in radians/s.
    # 
    # <param name="gy">
    # Gyroscope y axis measurement in radians/s.
    # 
    # <param name="gz">
    # Gyroscope z axis measurement in radians/s.
    # 
    # <param name="ax">
    # Accelerometer x axis measurement in any calibrated units.
    # 
    # <param name="ay">
    # Accelerometer y axis measurement in any calibrated units.
    # 
    # <param name="az">
    # Accelerometer z axis measurement in any calibrated units.
    # 
    # <remarks>
    # Optimised for minimal arithmetic.
    # Total +-: 45
    # Total *: 85
    # Total /: 3
    # Total sqrt: 3
    # 
    def update_without_mag(self, gx, gy, gz, ax, ay, az):
    
        q0, q1, q2, q3 = self.quaternion   # short name local variable for readability 
        recipNorm = None
        s0, s1, s2, s3 = None, None, None, None
        qDot1, qDot2, qDot3, qDot4 = None, None, None, None
        _2q0, _2q1, _2q2, _2q3 = None, None, None, None 
        _4q0, _4q1, _4q2, _8q1, _8q2 = None, None, None, None, None 
        q0q0, q1q1, q2q2, q3q3 = None, None, None, None

        # Rate of change of quaternion from gyroscope
        qDot1 = 0.5 * (-q1 * gx - q2 * gy - q3 * gz)
        qDot2 = 0.5 * (q0 * gx + q2 * gz - q3 * gy)
        qDot3 = 0.5 * (q0 * gy - q1 * gz + q3 * gx)
        qDot4 = 0.5 * (q0 * gz + q1 * gy - q2 * gx)

        # Compute feedback only if accelerometer measurement valid (as NaN in accelerometer normalisation)
        if not ((ax == 0.0) and (ay == 0.0) and (az == 0.0)):
        
            # Normalise accelerometer measurement
            recipNorm = 1.0 / math.sqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            # Auxiliary variables to a repeated arithmetic
            _2q0 = 2.0 * q0
            _2q1 = 2.0 * q1
            _2q2 = 2.0 * q2
            _2q3 = 2.0 * q3
            _4q0 = 4.0 * q0
            _4q1 = 4.0 * q1
            _4q2 = 4.0 * q2
            _8q1 = 8.0 * q1
            _8q2 = 8.0 * q2
            q0q0 = q0 * q0
            q1q1 = q1 * q1
            q2q2 = q2 * q2
            q3q3 = q3 * q3

            # Gradient decent algorithm corrective step
            s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay
            s1 = _4q1 * q3q3 - _2q3 * ax + 4.0 * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az
            s2 = 4.0 * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az
            s3 = 4.0 * q1q1 * q3 - _2q1 * ax + 4.0 * q2q2 * q3 - _2q2 * ay
            recipNorm = 1.0 / math.sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3) # normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            # Apply feedback step
            qDot1 -= self.beta * s0
            qDot2 -= self.beta * s1
            qDot3 -= self.beta * s2
            qDot4 -= self.beta * s3
        

        # Integrate rate of change of quaternion to yield quaternion
        q0 += qDot1 * self.samplePeriod
        q1 += qDot2 * self.samplePeriod
        q2 += qDot3 * self.samplePeriod
        q3 += qDot4 * self.samplePeriod

        # Normalise quaternion
        recipNorm = 1.0 / math.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
        self.quaternion[0] = q0 * recipNorm
        self.quaternion[1] = q1 * recipNorm
        self.quaternion[2] = q2 * recipNorm
        self.quaternion[3] = q3 * recipNorm
        self._update_angles()