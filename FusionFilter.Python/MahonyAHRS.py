import math
import FusionFilter

# <summary>
# MahonyAHRS class. Madgwick's implementation of Mayhony's AHRS algorithm.
# 
# <remarks>
# See: http:#www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
# 
class MahonyAHRS(FusionFilter.FusionFilter):

    # <summary>
    # Gets or sets the algorithm proportional gain.
    # 
    kp = None
    
    # <summary>
    # Gets or sets the algorithm integral gain.
    # 
    ki = None    
    
    # <summary>
    # Gets or sets the integral error.
    # 
    integralFB = None
    
    # <summary>
    # Initializes a new instance of the <see cref="MadgwickAHRS"/> class.
    # 
    # <param name="samplePeriod">
    # Sample period.
    # 
    # <param name="kp">
    # Algorithm proportional gain.
    # 
    # <param name="ki">
    # Algorithm integral gain.
    # 
    def __init__(self, samplePeriod, kp = 1.0, ki = 0.0):
        super(MahonyAHRS, self).__init__(samplePeriod)
        self.kp = kp
        self.ki = ki
        self.integralFB = [0.0, 0.0, 0.0] 
    
    
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
    #  
    def update(self, gx, gy, gz, ax, ay, az, mx, my, mz):
    
        q0, q1, q2, q3 = self.quaternion   # short name local variable for readability
        recipNorm = None
        q0q0, q0q1, q0q2, q0q3 = None, None, None, None 
        q1q1, q1q2, q1q3, q2q2, q2q3, q3q3 = None, None, None, None, None, None
        hx, hy, bx, bz = None, None, None, None
        halfvx, halfvy, halfvz, halfwx, halfwy, halfwz = None, None, None, None, None, None
        halfex, halfey, halfez = None, None, None
        qa, qb, qc = None, None, None

        # Use IMU algorithm if magnetometer measurement invalid (as NaN in magnetometer normalisation)
        if (mx == 0.0) and (my == 0.0) and (mz == 0.0):
        
            self.update_without_mag(gx, gy, gz, ax, ay, az)
            return
        

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
            hx = 2.0 * (mx * (0.5 - q2q2 - q3q3) + my * (q1q2 - q0q3) + mz * (q1q3 + q0q2))
            hy = 2.0 * (mx * (q1q2 + q0q3) + my * (0.5 - q1q1 - q3q3) + mz * (q2q3 - q0q1))
            bx = math.sqrt(hx * hx + hy * hy)
            bz = 2.0 * (mx * (q1q3 - q0q2) + my * (q2q3 + q0q1) + mz * (0.5 - q1q1 - q2q2))

            # Estimated direction of gravity and magnetic field
            halfvx = q1q3 - q0q2
            halfvy = q0q1 + q2q3
            halfvz = q0q0 - 0.5 + q3q3
            halfwx = bx * (0.5 - q2q2 - q3q3) + bz * (q1q3 - q0q2)
            halfwy = bx * (q1q2 - q0q3) + bz * (q0q1 + q2q3)
            halfwz = bx * (q0q2 + q1q3) + bz * (0.5 - q1q1 - q2q2)

            # Error is sum of cross product between estimated direction and measured direction of field vectors
            halfex = (ay * halfvz - az * halfvy) + (my * halfwz - mz * halfwy)
            halfey = (az * halfvx - ax * halfvz) + (mz * halfwx - mx * halfwz)
            halfez = (ax * halfvy - ay * halfvx) + (mx * halfwy - my * halfwx)

            # Compute and apply integral feedback if enabled
            if self.ki > 0.0:
            
                self.integralFB[0] += self.ki * halfex * self.samplePeriod    # integral error scaled by Ki
                self.integralFB[1] += self.ki * halfey * self.samplePeriod
                self.integralFB[2] += self.ki * halfez * self.samplePeriod
                gx += self.integralFB[0]  # apply integral feedback
                gy += self.integralFB[1]
                gz += self.integralFB[2]
            
            else:
            
                self.integralFB[0] = 0.0 # prevent integral windup
                self.integralFB[1] = 0.0
                self.integralFB[2] = 0.0
            

            # Apply proportional feedback
            gx += self.kp * halfex
            gy += self.kp * halfey
            gz += self.kp * halfez
        

        # Integrate rate of change of quaternion
        gx *= (0.5 * self.samplePeriod)     # pre-multiply common factors
        gy *= (0.5 * self.samplePeriod)
        gz *= (0.5 * self.samplePeriod)
        qa = q0
        qb = q1
        qc = q2
        q0 += (-qb * gx - qc * gy - q3 * gz)
        q1 += (qa * gx + qc * gz - q3 * gy)
        q2 += (qa * gy - qb * gz + q3 * gx)
        q3 += (qa * gz + qb * gy - qc * gx)

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
    def update_without_mag(self, gx, gy, gz, ax, ay, az):
    
        q0, q1, q2, q3 = self.quaternion   # short name local variable for readability
        recipNorm = None
        halfvx, halfvy, halfvz = None, None, None
        halfex, halfey, halfez = None, None, None
        qa, qb, qc = None, None, None

        # Compute feedback only if accelerometer measurement valid (as NaN in accelerometer normalisation)
        if not ((ax == 0.0) and (ay == 0.0) and (az == 0.0)):
        
            # Normalise accelerometer measurement
            recipNorm = 1.0 / math.sqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            # Estimated direction of gravity and vector perpendicular to magnetic flux
            halfvx = q1 * q3 - q0 * q2
            halfvy = q0 * q1 + q2 * q3
            halfvz = q0 * q0 - 0.5 + q3 * q3

            # Error is sum of cross product between estimated and measured direction of gravity
            halfex = (ay * halfvz - az * halfvy)
            halfey = (az * halfvx - ax * halfvz)
            halfez = (ax * halfvy - ay * halfvx)

            # Compute and apply integral feedback if enabled
            if self.ki > 0.0:
            
                self.integralFB[0] += self.ki * halfex * self.samplePeriod    # integral error scaled by Ki
                self.integralFB[1] += self.ki * halfey * self.samplePeriod
                self.integralFB[2] += self.ki * halfez * self.samplePeriod
                gx += self.integralFB[0]  # apply integral feedback
                gy += self.integralFB[1]
                gz += self.integralFB[2]
            
            else:
            
                self.integralFB[0] = 0.0 # prevent integral windup
                self.integralFB[1] = 0.0
                self.integralFB[2] = 0.0
            

            # Apply proportional feedback
            gx += self.kp * halfex
            gy += self.kp * halfey
            gz += self.kp * halfez
        

        # Integrate rate of change of quaternion
        gx *= (0.5 * self.samplePeriod)     # pre-multiply common factors
        gy *= (0.5 * self.samplePeriod)
        gz *= (0.5 * self.samplePeriod)
        qa = q0
        qb = q1
        qc = q2
        q0 += (-qb * gx - qc * gy - q3 * gz)
        q1 += (qa * gx + qc * gz - q3 * gy)
        q2 += (qa * gy - qb * gz + q3 * gx)
        q3 += (qa * gz + qb * gy - qc * gx)

        # Normalise quaternion
        recipNorm = 1.0 / math.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
        self.quaternion[0] = q0 * recipNorm
        self.quaternion[1] = q1 * recipNorm
        self.quaternion[2] = q2 * recipNorm
        self.quaternion[3] = q3 * recipNorm
        self._update_angles()