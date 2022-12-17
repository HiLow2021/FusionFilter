using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FusionFilter
{
    /// <summary>
    /// MahonyAHRS class. Madgwick's implementation of Mayhony's AHRS algorithm.
    /// </summary>
    /// <remarks>
    /// See: http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
    /// </remarks>
    public class MahonyAHRS : FusionFilter
    {
        /// <summary>
        /// Gets or sets the algorithm proportional gain.
        /// </summary>
        public float Kp { get; set; }

        /// <summary>
        /// Gets or sets the algorithm integral gain.
        /// </summary>
        public float Ki { get; set; }

        /// <summary>
        /// Gets or sets the integral error.
        /// </summary>
        private float[] IntegralFB { get; set; }

        /// <summary>
        /// Initializes a new instance of the <see cref="MadgwickAHRS"/> class.
        /// </summary>
        /// <param name="samplePeriod">
        /// Sample period.
        /// </param>
        public MahonyAHRS(float samplePeriod) : this(samplePeriod, 1f, 0f)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="MadgwickAHRS"/> class.
        /// </summary>
        /// <param name="samplePeriod">
        /// Sample period.
        /// </param>
        /// <param name="kp">
        /// Algorithm proportional gain.
        /// </param> 
        public MahonyAHRS(float samplePeriod, float kp) : this(samplePeriod, kp, 0f)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="MadgwickAHRS"/> class.
        /// </summary>
        /// <param name="samplePeriod">
        /// Sample period.
        /// </param>
        /// <param name="kp">
        /// Algorithm proportional gain.
        /// </param>
        /// <param name="ki">
        /// Algorithm integral gain.
        /// </param>
        public MahonyAHRS(float samplePeriod, float kp, float ki) : base(samplePeriod)
        {
            Kp = kp;
            Ki = ki;
            IntegralFB = new float[] { 0f, 0f, 0f };
        }

        /// <summary>
        /// Algorithm AHRS update method. Requires only gyroscope and accelerometer data.
        /// </summary>
        /// <param name="gx">
        /// Gyroscope x axis measurement in radians/s.
        /// </param>
        /// <param name="gy">
        /// Gyroscope y axis measurement in radians/s.
        /// </param>
        /// <param name="gz">
        /// Gyroscope z axis measurement in radians/s.
        /// </param>
        /// <param name="ax">
        /// Accelerometer x axis measurement in any calibrated units.
        /// </param>
        /// <param name="ay">
        /// Accelerometer y axis measurement in any calibrated units.
        /// </param>
        /// <param name="az">
        /// Accelerometer z axis measurement in any calibrated units.
        /// </param>
        /// <param name="mx">
        /// Magnetometer x axis measurement in any calibrated units.
        /// </param>
        /// <param name="my">
        /// Magnetometer y axis measurement in any calibrated units.
        /// </param>
        /// <param name="mz">
        /// Magnetometer z axis measurement in any calibrated units.
        /// </param>
        /// <remarks>
        /// Optimised for minimal arithmetic.
        /// </remarks> 
        public override void Update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz)
        {
            float q0 = Quaternion[0], q1 = Quaternion[1], q2 = Quaternion[2], q3 = Quaternion[3];   // short name local variable for readability
            float recipNorm;
            float q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;
            float hx, hy, bx, bz;
            float halfvx, halfvy, halfvz, halfwx, halfwy, halfwz;
            float halfex, halfey, halfez;
            float qa, qb, qc;

            // Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
            if ((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f))
            {
                Update(gx, gy, gz, ax, ay, az);
                return;
            }

            // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
            if (!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f)))
            {

                // Normalise accelerometer measurement
                recipNorm = 1.0f / (float)Math.Sqrt(ax * ax + ay * ay + az * az);
                ax *= recipNorm;
                ay *= recipNorm;
                az *= recipNorm;

                // Normalise magnetometer measurement
                recipNorm = 1.0f / (float)Math.Sqrt(mx * mx + my * my + mz * mz);
                mx *= recipNorm;
                my *= recipNorm;
                mz *= recipNorm;

                // Auxiliary variables to avoid repeated arithmetic
                q0q0 = q0 * q0;
                q0q1 = q0 * q1;
                q0q2 = q0 * q2;
                q0q3 = q0 * q3;
                q1q1 = q1 * q1;
                q1q2 = q1 * q2;
                q1q3 = q1 * q3;
                q2q2 = q2 * q2;
                q2q3 = q2 * q3;
                q3q3 = q3 * q3;

                // Reference direction of Earth's magnetic field
                hx = 2.0f * (mx * (0.5f - q2q2 - q3q3) + my * (q1q2 - q0q3) + mz * (q1q3 + q0q2));
                hy = 2.0f * (mx * (q1q2 + q0q3) + my * (0.5f - q1q1 - q3q3) + mz * (q2q3 - q0q1));
                bx = (float)Math.Sqrt(hx * hx + hy * hy);
                bz = 2.0f * (mx * (q1q3 - q0q2) + my * (q2q3 + q0q1) + mz * (0.5f - q1q1 - q2q2));

                // Estimated direction of gravity and magnetic field
                halfvx = q1q3 - q0q2;
                halfvy = q0q1 + q2q3;
                halfvz = q0q0 - 0.5f + q3q3;
                halfwx = bx * (0.5f - q2q2 - q3q3) + bz * (q1q3 - q0q2);
                halfwy = bx * (q1q2 - q0q3) + bz * (q0q1 + q2q3);
                halfwz = bx * (q0q2 + q1q3) + bz * (0.5f - q1q1 - q2q2);

                // Error is sum of cross product between estimated direction and measured direction of field vectors
                halfex = (ay * halfvz - az * halfvy) + (my * halfwz - mz * halfwy);
                halfey = (az * halfvx - ax * halfvz) + (mz * halfwx - mx * halfwz);
                halfez = (ax * halfvy - ay * halfvx) + (mx * halfwy - my * halfwx);

                // Compute and apply integral feedback if enabled
                if (Ki > 0.0f)
                {
                    IntegralFB[0] += Ki * halfex * SamplePeriod;    // integral error scaled by Ki
                    IntegralFB[1] += Ki * halfey * SamplePeriod;
                    IntegralFB[2] += Ki * halfez * SamplePeriod;
                    gx += IntegralFB[0];  // apply integral feedback
                    gy += IntegralFB[1];
                    gz += IntegralFB[2];
                }
                else
                {
                    IntegralFB[0] = 0.0f; // prevent integral windup
                    IntegralFB[1] = 0.0f;
                    IntegralFB[2] = 0.0f;
                }

                // Apply proportional feedback
                gx += Kp * halfex;
                gy += Kp * halfey;
                gz += Kp * halfez;
            }

            // Integrate rate of change of quaternion
            gx *= (0.5f * SamplePeriod);     // pre-multiply common factors
            gy *= (0.5f * SamplePeriod);
            gz *= (0.5f * SamplePeriod);
            qa = q0;
            qb = q1;
            qc = q2;
            q0 += (-qb * gx - qc * gy - q3 * gz);
            q1 += (qa * gx + qc * gz - q3 * gy);
            q2 += (qa * gy - qb * gz + q3 * gx);
            q3 += (qa * gz + qb * gy - qc * gx);

            // Normalise quaternion
            recipNorm = 1.0f / (float)Math.Sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
            Quaternion[0] = q0 * recipNorm;
            Quaternion[1] = q1 * recipNorm;
            Quaternion[2] = q2 * recipNorm;
            Quaternion[3] = q3 * recipNorm;
            UpdateAngles();
        }

        /// <summary>
        /// Algorithm IMU update method. Requires only gyroscope and accelerometer data.
        /// </summary>
        /// <param name="gx">
        /// Gyroscope x axis measurement in radians/s.
        /// </param>
        /// <param name="gy">
        /// Gyroscope y axis measurement in radians/s.
        /// </param>
        /// <param name="gz">
        /// Gyroscope z axis measurement in radians/s.
        /// </param>
        /// <param name="ax">
        /// Accelerometer x axis measurement in any calibrated units.
        /// </param>
        /// <param name="ay">
        /// Accelerometer y axis measurement in any calibrated units.
        /// </param>
        /// <param name="az">
        /// Accelerometer z axis measurement in any calibrated units.
        /// </param>
        public override void Update(float gx, float gy, float gz, float ax, float ay, float az)
        {
            float q0 = Quaternion[0], q1 = Quaternion[1], q2 = Quaternion[2], q3 = Quaternion[3];   // short name local variable for readability
            float recipNorm;
            float halfvx, halfvy, halfvz;
            float halfex, halfey, halfez;
            float qa, qb, qc;

            // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
            if (!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f)))
            {

                // Normalise accelerometer measurement
                recipNorm = 1.0f / (float)Math.Sqrt(ax * ax + ay * ay + az * az);
                ax *= recipNorm;
                ay *= recipNorm;
                az *= recipNorm;

                // Estimated direction of gravity and vector perpendicular to magnetic flux
                halfvx = q1 * q3 - q0 * q2;
                halfvy = q0 * q1 + q2 * q3;
                halfvz = q0 * q0 - 0.5f + q3 * q3;

                // Error is sum of cross product between estimated and measured direction of gravity
                halfex = (ay * halfvz - az * halfvy);
                halfey = (az * halfvx - ax * halfvz);
                halfez = (ax * halfvy - ay * halfvx);

                // Compute and apply integral feedback if enabled
                if (Ki > 0.0f)
                {
                    IntegralFB[0] += Ki * halfex * SamplePeriod;    // integral error scaled by Ki
                    IntegralFB[1] += Ki * halfey * SamplePeriod;
                    IntegralFB[2] += Ki * halfez * SamplePeriod;
                    gx += IntegralFB[0];  // apply integral feedback
                    gy += IntegralFB[1];
                    gz += IntegralFB[2];
                }
                else
                {
                    IntegralFB[0] = 0.0f; // prevent integral windup
                    IntegralFB[1] = 0.0f;
                    IntegralFB[2] = 0.0f;
                }

                // Apply proportional feedback
                gx += Kp * halfex;
                gy += Kp * halfey;
                gz += Kp * halfez;
            }

            // Integrate rate of change of quaternion
            gx *= (0.5f * SamplePeriod);     // pre-multiply common factors
            gy *= (0.5f * SamplePeriod);
            gz *= (0.5f * SamplePeriod);
            qa = q0;
            qb = q1;
            qc = q2;
            q0 += (-qb * gx - qc * gy - q3 * gz);
            q1 += (qa * gx + qc * gz - q3 * gy);
            q2 += (qa * gy - qb * gz + q3 * gx);
            q3 += (qa * gz + qb * gy - qc * gx);

            // Normalise quaternion
            recipNorm = 1.0f / (float)Math.Sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
            Quaternion[0] = q0 * recipNorm;
            Quaternion[1] = q1 * recipNorm;
            Quaternion[2] = q2 * recipNorm;
            Quaternion[3] = q3 * recipNorm;
            UpdateAngles();
        }
    }
}