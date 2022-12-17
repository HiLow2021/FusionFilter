using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FusionFilter
{
    public abstract class FusionFilter
    {
        public static readonly float radiansToDegrees = 57.29578f;
        public static readonly float degreesToRadians = 0.0174533f;

        /// <summary>
        /// Gets or sets the sample period.
        /// </summary>
        public virtual float SamplePeriod { get; set; }

        /// <summary>
        /// Gets or sets the Quaternion output.
        /// </summary>
        public virtual float[] Quaternion { get; set; }

        public virtual float Roll { get; private set; }
        public virtual float Pitch { get; private set; }
        public virtual float Yaw { get; private set; }

        public FusionFilter(float samplePeriod)
        {
            SamplePeriod = samplePeriod;
            Quaternion = new float[] { 1f, 0f, 0f, 0f };
        }

        public abstract void Update(float gx, float gy, float gz, float ax, float ay, float az);
        public abstract void Update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);

        public virtual void GetFusion(out float roll, out float pitch, out float yaw)
        {
            GetFusionInRadians(out roll, out pitch, out yaw);
            roll *= radiansToDegrees;
            pitch *= radiansToDegrees;
            yaw = yaw * radiansToDegrees + 180.0f;
        }

        public virtual void GetFusionInRadians(out float roll, out float pitch, out float yaw)
        {
            roll = Roll;
            pitch = Pitch;
            yaw = Yaw;
        }

        protected virtual void UpdateAngles()
        {
            float q0 = Quaternion[0], q1 = Quaternion[1], q2 = Quaternion[2], q3 = Quaternion[3];

            Roll = (float)Math.Atan2(q0 * q1 + q2 * q3, 0.5f - q1 * q1 - q2 * q2);
            Pitch = (float)Math.Asin(-2.0f * (q1 * q3 - q0 * q2));
            Yaw = (float)Math.Atan2(q1 * q2 + q0 * q3, 0.5f - q2 * q2 - q3 * q3);
        }
    }
}
