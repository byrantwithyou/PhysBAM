using System;
using System.Collections.Generic;
using System.Text;

namespace AMCtoPMC
{
    struct Vector3
    {
        public float X;
        public float Y;
        public float Z;

        public Vector3(float x, float y, float z)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
        }

        public static Vector3 operator*(Vector3 v, float f)
        {
            return new Vector3(v.X * f, v.Y * f, v.Z * f);
        }

        public static Vector3 operator+(Vector3 v, Vector3 v2)
        {
            return new Vector3(v.X + v2.X, v.Y + v2.Y, v.Z + v2.Z);
        }
    }
}
