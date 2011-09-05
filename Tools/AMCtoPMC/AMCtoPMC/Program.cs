using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace AMCtoPMC
{
    class Program
    {
        static void Main(string[] args)
        {
            // Read skeleton ASF file
            Skeleton s = new Skeleton(args[0]);

            // Read AMC file and parse
            Motion m = new Motion(s, args[1]);

            // Transform out skeleton, use to compute PhysBAM joint angles and COM ?
            s.Export(args[0]+".conv");

            // Transform amc motion to set of joint angles
            m.Export(args[1] + ".conv");
        }
    }
}
