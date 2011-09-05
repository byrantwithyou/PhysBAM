using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Text.RegularExpressions;

namespace AMCtoPMC
{
    class Motion
    {
        // Skeleton we're built against
        Skeleton s;

        // Motion
        public Motion(Skeleton s, string fname)
        {
            // Store skeleton
            this.s = s;

            // Read in input data
            FileStream fs = new FileStream(fname, FileMode.Open, FileAccess.Read);
            StreamReader sr = new StreamReader(fs);

            // Convert from deg to rad            
            float conv = (float)(Math.PI / 180);

            int frame = 0;
            while (!sr.EndOfStream) 
            {
                string fullline = sr.ReadLine();
                if (fullline.TrimStart(' ')[0] == '#' || fullline.TrimStart(' ')[0] == ':') continue;

                // Parse
                string[] line = fullline.Split(' ');
                
                if (line.Length == 1)
                {
                    frame = int.Parse(line[0]);
                }
                else
                {
                    string name = line[0];

                    // Skip joints that aren't in our skeleton
                    if (!s.joints.ContainsKey(name)) continue;
                    Joint j = s.joints[name];

                    // Temporary storage
                    float x = 0, y = 0, z = 0;
                
                    // Parse motion specifiers based on DOF
                    for (int i = 0; i < j.dof.Length; i++) {
                        switch (j.dof[i])
                        {
                            case "rx":
                                x = conv * float.Parse(line[i + 1]);                                
                                break;
                            case "ry":
                                y = conv * float.Parse(line[i + 1]);                                
                                break;
                            case "rz":
                                z = conv * float.Parse(line[i + 1]);                                
                                break;
                        }
                    }

                    // Add euler angle key to this frame
                    j.angles_track.AddKey(new Key(frame, new Vector3(x, y, z)));
                }
            }

            sr.Close();
        }


        // Export to a file
        public void Export(string fname)
        {
            FileStream fs = new FileStream(fname, FileMode.Create, FileAccess.Write);
            StreamWriter sw = new StreamWriter(fs);

            // Write out bodies and information
            sw.WriteLine(s.joints.Count);
            foreach (Joint j in s.joints.Values)
            {
                sw.WriteLine(j.id);
                j.angles_track.Export(sw);             
            }

            // Close now
            sw.Close();
        }

    }
}
