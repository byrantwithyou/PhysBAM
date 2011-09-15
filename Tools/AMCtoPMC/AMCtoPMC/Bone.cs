using System;
using System.Collections.Generic;
using System.Text;


namespace AMCtoPMC
{
    class Bone
    {
        // Attributes
        public int id;
        public string name;
        public Vector3 direction;
        public float length;
        public string[] dof;

        public Vector3 axis;
        public string axisorder;

        public string[] limits;

        // Parsing parts
        public Bone()
        {
        }

        // Parsing parts
        public Bone(string[] words)
        {
            // Parse other parts
            List<string[]> parts = ParseHelper.SplitWordStream(words, 1, words.Length-1, "id", "name", "direction", "length", "dof", "axis", "limits");

            // Interpret words
            foreach (string[] part in parts)
            {
                switch (part[0])
                {
                    case "id":
                        this.id = int.Parse(part[1]);
                        break;
                    case "name":
                        this.name = part[1];
                        break;
                    case "direction":
                        this.direction = ParseHelper.ParseVector3(part[1], part[2], part[3]);                        
                        break;
                    case "length":
                        this.length = float.Parse(part[1]);
                        break;
                    case "axis":
                        this.axis = ParseHelper.ParseVector3(part[1], part[2], part[3]);                        
                        this.axisorder = part[4];
                        break;
                    case "dof":
                        this.dof = new string[part.Length - 1];
                        for (int i = 1; i < part.Length; i++) this.dof[i - 1] = part[i];
                        break;
                    case "limits":
                        this.limits = part;
                        break;
                }
            }
        }
    }
}
