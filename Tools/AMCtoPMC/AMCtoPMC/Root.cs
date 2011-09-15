using System;
using System.Collections.Generic;
using System.Text;


namespace AMCtoPMC
{
    class Root
    {
        // Axis
        public string[] order;
        public string axis;
        public Vector3 position;
        public Vector3 orientation;

        // Parsing parts
        public Root(string[] words)
        {
            // Parse other parts
            List<string[]> parts = ParseHelper.SplitWordStream(words, 1, words.Length-1, "order", "axis", "position", "orientation");

            // Interpret words
            foreach (string[] part in parts)
            {
                switch (part[0])
                {
                    case "order":
                        order = new string[part.Length - 1];
                        for (int i = 1; i < part.Length; i++) order[i - 1] = part[i];
                        break;
                    case "axis":
                        axis = part[1];
                        break;
                    case "position":
                        position = ParseHelper.ParseVector3(part[1], part[2], part[3]);
                        break;
                    case "orientation":
                        orientation = ParseHelper.ParseVector3(part[1], part[2], part[3]);
                        break;
                }
            }
        }

    }

}
