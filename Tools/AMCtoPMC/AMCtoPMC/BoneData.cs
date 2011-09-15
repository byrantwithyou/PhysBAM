using System;
using System.Collections.Generic;
using System.Text;

namespace AMCtoPMC
{
    class BoneData
    {
        // Bones
        public List<Bone> bonelist = new List<Bone>();

        // Hash to bones
        public Dictionary<string, Bone> bonemap = new Dictionary<string, Bone>();

        // Parsing parts
        public BoneData(string[] words)
        {
            // Parse other parts
            List<string[]> parts = ParseHelper.SplitWordStream(words, 1, words.Length-1, "begin");

            // Interpret words
            foreach (string[] part in parts)
            {
                Bone nb = new Bone(part);
                bonelist.Add(nb);
                bonemap.Add(nb.name, nb);
            }
        }

    }
}
