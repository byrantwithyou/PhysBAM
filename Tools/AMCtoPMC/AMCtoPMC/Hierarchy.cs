using System;
using System.Collections.Generic;
using System.Text;

namespace AMCtoPMC
{
    class Hierarchy
    {
        // Trees that we parsed
        public List<Tree> treelist = new List<Tree>();

        // Constructor - note it is line sensitive!
        public Hierarchy(string[] lines)
        {        
            // Interpret words
            for (int i = 0; i < lines.Length; i++)
            {
                if (lines[i] == "") continue;


                if (lines[i] == "begin")
                {
                    i++;
                    Tree t = new Tree();
                    for (; i < lines.Length && lines[i] != "end"; i++) t.links.Add(lines[i].Split(' '));
                    treelist.Add(t);
                }

            }
        }
    }
}
