using System;
using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;


namespace AMCtoPMC
{
    class ParseHelper
    {
        public static Vector3 ParseVector3(string x, string y, string z)
        {
            return new Vector3(float.Parse(x), float.Parse(y), float.Parse(z));
        }

        public static List<string> CreateWordStream(string contents)
        {
            // Clean each line
            Regex c = new Regex(@"\r\n");
            string[] cleanlines = c.Split(contents);
            for (int i = 0; i < cleanlines.Length; i++)
            {
                // Trim leading spaces
                cleanlines[i] = cleanlines[i].TrimStart(' ');
                
                // Remove any comments on line
                Regex comment = new Regex(@"(.*)#.*");
                Match m = comment.Match(cleanlines[i]);
                if (m.Success)
                {
                    cleanlines[i] = m.Groups[1].Value;
                }
            }

            // Make into string of words
            List<string> words = new List<string>();
            for (int i = 0; i < cleanlines.Length; i++)
            {
                // Add spaces around parens, etc. so that symbols
                // are properly split during tokenization 
                cleanlines[i].Replace("(", " ( ");
                cleanlines[i].Replace(")", " ) ");

                // Split by spaces now
                string[] ws = cleanlines[i].Split(' ');
                foreach (string word in ws)
                {
                    if (word != "") words.Add(word);
                }
            }

            return words;
        }

        public static List<String[]> SplitWordStream(string[] words, params string[] keywords)
        {
            return SplitWordStream(words, 0, words.Length-1, keywords);
        }

        public static List<String[]> SplitWordStream(string[] words, int start, int end, params string[] keywords)
        {
            Regex r = new Regex(string.Join("|", keywords));
            List<string[]> chunks = new List<String[]>();

            int cstart = -1;
            for (int i = start; i <= end+1; i++)
            {
                // If we matched, add that chunk on
                if (i == end+1 || r.IsMatch(words[i]))
                {
                    if (cstart != -1)
                    {
                        string[] cword = new string[i - cstart];
                        for (int j = cstart; j < i; j++) cword[j - cstart] = words[j];
                        chunks.Add(cword);
                    }
                    cstart = i;
                }
            }

            return chunks;
        }


        public static string[] CreateLineStream(string contents)
        {
            // Clean each line
            Regex c = new Regex(@"\r\n");
            string[] cleanlines = c.Split(contents);
            for (int i = 0; i < cleanlines.Length; i++)
            {
                // Trim leading spaces
                cleanlines[i] = cleanlines[i].TrimStart(' ');

                // Remove any comments on line
                Regex comment = new Regex(@"(.*)#.*");
                Match m = comment.Match(cleanlines[i]);
                if (m.Success)
                {
                    cleanlines[i] = m.Groups[1].Value;
                }
            }

            return cleanlines;
        }
    }

}
