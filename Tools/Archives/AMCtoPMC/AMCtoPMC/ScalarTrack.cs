using System;
using System.IO;
using System.Collections.Generic;
using System.Text;

namespace AMCtoPMC
{
    class Track
    {
        public List<Key> keys = new List<Key>();

        public Track()
        {

        }

        public void AddKey(Key k)
        {
            keys.Add(k);            
        }        

        public void Export(StreamWriter s)
        {
            keys.Sort();
            s.WriteLine(keys.Count);
            for (int i = 0; i < keys.Count; i++)
            {
                s.WriteLine("{0} {1} {2} {3}", keys[i].time, keys[i].value.X, keys[i].value.Y, keys[i].value.Z);                
            }
        }
    }

    class Key : IComparable
    {
        public float time;
        public Vector3 value;

        public Key(float time, Vector3 value)
        {
            this.time = time;
            this.value = value;
        }

        public int CompareTo(object other)
        {
            return time.CompareTo(((Key)other).time);
        }
    
    }
}
