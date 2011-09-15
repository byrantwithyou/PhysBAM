using System;
using System.Collections.Generic;
using System.Text;

namespace AMCtoPMC
{
    class Joint
    {
        // Type
        public string[] dof = new string[0];

        // Id
        public int id;
        public string name;

        // Bone
        public Bone bone;

        // World space point
        public Vector3 world_position;
        public Vector3 axis;

        // Parent and child bodies
        public RigidBody parent;
        public RigidBody child;

        // Motion track
        public Track angles_track = new Track();        
            
        // Constructor
        public Joint()
        {

        }
    }
}
