using System;
using System.Collections.Generic;
using System.Text;


namespace AMCtoPMC
{
    class RigidBody
    {
        // Id
        public int id;

        // Name of body
        public string name;

        // Mass, etc.
        public bool initialized = false;
        //public Matrix world_to_local;
        public Vector3 world_bone_offset;
        public Vector3 world_center_offset;

        // Parent and children
        public List<RigidBody> children = new List<RigidBody>();
        public RigidBody parent;

        // Bone that we represent
        public Bone bone;

        // Constructor
        public RigidBody()
        {

        }

        // Tell the body to compute matrix
        public void InitializeTransform()
        {
            // If we're initialized, we must be done!
            if (this.initialized) return;

            // Traverse up, to make sure everyone is initialized
            if (parent != null)
            {
                parent.InitializeTransform();
                if (parent.bone != null)
                {
                    this.world_bone_offset = parent.world_bone_offset + parent.bone.direction * parent.bone.length;
                }
                if (this.bone != null)
                {
                    this.world_center_offset = this.world_bone_offset + this.bone.direction * (0.5f * this.bone.length);
                }
            }
            else
            {
                this.world_bone_offset = new Vector3(0.0f, 0.0f, 0.0f);
                this.world_center_offset = new Vector3(0.0f, 0.0f, 0.0f);
            }
            
            // Set init
            this.initialized = true;
            
        }
    }
}
