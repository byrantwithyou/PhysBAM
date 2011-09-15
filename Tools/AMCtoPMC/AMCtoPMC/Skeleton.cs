using System;
using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;
using System.IO;
//using Microsoft.DirectX;

namespace AMCtoPMC
{
    class Skeleton
    {
        // Output is a set of rigid bodies and articulation points
        public Root root;
        public BoneData bonedata;
        public Hierarchy hierarchy;

        // Computed rigid bodies and articulations
        public Dictionary<string,RigidBody> bodies = new Dictionary<string,RigidBody>();
        public Dictionary<string,Joint> joints = new Dictionary<string,Joint>();

        // Load skeleton from AMC file
        public Skeleton(string fname)
        {
            // Open file
            FileStream fs = new FileStream(fname, FileMode.Open, FileAccess.Read);
            StreamReader sr = new StreamReader(fs);

            // Read contents
            string contents = sr.ReadToEnd();
            sr.Close();
            
            // Make into string of words
            List<string> words = ParseHelper.CreateWordStream(contents);

            // Split into chunks
            List<string[]> chunks = ParseHelper.SplitWordStream(words.ToArray(), ":.*");
            
            // Now do a stream parsing
            foreach (string[] chunk in chunks)
            {
                // Look for key segments defined using :
                switch (chunk[0])
                {
                    case ":version":
                        break;
                    case ":name":
                        break;
                    case ":units":
                        break;
                    case ":documentation":
                        break;
                    case ":root":
                        root = new Root(chunk);
                        break;
                    case ":bonedata":
                        bonedata = new BoneData(chunk);
                        break;
                    case ":hierarchy":
                        //Note - unlike other regions, this is sensitive to newlines - thus, we
                        // need to reparse this chunk (unfortunately!)
                        Regex linematch = new Regex(@":hierarchy(.*)", RegexOptions.Singleline);
                        Match m = linematch.Match(contents);
                        string[] lines = ParseHelper.CreateLineStream(m.Groups[1].Value);
                        hierarchy = new Hierarchy(lines);
                        break;
                    default:
                        Console.WriteLine("Unknown part {0}, ignoring.", chunk[0]);
                        break;
                }
            }

            

            
            // Create root
            RigidBody rootbody = new RigidBody();            
            rootbody.name = "root";
            rootbody.bone = null;
            bodies.Add(rootbody.name, rootbody);

            // Build bodies from bones
            foreach (Bone bone in bonedata.bonelist)
            {
                // Ignore some bones for speed
                if (bone.name == "lfingers" || bone.name == "rfingers" ||
                    bone.name == "lthumb" || bone.name == "rthumb" ||
                    bone.name == "lhand" || bone.name == "rhand") continue;

                // Name segment after this
                RigidBody r = new RigidBody();                
                r.name = bone.name;
                r.bone = bone;                         

                bodies.Add(r.name,r);                
            }

            // Traverse hierarchy to build parent/child relationships
            foreach (Tree t in hierarchy.treelist)
            {
                foreach (string[] link in t.links)
                {
                    string parent = link[0];                    
                    for (int i = 1; i < link.Length; i++)
                    {
                        if (!bodies.ContainsKey(parent) || !bodies.ContainsKey(link[i])) continue;
                        bodies[parent].children.Add(bodies[link[i]]);
                        bodies[link[i]].parent = bodies[parent];
                    }
                }
            }

            // We can only compute global bone offsets, etc. by traversing hierarchy
            Stack<RigidBody> stack = new Stack<RigidBody>();                       
            stack.Push(rootbody);
            while (stack.Count > 0)
            {
                RigidBody b = stack.Pop();
                b.InitializeTransform();
                foreach (RigidBody c in b.children) stack.Push(c);
            }

            // Perform DFS collapse
            stack.Push(rootbody);
            while (stack.Count > 0)
            {
                RigidBody b = stack.Pop();

                // Non root elements w/o dof are collapsed into parent
                if (b.bone != null && b.bone.dof == null && b.bone.name != "root")
                {
                    // add out children to the parent
                    foreach (RigidBody child in b.children)
                    {
                        b.parent.children.Add(child);
                        child.parent = b.parent;
                    }

                    // remove us from the parent list
                    b.parent.children.Remove(b);
                }
                
                // Continue with dfs                
                foreach (RigidBody c in b.children) stack.Push(c);
            }

            
            // Set ids on segments, etc.
            int id = 1;
            stack.Push(rootbody);
            while (stack.Count > 0)
            {
                RigidBody b = stack.Pop();
                b.id = id++;                 
                foreach (RigidBody rb in b.children) stack.Push(rb);
            }
            
            // Articulation points are placed in world space at the segment start/ends, and transform them
            // into both parent and child space to set for the joint.
            // Also, all root elements shoudl be rigid
            int jid = 1;
            stack.Push(rootbody);
            while (stack.Count > 0)
            {
                RigidBody rb = stack.Pop();
                if (rb.parent != null)
                {
                    Joint j = new Joint();
                    j.world_position = rb.world_bone_offset;
                    j.bone = rb.bone;
                    j.parent = rb.parent;
                    j.child = rb;
                    j.dof = rb.bone.dof;
                    j.id = jid++;
                    j.name = rb.bone.name;

                    // Now, also print out joint's rotational axis
                    j.axis = rb.bone.axis;

                    // Important - takes on name of child, not parent
                    // This is due to the bone -> rb + joint mapping
                    joints.Add(rb.name, j);
                }
                foreach (RigidBody c in rb.children) stack.Push(c);
            }         
            
        }

        // Writes the skeleton out to a more PhysBAM-friendly format
        public void Export(string str)
        {
            FileStream fs = new FileStream(str, FileMode.Create, FileAccess.Write);
            StreamWriter sw = new StreamWriter(fs);

            // Convert from deg to rad            
            float conv = (float)(Math.PI / 180);

            // Count bodies in tree
            Stack<RigidBody> stack = new Stack<RigidBody>();
            stack.Push(bodies["root"]);
            int count = 0;
            while (stack.Count > 0)
            {
                RigidBody rb = stack.Pop();
                count++;
                foreach (RigidBody c in rb.children) stack.Push(c);
            }

            // Write out bodies and information
            sw.WriteLine(count);
            
            // Now write out bodies
            stack.Push(bodies["root"]);
            while (stack.Count > 0)
            {
                RigidBody rb = stack.Pop();
                sw.WriteLine(rb.name);
                sw.WriteLine(rb.world_center_offset.X);
                sw.WriteLine(rb.world_center_offset.Y);
                sw.WriteLine(rb.world_center_offset.Z);
                sw.WriteLine(rb.bone != null ? rb.bone.axis.X * conv: 0.0f);
                sw.WriteLine(rb.bone != null ? rb.bone.axis.Y * conv: 0.0f);
                sw.WriteLine(rb.bone != null ? rb.bone.axis.Z * conv: 0.0f);
                sw.WriteLine(rb.bone != null ? rb.bone.length : 1.0f);
                sw.WriteLine(rb.bone != null ? rb.bone.direction.X: 0.0f);
                sw.WriteLine(rb.bone != null ? rb.bone.direction.Y: 0.0f);
                sw.WriteLine(rb.bone != null ? rb.bone.direction.Z: 0.0f);
                foreach (RigidBody c in rb.children) stack.Push(c);
            }

            // Write out articulation information - degrees of freedom, etc.
            sw.WriteLine(joints.Count);
            foreach (Joint j in joints.Values)
            {
                sw.WriteLine(j.id);
                sw.WriteLine(j.name);
                sw.WriteLine(j.world_position.X);
                sw.WriteLine(j.world_position.Y);
                sw.WriteLine(j.world_position.Z);
                sw.WriteLine(j.axis.X * conv);
                sw.WriteLine(j.axis.Y * conv);
                sw.WriteLine(j.axis.Z * conv);
                sw.WriteLine(j.parent.id);
                sw.WriteLine(j.child.id);

                // write dof for this joint so we can generate constraints
                if (j.dof == null)
                {
                    sw.WriteLine(0);
                }
                else
                {
                    sw.WriteLine(j.dof.Length);
                    for (int k = 0; k < j.dof.Length; k++) sw.WriteLine(j.dof[k]);
                }
                
            }

            // Close now
            sw.Close();
        }

    }
}
