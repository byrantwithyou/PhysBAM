#!/usr/bin/python
#!/usr/bin/python

import physbam
import math
def subset(data,items):
    return map(lambda x: data[x],items)

class EMBEDDING_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.tests=physbam.SOLIDS_STANDARD_TESTS_Vf3(self)


        self.last_frame=1000
        self.frame_rate=24
        self.solids_parameters.cfl=.9
        self.solids_parameters.perform_self_collision=True
        self.output_directory="straight"
        self.printed_header=False
        self.g=9.8
        #self.Unite_Fragments()
        self.data_directory="/data/aselle/hair_modeling"

    def Initialize_Bodies(self):
        # helper references
        deformable_object=self.solids_parameters.deformable_object
        particles=deformable_object.particles
        rigid_particles=deformable_object.rigid_body_particles
        binding_list=deformable_object.binding_list
        soft_bindings=deformable_object.soft_bindings
            

        ################################################################
        # Geometry Phase
        ################################################################
        # spiral
        def X_spiral(s):
            return physbam.Vf3(math.cos(50*s),20*s,math.sin(50*s))

        #curve=self.tests.Create_Segmented_Curve("hairs.curve.gz",physbam.RIGID_BODY_STATE_Vf3(),False,False)
        self.fixed_nodes=physbam.LA_i()
        self.fixed_values=physbam.LA_Vf3()
        for i in range(1,5):
            self.fixed_nodes.Append(i)

        #physbam.Read_From_File("float","fixed_nodes",self.fixed_nodes)
        self.edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        self.extra_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        self.bending_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        self.torsion_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        self.volume=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)


        physbam.Read_From_File("float","particles",particles)
        physbam.Read_From_File("float","edges.curve",self.edges.mesh)
        physbam.Read_From_File("float","extra_edges.curve",self.extra_edges.mesh)
        physbam.Read_From_File("float","bending_edges.curve",self.bending_edges.mesh)
        physbam.Read_From_File("float","torsion_edges.curve",self.torsion_edges.mesh)
        physbam.Read_From_File("float","torsion_edges.curve",self.torsion_edges.mesh)
        physbam.Read_From_File("float","tets.gz",self.volume.mesh)
        particles.Store_Velocity(True)

        self.count=len(particles)


        self.edges.Set_Mass_Of_Particles(True,True)
        self.extra_edges.Set_Mass_Of_Particles(True,True)

        for i in [self.edges,self.extra_edges,self.bending_edges,self.torsion_edges]:
            i.Update_Number_Nodes()
        # todo fix mass
        for i in range(1,len(particles)+1):
            particles.mass[i]=1
        
        
        
        #volume=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)
        i=1
        #for i in range(1,len(self.volume.mesh.elements)+1):
        for nodes in self.volume.mesh.elements:
            #print i
            #print physbam.TETRAHEDRON_f(*map(lambda i:particles.X[i],nodes)).Signed_Volume()
            if physbam.TETRAHEDRON_f(*map(lambda i:particles.X[i],nodes)).Signed_Volume()<0:
                nodes[1],nodes[2],nodes[3],nodes[4]=(nodes[1],nodes[2],nodes[4],nodes[3])
            i+=1
            #print physbam.TETRAHEDRON_f(*map(lambda i:particles.X[i],nodes)).Signed_Volume()
        for nodes in self.volume.mesh.elements:
            #print i
            #print physbam.TETRAHEDRON_f(*map(lambda i:particles.X[i],nodes)).Signed_Volume()
            if physbam.TETRAHEDRON_f(*map(lambda i:particles.X[i],nodes)).Signed_Volume()<0:
                print "FUCK"
            #print physbam.TETRAHEDRON_f(*map(lambda i:particles.X[i],nodes)).Signed_Volume()
        self.volume.Update_Number_Nodes()
        
        ################################################################
        # Mass Update Phase
        ################################################################
        binding_list.Distribute_Mass_To_Parents(particles.mass.array)
        binding_list.Clear_Hard_Bound_Particles(particles.mass.array)
        particles.mass.Compute_Auxiliary_Attributes(soft_bindings)
        deformable_object.soft_bindings.Set_Mass_From_Effective_Mass()

        ################################################################
        # Forces
        ################################################################
        # Create forces
        tet_node_list=physbam.LA_i()
        deformable_object.Add_Force(physbam.GRAVITY_Vf3(particles,rigid_particles,True,True,self.g,physbam.Vf3(0,-1,0)))
        overdamping_fraction=40
        springs=[
            physbam.Create_Edge_Springs(self.edges,1e4,overdamping_fraction,True,.1,True,0.,True,True),
            physbam.Create_Edge_Springs(self.extra_edges,1e4,overdamping_fraction,True,.1,True,0.,True,True),
            physbam.Create_Edge_Springs(self.bending_edges,1e4,overdamping_fraction,True,.1,True,0.,True,True),
            physbam.Create_Edge_Springs(self.torsion_edges,1e4,overdamping_fraction,True,.1,True,0.,True,True),
            physbam.Create_Tet_Springs(self.volume,1e3,overdamping_fraction,False,.1,True,.1,True,0.,True,True)
            ]
        
        for force in springs:
            force.Clamp_Restlength(.01)
            deformable_object.Add_Force(force)

        #deformable_object.Add_Force(physbam.Create_Tet_Springs(self.volume,1e4,2.,False,.1,True,.1,True,0.,True))
        #deformable_object.Add_Force(physbam.Create_Altitude_Springs(volume,1e4,20.,False,.1,True,.1,True,0.,True))

        deformable_object.Add_Structure(self.edges)
        #deformable_object.Add_Structure(self.extra_edges)
        #deformable_object.Add_Structure(self.bending_edges)
        #deformable_object.Add_Structure(self.torsion_edges)
        #deformable_object.Add_Structure(self.volume)
        #deformable_object.Add_Structure(torsion_curve)
        
        # update fragments
        deformable_object.Update_Fragments()

        #self.tests.Add_Ground(.1,-2,0.5)
        #self.tests.Add_Rigid_Body("efty_100k_nose_filled",1,1,True)
        #center of mass = -0.00479698 14.673 8.13414
        #orientation = 0.999839 -0.0162252 -0.00339366 0.00689811
    

        #self.Collide_Segmented_Curve(self.edges)
        self.Collide_All_Bodies()

    def Set_External_Velocities(self,V,velocity_time,current_position_time,super_fragment):
        #self.Zero_Nodes(self.fixed_nodes,V)

        #V[1]=0
        #if velocity_time<3.9: V[self.count]=physbam.Vf3(-.5,0,0)

        X=self.solids_parameters.deformable_object.particles.X
        self.Zero_Nodes(self.fixed_nodes,V)
        offset=physbam.Vf3(-.1,0,0)
        if velocity_time>2.: offset=physbam.Vf3(0,0,0)
        for i in range(self.count-4,self.count+1):
            V[i]=4*physbam.Vf3(1,0,0).Cross_Product(X[i]-physbam.Vf3(X[i][1],0,0)) + offset
        #for i in range(self.count-4,self.count+1):
        #    V[i]=physbam.Vf3(-1,0,0)
        
    def Zero_Out_Enslaved_Velocity_Nodes(self,V,velocity_time,current_position_time,super_fragment):
        #self.Zero_Nodes(self.fixed_nodes,V)
        
        #V[1]=0
        #if velocity_time<3.9: V[self.count]=0

        self.Zero_Nodes(self.fixed_nodes,V)
        for i in range(self.count-4,self.count+1):
            V[i]=0


    def Preprocess_Frame(self,frame): pass
    def Postprocess_Frame(self,frame):
        percent=float(frame)/float(self.last_frame)*100.
        percent_int=int(percent)/2
        print "\r%5d of %5d (%3.1f%%) %s"%(frame,self.last_frame,percent,("["+"*"*percent_int+" "*(50-percent_int)+"]"))
    def Update_Time_Varying_Material_Properties(self,time,super_fragment): pass
    def Update_Collision_Body_Positions_And_Velocities(self,time): pass
    def Preprocess_Solids_Substep(self,time,substep): pass
    def Postprocess_Solids_Substep(self,time,substep): pass
    def Update_Solids_Parameters(self,time): pass
    def Set_External_Positions(self,X,time,super_fragment): pass
    def Apply_Constraints(self,dt,time): pass
    def Add_External_Forces(self,F,time,super_fragment): pass
    def Add_External_Forces(self,wrench,time,super_fragment): pass
    def Limit_Solids_Dt(self,time): pass
    def Update_Fragments(self): pass
    def Postprocess_Fragments(self): pass
    def Postprocess_Super_Fragments(self,swap_pairs,rebuild,old_max): pass
    def Update_Super_Fragments(self,swap_pairs,rebuild,old_max): pass
    def Add_External_Super_Fragment_Connectivity(self,union_find): pass
    def Add_External_Fragment_Connectivity(self,union_find,particle_is_simulated): pass

physbam.LOG.Initialize_Logging(True,True,0,False)
physbam.Set_Floating_Point_Exception(True,True,True,True,False,False)
#physbam.LOG.Initialize_Logging(False,False,9999,False)
example=EMBEDDING_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
#del example
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

