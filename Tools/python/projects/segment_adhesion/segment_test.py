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
        self.output_directory="adhesion"
        self.printed_header=False
        self.use_adhesion=True
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

        #curve=self.tests.Create_Segmented_Curve("hairs.curve.gz",physbam.RIGID_BODY_STATE_Vf3(),False,False)
        self.fixed_nodes=physbam.LA_i()
        self.fixed_values=physbam.LA_Vf3()

        #physbam.Read_From_File("float","fixed_nodes",self.fixed_nodes)
        self.edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)

        n=10
        particles.Add_Particles(n)
        base=1
        self.fixed_nodes.Append(base)
        for i in range(n):
            particles.X[base+i]=physbam.Vf3(.5+i*.01,0,.5);particles.mass[base+i]=1
            if i!=n-1: self.edges.mesh.elements.Append(base+i,base+i+1)
        self.fixed_nodes.Append(base+9)
        base+=n
        self.fixed_nodes.Append(base)
        particles.Add_Particles(n)
        for i in range(n):
            particles.X[base+i]=physbam.Vf3(.52,.02,.49+i*.01);particles.mass[base+i]=1
            if i!=n-1: self.edges.mesh.elements.Append(base+i,base+i+1)

            
        self.edges.Update_Number_Nodes()

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
        overdamping_fraction=100
        springs=[
            physbam.Create_Edge_Springs(self.edges,1e4,overdamping_fraction,True,.1,True,0.,True),
            physbam.Create_Segment_Bending_Springs(self.edges,1e3,overdamping_fraction,True,.1,True,0.,True),
            ]
        if self.use_adhesion:
            self.adhesion=physbam.SEGMENT_ADHESION_Vf3(particles,self.edges.mesh)
            self.adhesion.Set_Parameters(10.,10.,.003,.005) # ym,overdamp,start dist, stop dist
            deformable_object.Add_Force(self.adhesion)
            self.adhesion.Write_State(physbam.Get_Stream_Type("float"),self.output_directory+"/adhesion.%d"%0)
        
        for force in springs:
            try:
                force.Clamp_Restlength(.01)
            except:
                pass
            deformable_object.Add_Force(force)

        # update fragments
        deformable_object.Update_Fragments()

        #self.tests.Add_Ground(.1,-2,0.5)
        #self.tests.Add_Rigid_Body("efty_100k_nose_filled",1,1,True)
        #center of mass = -0.00479698 14.673 8.13414
        #orientation = 0.999839 -0.0162252 -0.00339366 0.00689811

        deformable_object.Add_Structure(self.edges)


    def Set_External_Velocities(self,V,velocity_time,current_position_time,super_fragment):
        self.Zero_Nodes(self.fixed_nodes,V)
        if velocity_time>1.:
            V[11]=physbam.Vf3(.1,0,0)
            
    def Zero_Out_Enslaved_Velocity_Nodes(self,V,velocity_time,current_position_time,super_fragment):
        self.Zero_Nodes(self.fixed_nodes,V)

    def Preprocess_Frame(self,frame):
        self.adhesion.Update_Springs(True)
        
    def Postprocess_Frame(self,frame):
        self.adhesion.Write_State(physbam.Get_Stream_Type("float"),self.output_directory+"/adhesion.%d"%frame)
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

physbam.LOG.Initialize_Logging(False,False,9999,False)
physbam.Set_Floating_Point_Exception(True,True,True,True,False,False)
#physbam.LOG.Initialize_Logging(False,False,9999,False)
example=EMBEDDING_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
#del example
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

