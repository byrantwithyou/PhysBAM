#!/usr/bin/python

import physbam
import math
def subset(data,items):
    return map(lambda x: data[x],items)

class EMBEDDING_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.last_frame=800
        self.solids_parameters.cfl=8
        self.solids_parameters.perform_self_collision=False
        self.output_directory="hair"
        self.printed_header=False
        self.g=32
        self.tests=physbam.SOLIDS_STANDARD_TESTS_Vf3(self)
        #self.Unite_Fragments()
        self.data_directory="/data/mlentine/hair_modeling"

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
        curve=self.tests.Create_Segmented_Curve("hairs.curve.gz",physbam.RIGID_BODY_STATE_Vf3(),False,False)
        self.fixed_nodes=physbam.LA_i()
        physbam.Read_From_File("float","fixed_nodes",self.fixed_nodes)

        curve.density=1
        curve.Set_Mass_Of_Particles(True,True)
        curve.Update_Number_Nodes()

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
        deformable_object.Add_Force(physbam.GRAVITY_Vf3(particles,rigid_particles,self.g,True,True,physbam.Vf3(0,-1,0)))
        deformable_object.Add_Force(physbam.Create_Edge_Springs(curve,1e5,3.,True,.1,True,0.,True))
        deformable_object.Add_Force(physbam.Create_Segment_Bending_Springs(curve,1e5,3.,True,.1,True,0.,True))
        
        # update fragments
        deformable_object.Update_Fragments()

        #self.tests.Add_Ground(.1,-2,0.5)
        self.tests.Add_Rigid_Body("efty_100k_nose_filled",1,1,True)
        #center of mass = -0.00479698 14.673 8.13414
        #orientation = 0.999839 -0.0162252 -0.00339366 0.00689811
    

        #self.Setup_Collisions()

    def Set_External_Velocities(self,V,velocity_time,current_position_time,fragment_id):
        #if velocity_time>1 and velocity_time<1.5: vel=physbam.Vf3(.5,0,0)
        #else: vel=physbam.Vf3(0,0,0)
        #V[1]=vel
        #V[self.panels+1]=-1*vel
        self.Zero_Nodes(self.fixed_nodes,V)
        
    def Zero_Out_Enslaved_Velocity_Nodes(self,V,velocity_time,current_position_time,fragment_id):
        #V[1]=0
        #V[self.panels+1]=0
        self.Zero_Nodes(self.fixed_nodes,V)

    def Add_External_Fragment_Connectivity(self,union_find,particle_is_simulated):
        self.Unite_All_Fragments(union_find,particle_is_simulated)

    def Preprocess_Frame(self,frame): pass
    def Postprocess_Frame(self,frame): pass
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


#physbam.LOG.Initialize_Logging(True,True,0,False)
physbam.LOG.Initialize_Logging(False,False,9999,False)
example=EMBEDDING_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
#del example
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

