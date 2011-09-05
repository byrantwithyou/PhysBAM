#!/usr/bin/python

import physbam
import math
def subset(data,items):
    return map(lambda x: data[x],items)

class EMBEDDING_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.last_frame=400
        self.solids_parameters.cfl=4
        self.solids_parameters.perform_self_collision=False
        self.output_directory="cloth"
        self.printed_header=False
        self.g=9.8
        self.panels,self.aspect=80,1
        self.tests=physbam.SOLIDS_STANDARD_TESTS_Vf3(self)

    def Initialize_Bodies(self):
        # helper references
        deformable_object=self.solids_parameters.deformable_object
        particles=deformable_object.particles
        binding_list=deformable_object.binding_list
        soft_bindings=deformable_object.soft_bindings

        ################################################################
        # Geometry Phase
        ################################################################
        tri=self.tests.Create_Cloth_Panel(self.panels,1,self.aspect,physbam.RIGID_BODY_STATE_Vf3(),None)
        tri.Set_Mass_Of_Particles(True,True)
        tri.Update_Number_Nodes()

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
        deformable_object.Add_Force(physbam.GRAVITY_Vf3(particles,self.g,physbam.Vf3(0,-1,0)))
        deformable_object.Add_Force(physbam.Create_Edge_Springs(tri,1e3,3.,True,.1,True,0.,True))
        deformable_object.Add_Force(physbam.Create_Bending_Springs(tri,1e2,2.,True,.1,True,0.,True))
        
        # update fragments
        deformable_object.Update_Fragments()

        #self.tests.Add_Ground(.1,-2,0.5)
        self.tests.Add_Rigid_Body("sphere",.2,0,True).Frame().t=(0,-.5,-.5)

        self.Setup_Collisions()

    def Set_External_Velocities(self,V,velocity_time,current_position_time,fragment_id):
        if velocity_time>1 and velocity_time<1.5: vel=physbam.Vf3(.5,0,0)
        else: vel=physbam.Vf3(0,0,0)
        V[1]=vel
        if velocity_time<2.5:
            V[self.panels+1]=-1*vel
        
    def Zero_Out_Enslaved_Velocity_Nodes(self,V,velocity_time,current_position_time,fragment_id):
        V[1]=0
        if velocity_time<2.5:
            V[self.panels+1]=0

    def Preprocess_Frame(self,frame): pass
    def Postprocess_Frame(self,frame): print frame
    def Preprocess_Solids_Substep(self,time,substep): pass
    def Postprocess_Solids_Substep(self,time,substep): pass
    def Update_Solids_Parameters(self,dt): pass
    def Add_External_Forces(self,F,time,fragment_id): pass
    def Add_External_Impulses(self,V,time,dt,fragment_id): pass
    def Add_External_Impulse(self,V,node,time,dt,fragment_id): pass
    def Align_Deformable_Bodies_With_Rigid_Bodies(self,*f): pass
    def Apply_Constraints(self,dt,time): pass
    def Update_Time_Varying_Material_Properties(self,time,fragment_id): pass
    def Set_External_Positions(self,X,time,fragment_id): pass
    def Update_Collision_Body_Positions_And_Velocities(self,time): pass
    def Limit_Solids_Dt(self,time): pass

physbam.LOG.Initialize_Logging(True,True,0,False)
example=EMBEDDING_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
#del example
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

