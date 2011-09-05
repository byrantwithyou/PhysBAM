#!/usr/bin/python

import physbam
import math
def subset(data,items):
    return map(lambda x: data[x],items)

class REPULSION_TEST_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.tests=physbam.SOLIDS_STANDARD_TESTS_Vf3(self)

        self.last_frame=400
        self.solids_parameters.cfl-=4.0
        self.solids_parameters.cg_tolerance=.01
        self.output_directory="repulsion_test"
        self.printed_header=False
        self.g=9.8
        self.panels=8
        self.aspect=2
        self.panel_size=1./2

        # collision parameters
        self.solids_parameters.perform_self_collision=True

        self.solids_parameters.min_collision_loops,self.solids_parameters.max_collision_loops=16,16
        self.solids_parameters.turn_off_all_collisions=True
        self.solids_parameters.check_mesh_for_self_intersection=False
        #self.solids_parameters.clamp_repulsion_thickness=False
        #
        self.solids_parameters.perform_per_collision_step_repulsions=False
        #elf.solids_parameters.perform_per_time_step_repulsions=False
        #elf.solids_parameters.repulsion_thickness=.03
        #elf.solids_parameters.collisions_repulsion_clamp_fraction=1.
        self.solids_parameters.collisions_final_repulsion_limiter_fraction=.05
        self.solids_parameters.repulsions_limiter_fraction=10
        self.solids_parameters.collisions_final_repulsion_youngs_modulus=3e5
        self.solids_parameters.repulsions_youngs_modulus=.1
        
        #elf.solids_parameters.triangle_repulsions.hierarchy_repulsion_thickness_multiplier=2.
        #elf.solids_parameters.collisions_repulsion_spring_multiplier=1.
        #elf.solids_parameters.collisions_repulsion_spring_constant_over_mass_times_length=10000
        #elf.solids_parameters.spring_limiter_fraction=10000.
        #elf.solids_parameters.output_interaction_pairs=True
        # turn off edge/edge
        #self.solids_parameters.triangle_repulsions.compute_edge_edge_friction=False
        self.solids_parameters.triangle_repulsions.compute_edge_edge_inelastic_collision_repulsion=False
        self.solids_parameters.triangle_repulsions.compute_edge_edge_repulsion=False
        # turn off point/face
        #self.solids_parameters.triangle_repulsions.compute_point_face_friction=False
        self.solids_parameters.triangle_repulsions.compute_point_face_inelastic_collision_repulsion=False
        self.solids_parameters.triangle_repulsions.compute_point_face_repulsion=False
        

    def Initialize_Bodies(self):
        # helper references
        deformable_object=self.solids_parameters.deformable_object
        particles=deformable_object.particles
        binding_list=deformable_object.binding_list
        soft_bindings=deformable_object.soft_bindings

        ################################################################
        # Geometry Phase
        ################################################################
        #import random
        #r=random.Random()
        tri=self.tests.Create_Cloth_Panel(self.panels,self.panel_size,self.aspect,physbam.RIGID_BODY_STATE_Vf3(),None)
        #for i in range(1,len(particles)+1):
        #    particles.X[i]+=.02*physbam.Vf3(r.random(),r.random(),r.random())
        tri.Set_Mass_Of_Particles(False,True)
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
        deformable_object.Add_Force(physbam.Create_Edge_Springs(tri,1e2,3.,True,.1,True,0.,True))
        deformable_object.Add_Force(physbam.Create_Bending_Springs(tri,1e2,2.,True,.1,True,0.,True))
        
        # update fragments
        deformable_object.Update_Fragments()

        #self.tests.Add_Ground(.1,-2,0.5)
        #self.tests.Add_Rigid_Body("cylinder",.8,3,True).Frame().t=(0,-.5,0)

        self.Setup_Collisions()

    def Set_External_Velocities(self,V,velocity_time,current_position_time,fragment_id):
        #if velocity_time>1 and velocity_time<1.5: vel=physbam.Vf3(.5,0,0)
        #else: vel=physbam.Vf3(0,0,0)
        #V[1]=vel
        #V[self.panels+1]=-1*vel
        #V[1]=0
        #V[self.panels+1]=0
        m=int(self.aspect*self.panels)+1
        n=self.panels+1;
        i=m/2+1
        for j in range(1,n+1):
            V[i+m*(j-1)]=0

        
    def Zero_Out_Enslaved_Velocity_Nodes(self,V,velocity_time,current_position_time,fragment_id):
        m=int(self.aspect*self.panels)+1
        n=self.panels+1;
        i=m/2+1
        for j in range(1,n+1):
            V[i+m*(j-1)]=0

    def Preprocess_Frame(self,frame): pass
    def Postprocess_Frame(self,frame): print "Completed Frame %d"%frame
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

physbam.LOG.Initialize_Logging(False,True,9999,False)
example=REPULSION_TEST_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
#example.Set_Write_Substeps_Level(1000)
#del example
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

