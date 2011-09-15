#!/usr/bin/python

import physbam
import math
def subset(data,items):
    return map(lambda x: data[x],items)

class EMBEDDING_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.last_frame=400
        self.solids_parameters.cfl=1
        self.solids_parameters.perform_self_collision=False
        self.output_directory="output_tet"
        self.printed_header=False

        # my paramemters
        self.use_zero_length_implicit=True # If false then use linear for both
        self.g=9.8
        self.ym=50
        self.restlength_clamp=.1

    def Initialize_Bodies(self):
        # helper references
        deformable_object=self.solids_parameters.deformable_object
        particles=deformable_object.particles
        binding_list=deformable_object.binding_list
        soft_bindings=deformable_object.soft_bindings

        ################################################################
        # Geometry Phase
        ################################################################
        tet=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)
        tri=physbam.TRIANGULATED_SURFACE_f.Create(particles)
        curve=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        # triangle points
        particles.Add_Particles(3)
        particles.X[1]=physbam.Vf3(-0.668156,1.666096,1.202544);
        particles.X[2]=physbam.Vf3(0.104949,0.717676,0.022739);
        particles.X[3]=(0.306855,0.963611,1.181754)
        tri_points=1,2,3
        for i in tri_points: particles.mass[i]=0 #.3333334
        # tet points
        particles.Add_Particles(4)
        particles.X[4]=(-0.433013,-0.018620,-0.750000);particles.X[5]=(-0.433013,-0.018620,0.750000);particles.X[6]=(0.866025,-0.018620,0.000000);
        particles.X[7]=(0.000000,1.000000,0.000000);
        self.tet_points=4,5,6,7
        for i in self.tet_points: particles.mass[i]=.25
        tet.mesh.elements.Append(self.tet_points)
        tetrahedron=physbam.TETRAHEDRON_f(*subset(particles.X,self.tet_points))
        # embed point 3
        self.embedded_index=2
        particles.X[self.embedded_index]=(0,1,0)
        binding_list.Add_Binding(physbam.LINEAR_BINDING_Vf3_4(particles,self.embedded_index,physbam.Vi4(*self.tet_points),tetrahedron.Barycentric_Coordinates(particles.X[self.embedded_index])))
        self.soft_embedded_index=particles.Add_Particle()
        self.target_embedded_index=particles.Add_Particle()
        particles.X[self.target_embedded_index]=particles.X[self.soft_embedded_index]=particles.X[self.embedded_index]
        #particles.mass[self.target_embedded_index]=particles.mass[self.soft_embedded_index]=particles.mass[self.embedded_index]=1
        particles.mass[self.target_embedded_index]=1000

        soft_tri_points=1,self.soft_embedded_index,3
        tri.mesh.elements.Append(soft_tri_points)

        #for i in range(1,len(particles)+1):
        #    print "Mass(%d)=%f"%(i,particles.mass[i])

        # make a curve

        curve.mesh.elements.Append(self.soft_embedded_index,self.target_embedded_index)
        if not self.use_zero_length_implicit:
            curve.mesh.elements.Append(self.embedded_index,self.soft_embedded_index)
        else:
            soft_bindings.Add_Binding((self.soft_embedded_index,self.embedded_index),False)
        # TODO: START HERE ------------>
        #soft_bindings.Append(

        geoms=[tet,tri,curve]
        map(lambda x: x.Update_Number_Nodes(),geoms)

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
        for i in self.tet_points: tet_node_list.Append(i)
        deformable_object.Add_Force(physbam.GRAVITY_Vf3(particles,tet_node_list,self.g,physbam.Vf3(0,-1,0)))
        springs=physbam.Create_Edge_Springs(curve,self.ym,2.,False,.1,True,0.,True)
        springs.Clamp_Restlength(self.restlength_clamp)
        springs.Set_Overdamping_Fraction(1)
        deformable_object.Add_Force(springs)
        deformable_object.Add_Force(physbam.Create_Finite_Volume(tet,physbam.NEO_HOOKEAN_f3(3e6,.475,.05,.25)))
        if self.use_zero_length_implicit:
            soft_bindings.Initialize_Binding_Mesh(False)
            deformable_object.Add_Force(physbam.Create_Edge_Binding_Springs(particles,soft_bindings.binding_mesh,self.ym/self.restlength_clamp,4.,True)) # note we use ym/restlength_clamp as we are zero length

        # add structures
        map(deformable_object.Add_Structure,geoms)

        # update fragments
        deformable_object.Update_Fragments()

    def Set_External_Velocities(self,V,velocity_time,current_position_time,fragment_id):
        V[self.target_embedded_index]=0
        
    def Zero_Out_Enslaved_Velocity_Nodes(self,V,velocity_time,current_position_time,fragment_id):
        V[self.target_embedded_index]=0

    def Preprocess_Frame(self,frame): pass
    def Postprocess_Frame(self,frame):
        if frame%10==0:
            particles=self.solids_parameters.deformable_object.particles
            binding_list=self.solids_parameters.deformable_object.binding_list
            mass=reduce(lambda x,y:x+y,subset(particles.mass,self.tet_points))
            analytic_distance=self.restlength_clamp*mass*self.g/self.ym
            if self.use_zero_length_implicit:
                numeric_distance_soft=(particles.X[self.embedded_index]-particles.X[self.soft_embedded_index]).Magnitude()
                numeric_distance_target=(particles.X[self.soft_embedded_index]-particles.X[self.target_embedded_index]).Magnitude()
                analytic_distance2=self.restlength_clamp/self.ym*(particles.mass[self.soft_embedded_index]*self.g*mass/particles.mass[7]+self.g*mass)
                if not self.printed_header:
                    print "%5s %32s | %32s"%("------","------------ soft -----------","---------------- hard -----------")
                    print "%5s %10s %10s %10s | %10s %10s %10s"%("frame","numeric","analytic","error","numeric","analytic","error")
                print "%5d %10.4f %10.4f %10.4g | %10.4f %10.4f %10.4g"%(frame,numeric_distance_soft,analytic_distance,numeric_distance_soft-analytic_distance,numeric_distance_target,analytic_distance2,numeric_distance_target-analytic_distance2)

            else:
                numeric_distance=(particles.X[self.embedded_index]-particles.X[self.soft_embedded_index]).Magnitude()
                print "%10.4g %10.4g %10.6g"%(analytic_distance,numeric_distance,math.fabs(numeric_distance-analytic_distance))
                
            self.printed_header=True

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

physbam.LOG.Initialize_Logging(True,True,0,False)
example=EMBEDDING_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

