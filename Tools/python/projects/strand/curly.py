#!/usr/bin/python

import physbam
import math
def subset(data,items):
    return map(lambda x: data[x],items)

class EMBEDDING_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.tests=physbam.SOLIDS_STANDARD_TESTS_Vf3(self)


        self.last_frame=8000
        self.frame_rate=24
        self.solids_parameters.cfl=.9
        self.solids_parameters.perform_self_collision=False
        self.output_directory="output"
        self.printed_header=False
        self.g=32
        self.Unite_Fragments()
        self.data_directory="/data/aselle/hair_modeling"

    def Initialize_Bodies(self):
        # helper references
        deformable_object=self.solids_parameters.deformable_object
        particles=deformable_object.particles
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
        for i in range(1,5):
            self.fixed_nodes.Append(i)

        #physbam.Read_From_File("float","fixed_nodes",self.fixed_nodes)
        curve=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        torsion_curve=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        n_segments=50
        particles.Add_Particles(n_segments+1)
        print "len(particles)=%d"%len(particles)
        for i in range(1,n_segments+1):
            s=float(i-1)/(n_segments+1)
            particles.X[i]=X_spiral(s)
            print "(%d,%d)"%(i,i+1)
            curve.mesh.elements.Append(i,i+1)
            if i+3<=n_segments:
                torsion_curve.mesh.elements.Append(i,i+3)
        particles.X[n_segments+1]=X_spiral(1)
        curve.density=.1
        curve.Set_Mass_Of_Particles(True,True)
        curve.Update_Number_Nodes()
        torsion_curve.Update_Number_Nodes()
        
        volume=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)
        for i in range(1,n_segments-1):
            def tet_volume(e):
                return physbam.TETRAHEDRON_f(*map(lambda i:particles.X[i],e)).Signed_Volume()
            element=(i,i+1,i+2,i+3)
            if not tet_volume(element)>=0:
                element=(i,i+1,i+3,i+2)
            print "tet %s volume=%f"%(repr(element),tet_volume(element))
            volume.mesh.elements.Append(*element)
        volume.Update_Number_Nodes()
        

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
        deformable_object.Add_Force(physbam.Create_Edge_Springs(curve,1e5,2.,True,.1,True,0.,True))
        deformable_object.Add_Force(physbam.Create_Edge_Springs(torsion_curve,1e4,2.,True,.1,True,0.,True))
        deformable_object.Add_Force(physbam.Create_Segment_Bending_Springs(curve,1e6,2.,True,.1,True,0.,True))
        deformable_object.Add_Force(physbam.Create_Tet_Springs(volume,1e4,2.,False,.1,True,.1,True,0.,True))
        #deformable_object.Add_Force(physbam.Create_Altitude_Springs(volume,1e4,20.,False,.1,True,.1,True,0.,True))

        deformable_object.Add_Structure(curve)
        deformable_object.Add_Structure(volume)
        deformable_object.Add_Structure(torsion_curve)
        
        # update fragments
        deformable_object.Update_Fragments()

        #self.tests.Add_Ground(.1,-2,0.5)
        #self.tests.Add_Rigid_Body("efty_100k_nose_filled",1,1,True)
        #center of mass = -0.00479698 14.673 8.13414
        #orientation = 0.999839 -0.0162252 -0.00339366 0.00689811
    

        self.Setup_Collisions()

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

    def Preprocess_Frame(self,frame): pass
    def Postprocess_Frame(self,frame): print frame
    def Preprocess_Solids_Substep(self,time,dt,substep):  pass
    def Postprocess_Solids_Substep(self,time,dt,substep):  pass
    def Update_Solids_Parameters(self,dt): pass
    def Postprocess_Solids_Substep(self,time,dt,substep): pass
    def Add_External_Forces(self,F,time,fragment_id): pass
    def Add_External_Impulses(self,V,time,dt,fragment_id): pass
    def Add_External_Impulse(self,V,node,time,dt,fragment_id): pass
    def Align_Deformable_Bodies_With_Rigid_Bodies(self,*f): pass
    def Apply_Constraints(self,dt,time): pass
    def Update_Time_Varying_Material_Properties(self,time,fragment_id): pass
    def Set_External_Positions(self,X,time,fragment_id): pass
    def Update_Collision_Body_Positions_And_Velocities(self,time): pass
    def Limit_Solids_Dt(self): pass

physbam.LOG.Initialize_Logging(True,True,0,False)
#physbam.LOG.Initialize_Logging(False,False,9999,False)
example=EMBEDDING_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
#del example
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

