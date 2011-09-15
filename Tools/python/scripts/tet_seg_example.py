#!/usr/bin/python

import physbam

class TET_SEG_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.last_frame=200
        self.solids_parameters.cfl=4.0
        self.solids_parameters.perform_self_collision=False
        self.output_directory="output_tet_seg"
        
        #self.frame_rate=48
        pass
    def Initialize_Bodies(self):
        # helper references
        deformable_object=self.solids_parameters.deformable_object
        particles=deformable_object.particles
        # make geometry
        tet=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)
        tet.particles.Add_Particles(4)
        tet.particles.X[1]=(0,0,0);tet.particles.X[2]=(1,0,0);tet.particles.X[3]=(0,1,0);tet.particles.X[4]=(0,0,1)
        tet.mesh.elements.Append(1,2,3,4)
        seg=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        particles.Add_Particle();particles.X[5]=(-1,0,0);
        seg.mesh.elements.Append(4,5)
        for i in range(1,len(tet.particles)+1): tet.particles.mass[i]=1
        seg.Update_Number_Nodes();tet.Update_Number_Nodes()

        # figure out masses
        tet.Set_Mass_Of_Particles(False,True);
        particles.mass.Compute_Auxiliary_Attributes(deformable_object.soft_bindings)

        # Create forces
        deformable_object.Add_Force(physbam.GRAVITY_Vf3(particles,9.8,(0,-1,0)))
        deformable_object.Add_Force(physbam.Create_Edge_Springs(tet,2e3,1.,False,.1,True,0.,True))
        deformable_object.Add_Force(physbam.Create_Edge_Springs(seg,2e3,1.,False,.1,True,0.,True))
        
        # add structures
        self.solids_parameters.deformable_object.Add_Structure(tet)
        self.solids_parameters.deformable_object.Add_Structure(seg)

        # update fragments
        deformable_object.Update_Fragments()

    def Set_External_Velocities(self,V,velocity_time,current_position_time,fragment_id):
        if velocity_time>2: V[5]=(0,1,0)
        else: V[5]=(0,0,0)
        if velocity_time<1: V[1]=(0,0,0)
        
        
    def Zero_Out_Enslaved_Velocity_Nodes(self,V,velocity_time,current_position_time,fragment_id):
        V[5]=(0,0,0)
        if velocity_time<1: V[1]=(0,0,0)

    def Preprocess_Frame(self,frame):
        pass

    def Postprocess_Frame(self,frame):
        pass
        #print "Frame %d"%frame

physbam.LOG.Initialize_Logging(False,False,2,False)
example=TET_SEG_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

