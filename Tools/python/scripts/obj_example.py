#!/usr/bin/python
import os
import physbam

class OBJ_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.tests=physbam.SOLIDS_STANDARD_TESTS_Vf3(self)

        self.last_frame=100
        print "FUCK FUCK FUCK %f" % self.solids_parameters.cfl
        self.solids_parameters.cfl=30
        print "FUCK FUCK FUCK %f" % self.solids_parameters.cfl
        self.solids_parameters.verbose_dt=True
        self.solids_parameters.perform_self_collision=False
        self.output_directory="output_tri"
        
        #self.frame_rate=48
        pass
    def Initialize_Bodies(self):
        # helper references
        deformable_object=self.solids_parameters.deformable_object
        particles=deformable_object.particles
        # make geometry

        geom=[]
        count=0
        for f in ["Rigid_Bodies/Thin_Shells/boat_medres.tri"]: #,"Rigid_Bodies/Thin_Shells/boat_hires.tri"]:
            filename=os.path.join(os.environ["PHYSBAM"],"Public_Data",f)
            tri_temp=physbam.TRIANGULATED_SURFACE_f.Create()
            physbam.Read_From_File("float",filename,tri_temp)
            indices=physbam.LA_i()
            for i in range(1,len(tri_temp.particles)+1):
                tri_temp.particles.X[i]+=count*physbam.Vf3(1.5,0,0)
            geom.append(tri_temp.Append_Particles_And_Create_Copy(particles,indices))
            count+=1

        # update # nodes and masses
        for g in geom:
            g.Update_Number_Nodes()
            g.Set_Mass_Of_Particles(False,True)
        particles.mass.Compute_Auxiliary_Attributes(deformable_object.soft_bindings)

        # Create forces
        deformable_object.Add_Force(physbam.GRAVITY_Vf3(particles,9.8,(0,-1,0)))
        for g in geom:
            deformable_object.Add_Force(physbam.Create_Edge_Springs(g,2e5,3.,False,.1,True,0.,True))
            deformable_object.Add_Force(physbam.Create_Bending_Springs(g,2e5,3.,False,.1,True,0.,True))

        for g in geom:
            self.solids_parameters.deformable_object.Add_Structure(g)

        deformable_object.Update_Fragments()

        self.tests.Add_Ground(.1,-2,0.5)
        self.Setup_Collisions()


    def Set_External_Velocities(self,V,velocity_time,current_position_time,fragment_id):
        V[1]=physbam.Vf3()
        
    def Zero_Out_Enslaved_Velocity_Nodes(self,V,velocity_time,current_position_time,fragment_id):
        V[1]=physbam.Vf3()

    def Preprocess_Frame(self,frame): pass
    def Postprocess_Frame(self,frame): pass

example=OBJ_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

