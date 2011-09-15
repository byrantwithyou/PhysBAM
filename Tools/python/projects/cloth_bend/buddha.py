#!/usr/bin/python
import os
import physbam

class BEND_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.tests=physbam.SOLIDS_STANDARD_TESTS_Vf3(self)

        self.last_frame=100
        self.solids_parameters.cfl=8
        self.solids_parameters.verbose_dt=True
        self.solids_parameters.perform_self_collision=False
        self.output_directory="buddha"

    def Initialize_Bodies(self):
        deformable_object=self.solids_parameters.deformable_object
        particles=deformable_object.particles
        rigid_particles=deformable_object.rigid_body_particles

        filename=os.path.join(os.environ["PHYSBAM"],"Public_Data","Triangulated_Surfaces/buddha_fine.tri") # ,"Rigid_Bodies/Thin_Shells/boat_medres.tri")
        rot=physbam.ROTATION_Vf3(.8,physbam.Vf3(1,0,0))
        frame=physbam.FRAME_Vf3(physbam.Vf3(),rot)
        tri=self.tests.Create_Triangulated_Object(filename,physbam.RIGID_BODY_STATE_Vf3(frame),False,True,1.)

        print "initializing tet mesh"
        tet=physbam.TETRAHEDRALIZED_VOLUME_f.Create(tri.particles)
        tet.mesh.Initialize_Bending_Tetrahedrons(tri.mesh)
        # swap ordering if signed volume negative
        for i in tet.mesh.elements:
            t=physbam.TETRAHEDRON_f(*map(lambda i: tet.particles.X[i],i))
            if t.Signed_Volume()<0: i[2],i[3]=i[3],i[2]
        tet.Update_Number_Nodes()
        print "done."


        particles.mass.Compute_Auxiliary_Attributes(deformable_object.soft_bindings)

        deformable_object.Add_Force(physbam.GRAVITY_Vf3(particles,rigid_particles,9.8,(0,-1,0)))

        deformable_object.Add_Force(physbam.Create_Edge_Springs(tri,2e1,3.,False,.1,True,0.,True))
        deformable_object.Add_Force(physbam.Create_Bending_Springs(tri,2e2,3.,False,.1,True,0.,True))
        springs=physbam.Create_Tet_Springs(tet,2e2,3.,False,.1,True,.1,True,0.,True)
        springs.Clamp_Restlength(.005);
        self.tet_force_id=deformable_object.Add_Force(springs)
        deformable_object.Add_Structure(tet)

        self.tests.Add_Ground(.1,-1,0.5)
        deformable_object.Update_Fragments()

        self.Setup_Collisions()
        
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
                                   

#for frame in range(1,100):
#        last_frame=100
#        percent=float(frame)/float(last_frame)*100.
#        percent_int=int(percent)
#        foo="["+"*"*percent_int+" "*(100-percent_int)+"]"
#        print "\r%5d of %5d %.1f  %s"%(frame,last_frame,percent,foo)

physbam.LOG.Initialize_Logging(True,True,9999,False)
example=BEND_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

print "\nDone."
