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
        self.output_directory="bend"

    def Initialize_Bodies(self):
        deformable_object=self.solids_parameters.deformable_object
        particles=deformable_object.particles

        filename=os.path.join(os.environ["PHYSBAM"],"Public_Data","Triangulated_Surfaces/torus.tri") # ,"Rigid_Bodies/Thin_Shells/boat_medres.tri")
        rot=physbam.ROTATION_Vf3(.8,physbam.Vf3(1,0,0))
        frame=physbam.FRAME_Vf3(physbam.Vf3(),rot)
        tri=self.tests.Create_Triangulated_Object(filename,physbam.RIGID_BODY_STATE_Vf3(frame),False,True,1.)

        tet=physbam.TETRAHEDRALIZED_VOLUME_f.Create(tri.particles)
        tet.mesh.Initialize_Bending_Tetrahedrons(tri.mesh)
        # swap ordering if signed volume negative
        for i in tet.mesh.elements:
            t=physbam.TETRAHEDRON_f(*map(lambda i: tet.particles.X[i],i))
            if t.Signed_Volume()<0: i[2],i[3]=i[3],i[2]
        tet.Update_Number_Nodes()


        particles.mass.Compute_Auxiliary_Attributes(deformable_object.soft_bindings)

        deformable_object.Add_Force(physbam.GRAVITY_Vf3(particles,9.8,(0,-1,0)))

        deformable_object.Add_Force(physbam.Create_Edge_Springs(tri,2e2,3.,False,.1,True,0.,True))
        deformable_object.Add_Force(physbam.Create_Bending_Springs(tri,2e3,3.,False,.1,True,0.,True))
        self.tet_force_id=deformable_object.Add_Force(physbam.Create_Tet_Springs(tet,2e1,3.,False,.1,True,.1,True,0.,True))
        deformable_object.Add_Structure(tet)

        deformable_object.Update_Fragments()

        self.tests.Add_Ground(.1,-1,0.5)
        self.Setup_Collisions()
        
                                   

                                          

example=BEND_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

