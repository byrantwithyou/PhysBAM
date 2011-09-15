#!/usr/bin/python
import sys
import os
import physbam
import math
import shutil

class HAIR_EXAMPLE(physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3):
    def __init__(self,stream_type,regions,fluid_type):
        physbam.SOLIDS_FLUIDS_EXAMPLE_UNIFORM_f3.__init__(self,stream_type,regions,fluid_type)
        self.tests=physbam.SOLIDS_STANDARD_TESTS_Vf3(self)

        self.last_frame=1000
        self.frame_rate=24
        self.solids_parameters.cfl=6

        if sys.argv<4:
            print "Usage: %s <modeldir> <sim> {collide|nocollide}"
            sys.exit(0)

        collide_str=sys.argv[3]
        if collide_str=="collide": self.solids_parameters.perform_self_collision=True
        elif collide_str=="nocollide": self.solids_parameters.perform_self_collision=False
        else:
            print "Usage: %s <modeldir> <sim> {collide|nocollide}"
            sys.exit(0)
        
        self.data_directory=sys.argv[1] # "/data/"+os.environ["USER"]+"/PhysBAM/Public_Data"
        self.model_directory=os.path.join(sys.argv[1],sys.argv[2])
        self.output_directory="output_%s_%s_%s"%(sys.argv[1].split("/")[-1],sys.argv[2],collide_str)
        self.printed_header=False
        self.g=9.8
        #self.Unite_Fragments()
        self.solids_parameters.use_push_out=False

        if not os.path.exists(self.output_directory): os.mkdir(self.output_directory)
        shutil.copy(sys.argv[0],os.path.join(self.output_directory,"sim.py"))

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
        print "Reading data...",
        
        self.fixed_nodes=physbam.LA_i()
        self.fixed_positions=physbam.LA_Vf3()
        self.edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        self.extra_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        self.bending_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        self.torsion_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
        self.volume=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)
        physbam.Read_From_File("float",self.model_directory+"/fixed_nodes",self.fixed_nodes)
        physbam.Read_From_File("float",self.model_directory+"/edges.curve",self.edges.mesh)
        physbam.Read_From_File("float",self.model_directory+"/extra_edges.curve",self.extra_edges.mesh)
        physbam.Read_From_File("float",self.model_directory+"/bending_edges.curve",self.bending_edges.mesh)
        physbam.Read_From_File("float",self.model_directory+"/torsion_edges.curve",self.torsion_edges.mesh)
        physbam.Read_From_File("float",self.model_directory+"/torsion_edges.curve",self.torsion_edges.mesh)
        physbam.Read_From_File("float",self.model_directory+"/tets.gz",self.volume.mesh)
        physbam.Read_From_File("float",self.model_directory+"/particles.gz",particles)
        particles.Store_Velocity(True)
        print "     done reading data..."

        # Store object space positions
        for i in self.fixed_nodes:
            self.fixed_positions.Append(particles.X[i])


        # Set masses
        #self.edges.Set_Mass_Of_Particles(True,True)
        #self.extra_edges.Set_Mass_Of_Particles(True,True)
        masses=physbam.LA_f()
        physbam.Read_From_File("float",self.model_directory+"/masses.gz",masses)
        for i in range(1,len(particles)+1):
            #print "i=%d len(p.m)=%d len(masses.m)=%d"%(i,len(particles.mass),len(masses))
            particles.mass[i]=masses[i]
            

        # Set number nodes
        for i in [self.edges,self.extra_edges,self.bending_edges,self.torsion_edges]:
            i.Update_Number_Nodes()
        # TODO use more sane mass
        for i in range(1,len(particles)+1):
            particles.mass[i]=1
        
        # Orient tetrahedra
        print "Orienting tets",
        i=1
        for nodes in self.volume.mesh.elements:
            if physbam.TETRAHEDRON_f(*map(lambda i:particles.X[i],nodes)).Signed_Volume()<0:
                nodes[1],nodes[2],nodes[3],nodes[4]=(nodes[1],nodes[2],nodes[4],nodes[3])
            i+=1
        self.volume.Update_Number_Nodes()
        print "done."
        
        ################################################################
        # Mass Update Phase
        ################################################################
        print "Updating mass attributes"
        binding_list.Distribute_Mass_To_Parents()
        binding_list.Clear_Hard_Bound_Particles(particles.mass)
        particles.mass.Compute_Auxiliary_Attributes(soft_bindings)
        deformable_object.soft_bindings.Set_Mass_From_Effective_Mass()

        ################################################################
        # Forces and Fragme1nt Connectivity
        ################################################################
        print "Creating forces"
        tet_node_list=physbam.LA_i()
        deformable_object.Add_Force(physbam.GRAVITY_Vf3(particles,rigid_particles,True,True,self.g,physbam.Vf3(0,-1,0)))
        overdamping_fraction=40
        springs=[
            physbam.Create_Edge_Springs(self.edges,1e3,overdamping_fraction,True,.1,True,0.,True,False),
            physbam.Create_Edge_Springs(self.extra_edges,1e2,overdamping_fraction,True,.1,True,0.,True,False),
            physbam.Create_Edge_Springs(self.bending_edges,1e2,overdamping_fraction,True,.1,True,0.,True,False),
            physbam.Create_Edge_Springs(self.torsion_edges,1e2,overdamping_fraction,True,.1,True,0.,True,False),
            physbam.Create_Tet_Springs(self.volume,1e2,overdamping_fraction,False,.1,True,.1,True,0.,True,False)
            ]
        # clamp restlengths
        for force in springs:
            force.Clamp_Restlength(.01)
            deformable_object.Add_Force(force)
        # adhesion
        self.adhesion=None
        self.use_adhesion=True
        #self.adhesion=physbam.SEGMENT_ADHESION_Vf3(particles,self.edges.mesh)

        if self.use_adhesion:
            self.stiction_mesh=physbam.SEGMENT_MESH()
            
            union_find=physbam.UNION_FIND_int(len(particles))
            for e in self.edges.mesh.elements:
                self.stiction_mesh.elements.Append(e)
            for e in self.volume.mesh.elements:
                union_find.Union(e[1],e[2])
                union_find.Union(e[2],e[3])
                union_find.Union(e[1],e[2])
            
            self.particle_to_spring=physbam.LA_HAIR_ID()
            self.particle_to_spring.Resize(len(particles))
            next_hair=1
            for i in xrange(1,len(particles)+1):
                parent=union_find.Find(i)
                if self.particle_to_spring[parent]==0:
                    self.particle_to_spring[parent]=physbam.HAIR_ID(next_hair)
                    next_hair+=1
                self.particle_to_spring[i]=self.particle_to_spring[parent]
                
            
            self.adhesion=physbam.SEGMENT_ADHESION_Vf3(particles,self.stiction_mesh,self.particle_to_spring)
            # high volume
            #minimum_width=.01
            #maximum_width=.2
            #self.adhesion.Set_Parameters(5.,10.,minimum_width,maximum_width,50) # ym,overdamp,start dist, stop dist
            #self.adhesion.Set_Restlength(.2)
            # medium volume
            #minimum_width=.01
            #maximum_width=.1
            #self.adhesion.Set_Parameters(5.,10.,minimum_width,maximum_width,50) # ym,overdamp,start dist, stop dist
            #self.adhesion.Set_Restlength(.1)
            # small volume
            #minimum_width=.01
            #maximum_width=.05
            #self.adhesion.Set_Parameters(5.,10.,minimum_width,maximum_width,50) # ym,overdamp,start dist, stop dist
            #self.adhesion.Set_Restlength(.05)

            # small volume
            minimum_width=.01
            maximum_width=.1
            self.adhesion.Set_Parameters(10.,10.,minimum_width,maximum_width,50) # ym,overdamp,start dist, stop dist
            self.adhesion.Set_Restlength(.05)

            #uber
            #minimum_width=.01
            #maximum_width=.1
            #self.adhesion.Set_Parameters(10.,10.,minimum_width,maximum_width,50) # ym,overdamp,start dist, stop dist
            #self.adhesion.Set_Restlength(.4)
            self.changed_stiction_parameters=False

            deformable_object.Add_Force(self.adhesion)
            self.adhesion.Write_State(physbam.Get_Stream_Type("float"),self.output_directory+"/adhesion.%d"%0)
        
#        if self.adhesion:
#            self.adhesion.Set_Parameters(.01,10.,.02,.05) # ym,overdamp,start dist, stop dist
#            deformable_object.Add_Force(self.adhesion)
            #self.adhesion.Write_State(physbam.Get_Stream_Type("float"),self.output_directory+"/adhesion.%d"%0)
        # update fragments
        deformable_object.Update_Fragments()

        ################################################################
        # Collision Bodies
        ################################################################
        #self.tests.Add_Ground(.1,-2,0.5)
        self.body=self.tests.Add_Rigid_Body("sphere",1,.1,True,False)
        self.body.is_kinematic=True

        self.hand=self.tests.Add_Rigid_Body("hand",1,.1,True,False)
        self.hand.Frame().t=physbam.Vf3(0,-1,1)
        self.hand.Frame().r=physbam.ROTATION_Vf3(3*math.pi/4,physbam.Vf3(0,1,0))
        self.hand.is_static=True


        # Collisions
        self.Collide_Segmented_Curve(self.edges)
        self.Collide_Rigid_Bodies();

        ################################################################
        # Add Geometry to Deformable object (last to prevent early invalid)
        ################################################################
        # Add structures last
        deformable_object.Add_Structure(self.edges)
        if self.adhesion: self.adhesion.Update_Partitions(False,None,self.output_directory)
    
        deformable_object.Add_Structure(self.extra_edges)
        #deformable_object.Add_Structure(self.bending_edges)
        #deformable_object.Add_Structure(self.torsion_edges)
        #deformable_object.Add_Structure(self.volume)
        #deformable_object.Add_Structure(torsion_curve)
        
    def Set_External_Velocities(self,V,velocity_time,current_position_time):
        self.Set_Rigid_Velocities(self.fixed_nodes,self.fixed_positions,self.body,V)
        
    def Zero_Out_Enslaved_Velocity_Nodes(self,V,velocity_time,current_position_time,super_fragment):
        self.Zero_Nodes(self.fixed_nodes,V)

    def Set_External_Positions(self,X,time):
        self.Set_Rigid_Positions(self.fixed_nodes,self.fixed_positions,self.body,X)

    def Set_Kinematic_Positions(self,frame,time,id):
        #if time<2: frame.r=physbam.ROTATION_Vf3.From_Rotation_Vector(time*physbam.Vf3(0,1,0))
        #elif time<4: frame.r=physbam.ROTATION_Vf3.From_Rotation_Vector(2*physbam.Vf3(0,1,0))*physbam.ROTATION_Vf3.From_Rotation_Vector((time-2.)*physbam.Vf3(0,-1,0))
        #else:        frame.r=physbam.ROTATION_Vf3.From_Rotation_Vector(2*physbam.Vf3(0,1,0))*physbam.ROTATION_Vf3.From_Rotation_Vector(2*physbam.Vf3(0,-1,0))*physbam.ROTATION_Vf3.From_Rotation_Vector((time-4.)*physbam.Vf3(1,0,0))
        return
        
        
    def Set_Kinematic_Velocities(self,twist,time,id):
        #if time<2: twist.angular=physbam.Vf3(0,1,0)
        #elif time<4: twist.angular=physbam.Vf3(0,-1,0)
        #else: twist.angular=physbam.Vf3(1,0,0)
        return True
        
    def Postprocess_Frame(self,frame):
        #if self.adhesion: self.adhesion.Write_State(physbam.Get_Stream_Type("float"),self.output_directory+"/adhesion.%d"%frame)
        percent=float(frame)/float(self.last_frame)*100.
        percent_int=int(percent)/2
        print "\r%5d of %5d (%3.1f%%) %s"%(frame,self.last_frame,percent,("["+"*"*percent_int+" "*(50-percent_int)+"]"))

    def Add_External_Fragment_Connectivity(self,union_find,particle_is_simulated):
        self.Unite_All_Fragments(union_find,particle_is_simulated)

    def Add_External_Impulses_Before(self,V,time,dt):
        pass

    def Add_External_Impulses(self,V,time,dt):
        pass

    def Add_External_Impulse(self,V,node,time,dt):
        pass

    def Preprocess_Frame(self,frame):
        if self.adhesion:
            print "Hi"
            self.adhesion.Update_Springs(True)

    def Postprocess_Frame(self,frame):
        if self.adhesion:
            self.adhesion.Write_State(physbam.Get_Stream_Type("float"),self.output_directory+"/adhesion.%d"%frame)


    def Update_Solids_Parameters(self,time):
        if self.adhesion and time>6.5 and self.changed_stiction_parameters==False:
            self.changed_stiction_parameters=True
            minimum_width=.01
            maximum_width=.2
            self.adhesion.Set_Parameters(10.,10.,minimum_width,maximum_width,50) # ym,overdamp,start dist, stop dist
            self.adhesion.Set_Restlength(.2)

    def Preprocess_Solids_Substep(self,time,substep):
        pass

    # Unused callbacks
    def Update_Time_Varying_Material_Properties(self,time):        pass
    def Update_Collision_Body_Positions_And_Velocities(self,time): pass
    def Postprocess_Solids_Substep(self,time,substep): pass
    def Apply_Constraints(self,dt,time): pass
    def Add_External_Forces(self,F,time): pass
    def Add_External_Forces(self,wrench,time): pass
    def Limit_Solids_Dt(self,time): pass
    def Update_Fragments(self): pass
    def Postprocess_Fragments(self): pass
    def Postprocess_Super_Fragments(self,swap_pairs,rebuild,old_max): pass
    def Update_Super_Fragments(self,swap_pairs,rebuild,old_max): pass
    def Add_External_Super_Fragment_Connectivity(self,union_find): pass

physbam.Initialize_Logging(False,False,9999,False)
physbam.Set_Floating_Point_Exception(True,True,True,True,False,False)
#physbam.LOG.Initialize_Logging(False,False,9999,False)
example=HAIR_EXAMPLE(physbam.stream_type_float,0,physbam.FLUIDS_PARAMETERS_f3.NONE)
#del example
driver=physbam.SOLIDS_FLUIDS_DRIVER_UNIFORM_f3(example)
driver.Execute_Main_Program()

