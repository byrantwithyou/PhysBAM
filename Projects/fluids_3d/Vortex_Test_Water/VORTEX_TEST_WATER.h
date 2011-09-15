//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VORTEX_TEST_WATER__
#define __VORTEX_TEST_WATER__

#include <PhysBAM_Tools/Particles_Interpolation/SCATTERED_INTERPOLATION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template <class T,class RW=T>
class VORTEX_TEST_WATER:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::fluids_parameters;using BASE::first_frame;using BASE::last_frame;
    using BASE::frame_rate;using BASE::write_output_files;using BASE::output_directory;
    T source_vorticity_magnitude;
    T particle_vorticity_minimum_density;
    RANDOM_NR3 random;
    VORTEX_PARTICLE_EVOLUTION_3D<T,RW> vortex_particle_evolution;
    T rho,rho_bottom,rho_top,buoyancy_constant;
    BOX_3D<T> source_domain;
    
    VORTEX_TEST_WATER()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(1,fluids_parameters.WATER)
    {
        fluids_parameters.grid.Initialize(50,50,50,0,1,0,1,0,1);
        first_frame=0;
        last_frame=200;
        frame_rate=24;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=false;
        fluids_parameters.use_vorticity_confinement=false;fluids_parameters.confinement_parameter=(T).3;
        fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.kolmogorov=(T)0;
        write_output_files=true;fluids_parameters.write_debug_data=true;
        output_directory="Vortex_Test_Water/output";
        // Setup vortex particle parameters
        fluids_parameters.use_vorticity_confinement=false; // use our particle instead... 
        fluids_parameters.use_body_force=true;
        fluids_parameters.confinement_parameter=(T)0.3;
        vortex_particle_evolution.particle_confinement_parameter=(T)1.00;
        vortex_particle_evolution.renormalize_vorticity_after_stretching_tilting=true;
        particle_vorticity_minimum_density=(T).9;
        source_vorticity_magnitude=(T)10000;
        vortex_particle_evolution.Set_Radius((T).2);
        vortex_particle_evolution.Initialize(fluids_parameters.grid);
        vortex_particle_evolution.force_scaling=1e-2;
        VORTICITY_PARTICLES<T,VECTOR<T,3> >& particles=vortex_particle_evolution.vorticity_particles;
        // add two test particles
        //int index=particles.array_collection->Add_Element();
        //particles.X(index)=VECTOR<T,3>((T).75,(T).4,(T).5);
        //particles.vorticity(index)=VECTOR<T,3>(.5,1,-.5);
        //index=particles.array_collection->Add_Element();
        //particles.X(index)=VECTOR<T,3>((T).25,(T).4,(T).5);
        //particles.vorticity(index)=VECTOR<T,3>(-.5,-1,-.5);
        int index=particles.array_collection->Add_Element();
        particles.X(index)=VECTOR<T,3>((T).5,(T).45,(T).5);
        particles.vorticity(index)=VECTOR<T,3>(0,1,0);
    }
    
    ~VORTEX_TEST_WATER()
    {}
    
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=iterator.Location().y-.5;
}

//#####################################################################
// Function Update_Forces
//#####################################################################
void Get_Body_Force(ARRAY<T,FACE_INDEX<3> >& force,const T dt,const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,FACE_INDEX<3> > face_velocities_ghost(grid,3,false);
    fluids_parameters.fluid_boundary->Fill_Ghost_Cells_Face(grid,fluids_parameters.incompressible->projection.face_velocities,face_velocities_ghost,time);
    
    // compute force
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > cell_centered_force(grid,3);
    vortex_particle_evolution.Compute_Body_Force(face_velocities_ghost,cell_centered_force,dt,time);
    for(UNIFORM_GRID_ITERATOR_FACE_3D<T> iterator(grid,1);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        force.Component(axis)(iterator.Face_Index())=(T).5*(cell_centered_force(iterator.First_Cell_Index())[axis]+cell_centered_force(iterator.Second_Cell_Index())[axis]);}

    // advance Particles
    vortex_particle_evolution.Euler_Step(face_velocities_ghost,dt,time);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Write_Output_Files(frame);
    //FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/fluid_force.%d",output_directory.c_str(),frame),fluids_parameters.incompressible.force);
    vortex_particle_evolution.Write_Output_Files(output_directory,frame);
}
//#####################################################################
};      
}
#endif


