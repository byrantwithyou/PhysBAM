//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VORTEX_PARTICLES_STREAM__
#define __VORTEX_PARTICLES_STREAM__

#include <PhysBAM_Tools/Particles_Interpolation/SCATTERED_INTERPOLATION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "RIVER_PHI_BOUNDARY.h"

namespace PhysBAM{

template <class T,class RW=T>
class VORTEX_PARTICLES_STREAM:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::fluids_parameters;using BASE::first_frame;using BASE::last_frame;
    using BASE::frame_rate;using BASE::write_output_files;using BASE::output_directory;

    VORTEX_PARTICLE_EVOLUTION_3D<T,RW> vortex_particle_evolution;
    T vorticity_magnitude;
    T particle_distance_from_surface;
    int vortex_count;
    T inflow;
    RANDOM_NR3 random;    
    RIVER_PHI_BOUNDARY<T,T> custom_phi_boundary;
    bool has_initialized_vortex_particle_evolution;

    VORTEX_PARTICLES_STREAM()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(1,fluids_parameters.WATER),has_initialized_vortex_particle_evolution(0)
    {
        int scale=2;
        fluids_parameters.grid.Initialize(64*scale+1,8*scale+1,16*scale+1,0,1.6,0,.2,0,.4);
        first_frame=0;
        last_frame=200;
        frame_rate=60;
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=false;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.use_vorticity_confinement=false;fluids_parameters.confinement_parameter=(T).3;
        fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.write_particles=true;
        fluids_parameters.kolmogorov=(T)0;
        write_output_files=true;fluids_parameters.write_debug_data=true;
        output_directory="Vortex_Particles_Stream/output";
        // Setup vortex particle parameters
        fluids_parameters.use_body_force=true;
        fluids_parameters.confinement_parameter=(T)0;
        vorticity_magnitude=400; //1000; 
        particle_distance_from_surface=.2;
        BOX_3D<T> box=fluids_parameters.grid.Domain();
        //vortex_particle_evolution.Set_Radius(0.05);
        vortex_particle_evolution.force_scaling=3;
        vortex_particle_evolution.renormalize_vorticity_after_stretching_tilting=true;
        inflow=1;
        custom_phi_boundary.inflow_height=(T).075;
        fluids_parameters.phi_boundary=&custom_phi_boundary;
    }
    
    ~VORTEX_PARTICLES_STREAM()
    {}
    
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
    ARRAY<bool,FACE_INDEX<3> >& psi_N=fluids_parameters.incompressible->projection.elliptic_solver->psi_N;
    ARRAY<bool,VECTOR<int,3> > psi_D=fluids_parameters.incompressible->projection.elliptic_solver->psi_D;
    T rho=1;
    T inflow_height=custom_phi_boundary.inflow_height;

    for(FACE_ITERATOR iterator(fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,1);iterator.Valid();iterator.Next()){
            psi_N.Component(iterator.Axis())(iterator.Face_Index()+VECTOR<int,3>(1,0,0))=true;
            fluids_parameters.incompressible->projection.face_velocities.Component(iterator.Axis())(iterator.Face_Index()+VECTOR<int,3>(1,0,0))=inflow;
            psi_N.Component(2)(iterator.Face_Index())=true;
            fluids_parameters.incompressible->projection.face_velocities.Component(2)(iterator.Face_Index())=0;}
    for(CELL_ITERATOR iterator(fluids_parameters.grid,1,GRID<TV>::GHOST_REGION,1);iterator.Valid();iterator.Next()) 
        psi_D(iterator.Cell_Index())=false;
    for(CELL_ITERATOR iterator(fluids_parameters.grid,1,GRID<TV>::GHOST_REGION,2);iterator.Valid();iterator.Next()){
        psi_D(iterator.Cell_Index())=true;fluids_parameters.incompressible->projection.p(iterator.Cell_Index())=inflow_height-iterator.Location().y*fluids_parameters.gravity*rho;}        
}
//#####################################################################
// Function Update_Force
//#####################################################################
void Get_Body_Force(ARRAY<T,FACE_INDEX<3> >& force,const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(!has_initialized_vortex_particle_evolution){
        vortex_particle_evolution.Initialize(fluids_parameters.grid);has_initialized_vortex_particle_evolution=true;}

    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,FACE_INDEX<3> > face_velocities_ghost(grid,3,false);
    fluids_parameters.fluid_boundary->Fill_Ghost_Cells_Face(grid,fluids_parameters.incompressible->projection.face_velocities,face_velocities_ghost,time);

    // Compute force
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > cell_centered_force(grid,3);
    vortex_particle_evolution.Compute_Body_Force(face_velocities_ghost,cell_centered_force,dt,time);
    for(UNIFORM_GRID_ITERATOR_FACE_3D<T> iterator(grid,1);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        force.Component(axis)(iterator.Face_Index())=(T).5*(cell_centered_force(iterator.First_Cell_Index())[axis]+cell_centered_force(iterator.Second_Cell_Index())[axis]);}

    // Advance Particles
    vortex_particle_evolution.Euler_Step(face_velocities_ghost,dt,time);

    // delete particles that are outside the water
    VORTICITY_PARTICLES<T,VECTOR<T,3> >& vorticity_particles=vortex_particle_evolution.vorticity_particles;
    for(int p=vorticity_particles.array_collection->Size();p>=1;p--) if(fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.Phi(vorticity_particles.X(p))>0)
        vorticity_particles.Delete_Particle(p);
    // seed particles with rejection sampling... to get a location
    T min_cell_size=grid.min_dx_dy_dz;
    int seeded=0,number_to_seed=1;
    while(seeded<number_to_seed){
        VECTOR<T,3> seed_location(min_cell_size,random.Get_Uniform_Number((T)grid.ymin,(T)grid.ymax),random.Get_Uniform_Number((T)grid.zmin,(T)grid.zmax));
        if(fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.Phi(seed_location)<0){
            int particle_id=vorticity_particles.array_collection->Add_Element(); 
            vorticity_particles.X(particle_id)=seed_location;
            vorticity_particles.radius(particle_id)=random.Get_Uniform_Number(.5,1.5)*.05;
            int sign=random.Get_Uniform_Integer(0,1);if(sign==0)sign=-1;
            vorticity_particles.vorticity(particle_id)=VECTOR<T,3>(0,(T)sign*vorticity_magnitude,0);
            seeded++;
            std::cout<<"ADDING PARTICLES"<<std::endl;}}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(fluids_parameters.grid);iterator.Valid();iterator.Next())if(iterator.Axis()==1&&iterator.Location().y<=custom_phi_boundary.inflow_height){
        fluids_parameters.incompressible->projection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=inflow;}
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    T left_edge=3*grid.min_dx_dy_dz;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(iterator.Location().x<left_edge)
        phi(iterator.Cell_Index())=iterator.Location().y-custom_phi_boundary.inflow_height;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=iterator.Location().y-custom_phi_boundary.inflow_height;
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
// Function Read_Output_Files_Fluids
//#####################################################################
void Read_Output_Files_Fluids(const int frame) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Read_Output_Files_Fluids(frame);
    vortex_particle_evolution.Read_Output_Files(output_directory,frame);
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
{
    if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,3> >(fluids_parameters.grid,1);
    T padding=3*fluids_parameters.grid.DX().Max();
    for(CELL_ITERATOR iterator(fluids_parameters.grid);iterator.Valid();iterator.Next()){
        VECTOR<T,3> location=iterator.Location();
        if(location.y <= custom_phi_boundary.inflow_height+padding && location.x < (T).1) (*cell_centered_mask)(iterator.Cell_Index())=true;}
}
//#####################################################################
};      
}
#endif


