//#####################################################################
// Copyright 2004-2006, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPINNING_BAR
//#####################################################################
#ifndef __SPINNING_BAR__
#define __SPINNING_BAR__

#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SPINNING_BAR:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::viewer_dir;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::fluid_collection;using BASE::resolution;
    using BASE::user_last_frame;

    T angular_velocity;

    SPINNING_BAR(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,1,fluids_parameters.WATER)
    {
        parse_args.Parse();
        fluids_parameters.grid->Initialize(TV_INT(resolution*20,resolution*20),RANGE<TV>(TV(0,0),TV(1,1)));
        if(!user_last_frame) last_frame=2000;
        if(!this->user_frame_rate) frame_rate=24;
        restart=0;
        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Spinning_Bar/output_%d",resolution);
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false; 
        fluids_parameters.number_particles_per_cell=32;
        fluids_parameters.particle_half_bandwidth=1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fluids_parameters.reseeding_frame_rate=1;
        fluids_parameters.bias_towards_negative_particles=true;
        fluids_parameters.incompressible_iterations=20;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=true;fluids_parameters.write_ghost_values=true;
        fluids_parameters.delete_fluid_inside_objects=true;

        //fluids_parameters.levelset_substeps=.9;

        fluids_parameters.use_strain=true;fluids_parameters.write_strain=true;
        fluids_parameters.elastic_modulus=1000;
        fluids_parameters.plasticity_alpha=0;
        fluids_parameters.cfl/=10;
        fluids_parameters.gravity=TV();
        frame_rate*=5;

        angular_velocity=5;
        //use_external_velocity=true;

        fluids_parameters.viscosity=0;//(T)25;
        if(fluids_parameters.viscosity)
            if(!this->user_output_directory)
                viewer_dir.output_directory+=LOG::sprintf("_v%g",fluids_parameters.viscosity);
    }

//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    override
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    RANGE<TV> box(TV((T).25,(T).4),TV((T).75,(T).6));

    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++)
        phi(i,j)=box.Signed_Distance(grid.X(TV_INT(i,j)));
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() override
{
    GRID<TV>& grid=*fluids_parameters.grid;
    TV center=grid.Domain().Center();
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(iterator.Face_Index())=
            TV::Cross_Product(VECTOR<T,1>(angular_velocity),iterator.Location()-center)[axis];}
    //ARRAY<SYMMETRIC_MATRIX<T,2> ,VECTOR<int,2> >::copy(SYMMETRIC_MATRIX<T,2>(1,0,-1),e);
}
//#####################################################################
// Function Get_External_Velocity
//#####################################################################
void Get_External_Velocity(ARRAY<TV,VECTOR<int,2> >& V_blend,ARRAY<T,VECTOR<int,2> >& blend,const T time) override
{
    //TV center=grid.Domain().Center();
    //for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)
    //    V(i,j)=angular_velocity*(grid.X(TV_INT(i,j))-center).Rotate_Counterclockwise_90();
} 
//#####################################################################
};      
}
#endif
