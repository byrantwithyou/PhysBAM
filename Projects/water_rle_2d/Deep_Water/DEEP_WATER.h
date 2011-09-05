#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEEP_WATER
//#####################################################################
#ifndef __DEEP_WATER__
#define __DEEP_WATER__

#include <PhysBAM_Tools/Read_Write/Grids_RLE/READ_WRITE_RLE_GRID.h>
#include <PhysBAM_Geometry/Fourier_Transforms_Calculations/DEEP_WATER_EVOLUTION_HEIGHTS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_RLE.h>
#include "DEEP_WATER_TESTS_2D.h"
#include "DEEP_WATER_TESTS_3D.h"
namespace PhysBAM{

template<class T> class DEEP_WATER_TESTS_3D;

template<class T_GRID>
class DEEP_WATER:public SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::HORIZONTAL_GRID T_HORIZONTAL_GRID;typedef typename GRID_ARRAYS_POLICY<T_HORIZONTAL_GRID>::ARRAYS_SCALAR T_ARRAYS_HORIZONTAL_T;
    typedef typename T_GRID::FACE_Y_ITERATOR FACE_Y_ITERATOR;typedef typename T_HORIZONTAL_GRID::CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;
    typedef typename TV_HORIZONTAL::template REBIND<int>::TYPE TV_HORIZONTAL_INT;

    typedef typename IF<TV::m==2,DEEP_WATER_TESTS_2D<T>,DEEP_WATER_TESTS_3D<T> >::TYPE T_STANDARD_TESTS;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID> BASE;
    using BASE::fluids_parameters;using BASE::output_directory;using BASE::write_frame_title;using BASE::particle_levelset;using BASE::incompressible;
    using BASE::data_directory;using BASE::Adjust_Phi_With_Sources;using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::solids_parameters;using BASE::refine_all_water;
    using BASE::vertical_refinement_depth;using BASE::use_phi_for_vertical_refinement;using BASE::enforce_refinement_slope;using BASE::frame_rate;using BASE::last_frame;
    using BASE::use_incompressible_cfl;using BASE::clamp_long_velocities_to_short_velocities;using BASE::ground_j;using BASE::deep_water;using BASE::use_deep_water;

    T_STANDARD_TESTS* tests; 

    DEEP_WATER(const STREAM_TYPE stream_type,const T optical_depth)
        :SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>(stream_type)
    {
        T_GRID& grid=*fluids_parameters.grid;
        tests=new T_STANDARD_TESTS(*this,fluids_parameters,solid_body_collection.rigid_body_collection,test_number,resolution);

        particle_levelset.levelset.clamp_phi_with_collision_bodies=false;
        particle_levelset.use_removed_negative_particles_in_long_cells=true;
        clamp_long_velocities_to_short_velocities=false;

        if(!optical_depth) refine_all_water=true;

        if(refine_all_water) output_directory+="_uniform";
        else output_directory+=STRING_UTILITIES::string_sprintf("_Depth_%g",optical_depth);

        LOG::cout<<"uniform grid = "<<tests->grid<<std::endl;
        LOG::cout<<"optical depth = "<<optical_depth<<std::endl;
        grid.Set_Uniform_Grid(tests->grid);
        grid.Set_Positive_Bandwidth_In_Cells(3);
        if(optical_depth<5*grid.Minimum_Edge_Length()){
            if(!refine_all_water) PHYSBAM_FATAL_ERROR();
            grid.Set_Negative_Bandwidth_In_Cells(5);}
        else grid.Set_Negative_Bandwidth_From_Optical_Depth(optical_depth);
        grid.Set_Linear_Pressure_And_Linear_Horizontal_Velocity_Profile();
        LOG::cout<<"negative bandwidth = "<<grid.negative_bandwidth<<", positive_bandwidth = "<<grid.positive_bandwidth<<std::endl;

        if(T_GRID::dimension==3) fluids_parameters.number_particles_per_cell=32;
        use_incompressible_cfl=true;
        fluids_parameters.incompressible_iterations=200;

        fluids_parameters.collidable_phi_replacement_value=grid.positive_bandwidth;
        tests->Initialize_Bodies();
        tests->Initialize_Advection();

        // deep water setup
        fluids_parameters.periodic.Fill(true);fluids_parameters.periodic.y=false;
        use_deep_water=tests->use_deep_water;
        deep_water.phillips_spectrum.amplitude=0;

        //frame_rate*=5;
        //last_frame*=5;
    }

    ~DEEP_WATER()
    {delete tests;}

    // Unused callbacks
    void Initialize_Velocities() PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Grids_Extra() PHYSBAM_OVERRIDE {}
    void Modify_Grid_After_Rebuild(T_GRID& new_grid,const T time) PHYSBAM_OVERRIDE {}
    void Transfer_Extra_State(const T_GRID& new_grid) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

    T Initial_Ground(const TV_HORIZONTAL& X) const PHYSBAM_OVERRIDE
    {return tests->deep_water_level;}

    T Initial_Phi(const TV& X) const PHYSBAM_OVERRIDE
    {return max(tests->Initial_Phi(X),-tests->Initial_Phi_Object(X));}

    T Initial_Phi_Object(const TV& X) const PHYSBAM_OVERRIDE
    {return tests->Initial_Phi_Object(X);}

    void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Update_Fluid_Parameters(dt,time);tests->Update_Sources(time);}

    void Set_Kinematic_Positions(TV& X, ROTATION<TV>& rotation,const T time,const int id) PHYSBAM_OVERRIDE
    {tests->Set_Kinematic_Positions(X,rotation,time,id);}

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
    {return tests->Set_Kinematic_Velocities(twist,time,id);}

    bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
    {for(int s=1;s<=tests->sources.m;s++)Adjust_Phi_With_Source(tests->sources(s),tests->world_to_source(s));return false;}

    void Get_Source_Reseed_Mask(ARRAY<bool>*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
    {for(int s=1;s<=tests->sources.m;s++)Get_Source_Reseed_Mask(tests->sources(s),tests->world_to_source(s),cell_centered_mask,s==1);}

    void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
    {for(int s=1;s<=tests->sources.m;s++)Get_Source_Velocities(tests->sources(s),tests->world_to_source(s),tests->source_velocity(s));}

    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
    {tests->Limit_Dt(dt,time);}

    void Get_Cell_Should_Be_Long(ARRAY<bool>& cell_should_be_long,const T time) const PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Get_Cell_Should_Be_Long(cell_should_be_long,time);
    for(int s=1;s<=tests->sources.m;s++)Get_Source_Cell_Should_Be_Long(tests->sources(s),tests->world_to_source(s),cell_should_be_long);}

    void Get_Neumann_And_Dirichlet_Boundary_Conditions(const T dt,const T time) PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Get_Neumann_And_Dirichlet_Boundary_Conditions(dt,time);
    LOG::Time("getting top surface pressure");
    tests->Get_Surface_Pressure(dt,time);}

//#####################################################################
};
}
#endif
#endif
