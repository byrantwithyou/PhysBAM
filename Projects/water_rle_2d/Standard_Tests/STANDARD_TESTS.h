#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_RLE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_RLE.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_GRID>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename IF<TV::m==2,WATER_STANDARD_TESTS_2D<T_GRID>,WATER_STANDARD_TESTS_3D<T_GRID> >::TYPE T_STANDARD_TESTS;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID> BASE;
    using BASE::output_directory;using BASE::write_frame_title;using BASE::particle_levelset;using BASE::incompressible;using BASE::data_directory;
    using BASE::Adjust_Phi_With_Sources;using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::solids_parameters;using BASE::refine_all_water;
    using BASE::vertical_refinement_depth;using BASE::use_phi_for_vertical_refinement;using BASE::enforce_refinement_slope;using BASE::frame_rate;using BASE::last_frame;
    using BASE::use_incompressible_cfl;using BASE::fluids_parameters;using BASE::clamp_long_velocities_to_short_velocities;

    T_STANDARD_TESTS* tests;
    T ground_level; // 0 for normal standard tests

    STANDARD_TESTS(const STREAM_TYPE stream_type,const T optical_depth,const T vertical_refinement_depth_input,
        const bool enforce_refinement_slope_input)
        :SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>(stream_type),ground_level(0)
    {
        T_GRID& grid=*fluids_parameters.grid;
        tests=new T_STANDARD_TESTS(*this,fluids_parameters,solid_body_collection.rigid_body_collection,test_number,resolution);

        vertical_refinement_depth=vertical_refinement_depth_input;
        use_phi_for_vertical_refinement=false;
        particle_levelset.levelset.clamp_phi_with_collision_bodies=false;
        enforce_refinement_slope=enforce_refinement_slope_input;
        particle_levelset.use_removed_negative_particles_in_long_cells=true;
        clamp_long_velocities_to_short_velocities=true;

        output_directory+=STRING_UTILITIES::string_sprintf("_Depth_%g",optical_depth);
        if(vertical_refinement_depth) output_directory+=STRING_UTILITIES::string_sprintf("_Vertical_%g",vertical_refinement_depth);
        if(enforce_refinement_slope) output_directory+="_slope";
        if(refine_all_water) output_directory+="_uniform";

        LOG::cout<<"uniform grid = "<<tests->grid<<std::endl;
        grid.Set_Uniform_Grid(tests->grid);
        grid.Set_Positive_Bandwidth_In_Cells(3);
        if(optical_depth<5*grid.Minimum_Edge_Length()){
            LOG::cerr<<"Warning: optical depth too small, switching to bandwidth of 5 cells."<<std::endl;
            grid.Set_Negative_Bandwidth_In_Cells(5);}
        else grid.Set_Negative_Bandwidth_From_Optical_Depth(optical_depth);
        grid.Set_Linear_Pressure_And_Linear_Horizontal_Velocity_Profile();
        LOG::cout<<"negative bandwidth = "<<grid.negative_bandwidth<<", positive_bandwidth = "<<grid.positive_bandwidth<<std::endl;

        std::string suffix=T_GRID::dimension==2?"_2D":"";
        tests->rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies"+suffix+"/ground",1,true,true,false);
        PHYSBAM_FATAL_ERROR("refine nearby fluid has moved to the wrapper");
        //tests->rigid_body_collection.Rigid_Body(1).refine_nearby_fluid=false;
        tests->rigid_body_collection.rigid_body_particle.X(1).y=ground_level;
        tests->Initialize_Bodies();

        if(T_GRID::dimension==3) fluids_parameters.number_particles_per_cell=32;
        use_incompressible_cfl=true;
        fluids_parameters.incompressible_iterations=200;

        fluids_parameters.collidable_phi_replacement_value=grid.positive_bandwidth;
        fluids_parameters.solid_affects_fluid=fluids_parameters.fluid_affects_solid=true;
        tests->Initialize_Advection();

        if(test_number==4){
            frame_rate=96;
            last_frame=800;}
    }

    ~STANDARD_TESTS()
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

    T Initial_Ground(const TV_HORIZONTAL& X) const PHYSBAM_OVERRIDE
    {return ground_level;}

    T Initial_Phi(const TV& X) const PHYSBAM_OVERRIDE
    {return max(tests->Initial_Phi(X),-X.y+ground_level);}

    T Initial_Phi_Object(const TV& X) const PHYSBAM_OVERRIDE
    {return min(X.y-ground_level,tests->Initial_Phi_Object(X));}

    void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Update_Fluid_Parameters(dt,time);tests->Update_Sources(time);}

    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
    {tests->Set_Kinematic_Positions(frame,time,id);}

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

//#####################################################################
};
}
#endif
#endif
