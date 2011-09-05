#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIVER
//#####################################################################
#ifndef __RIVER__
#define __RIVER__

#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Geometry/Implicit_Objects_RLE/RLE_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_RLE/READ_WRITE_RLE_IMPLICIT_OBJECT.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_RLE.h>
namespace PhysBAM{

template<class T_input>
class RIVER:public SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_3D<T_input> >
{
    typedef T_input T;
    typedef RLE_GRID_3D<T> T_GRID;
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::HORIZONTAL_GRID::CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID> BASE;
    using BASE::output_directory;using BASE::write_frame_title;using BASE::particle_levelset;using BASE::incompressible;using BASE::fluids_parameters;
    using BASE::data_directory;using BASE::Adjust_Phi_With_Sources;using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::solids_parameters;
    using BASE::refine_all_water;using BASE::frame_rate;using BASE::first_frame;using BASE::last_frame;using BASE::initial_time;using BASE::use_incompressible_cfl;
    using BASE::clamp_long_velocities_to_short_velocities;

    T_GRID& grid;
    GRID<VECTOR<T,2> > river_height_grid;
    ARRAY<T,VECTOR<int,2> > river_height;
    LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> > ,T> river_height_interpolation;
    T_GRID river_grid;
    ARRAY<T> river_phi;
    RLE_IMPLICIT_OBJECT<TV> river_implicit_surface;

    T initial_height,initial_slope;
    PLANE<T> initial_plane;
    ARRAY<BOX<TV> > sources;
    ARRAY<MATRIX<T,4> > world_to_source;
    ARRAY<TV> source_velocity;

    RIVER(STREAM_TYPE stream_type,const int resolution,const T optical_depth)
        :SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_3D<T> >(stream_type),grid(*fluids_parameters.grid),river_implicit_surface(river_grid,river_phi),
        initial_height((T).3),initial_slope(-(T).2),initial_plane(TV(-initial_slope,1,0).Normalized(),TV(0,initial_height,0))
    {
        T_GRID& grid=*fluids_parameters.grid;
        LOG::cout<<"running river at resolution "<<resolution<<", optical depth "<<optical_depth<<std::endl;

        // frames
        frame_rate=96;
        first_frame=0;
        last_frame=1000;

        // set up the standard fluid environment
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_debug_data=true;
        fluids_parameters.write_ghost_values=true;fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.incompressible_iterations=40;
        fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.second_order_cut_cell_method=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_vorticity_confinement_fuel=false;

        // override standard fluid environment
        fluids_parameters.incompressible_iterations=200;
        fluids_parameters.solve_neumann_regions=false;
        clamp_long_velocities_to_short_velocities=true;

        // read river data
        std::string river_name="river3_extended";
        //std::string river_name="river2_extended";
        FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/RLE/"+river_name+"_grid.gz",river_height_grid);
        FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/RLE/"+river_name+"_height.gz",river_height);
        LOG::cout<<"river height grid = "<<river_height_grid<<std::endl;

        // setup grid
        //grid.Set_Uniform_Grid(GRID<TV>(50*resolution+1,30*resolution+1,10*resolution+1,BOX<TV>(0,5,-1,2,0,1)));
        grid.Set_Uniform_Grid(GRID<TV>(100*resolution+1,30*resolution+1,10*resolution+1,BOX<TV>(0,10,-1,2,0,1)));
        assert(river_height_grid.Domain()==grid.horizontal_grid.Domain());
        LOG::cout<<"uniform_grid = "<<grid.uniform_grid<<std::endl;
        
        // set bandwidths
        grid.Set_Positive_Bandwidth_In_Cells(3);
        if(optical_depth<5*grid.Minimum_Edge_Length()){
            LOG::cout<<"Warning: optical depth too small, switching to bandwidth of 5 cells."<<std::endl;
            grid.Set_Negative_Bandwidth_In_Cells(5);}
        else grid.Set_Negative_Bandwidth_From_Optical_Depth(optical_depth);
        grid.Set_Linear_Pressure_And_Linear_Horizontal_Velocity_Profile();

        // set output directory
        output_directory=STRING_UTILITIES::string_sprintf("River/%s_resolution_%d_%d_%d_depth_%g",river_name.c_str(),(grid.uniform_grid.counts.x-1),(grid.uniform_grid.counts.y-1),(grid.uniform_grid.counts.z-1),optical_depth);
        if(refine_all_water) output_directory+="_uniform";

        // add river object
        LOG::Time("loading river surface");
        //solid_body_collection.rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/RLE/"+river_name+"_rle_levelset",1,true,false,false,false);
        solid_body_collection.rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/RLE/"+river_name,1,true,false,false,false);
        LOG::Time("loading river levelset");
        FILE_UTILITIES::Read_From_File<float>(data_directory+"/RLE/"+river_name+"_rle_grid_and_levelset",river_implicit_surface);
        LOG::Stop_Time();
        solid_body_collection.rigid_body_collection.Rigid_Body(1).Add_Structure(river_implicit_surface);
        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(solid_body_collection.rigid_body_collection.rigid_geometry_collection);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Get_Collision_Geometry(1)->refine_nearby_fluid=false;

        // nonstandard settings
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=false;
        fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;
        particle_levelset.levelset.clamp_phi_with_collision_bodies=false;
        particle_levelset.use_removed_negative_particles_in_long_cells=true;
        fluids_parameters.collidable_phi_replacement_value=grid.positive_bandwidth;
        use_incompressible_cfl=false;
        incompressible.donor_cell_advection.clamp_divergence_fix=true;

        // setup sources
        sources.Append(BOX<TV>(-1,(T).1,-1,(T).35,-1,2));
        world_to_source.Append(MATRIX<T,4>::Identity_Matrix());
        source_velocity.Append(TV((T)2.5,0,0));

        // set coupling parameters
        fluids_parameters.Use_Fluid_Coupling_Defaults();
    }

    ~RIVER()
    {}

    // Unused callbacks
    void Initialize_Velocities() PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Grids_Extra() PHYSBAM_OVERRIDE {}

    T Initial_Ground(const TV_HORIZONTAL& X) const PHYSBAM_OVERRIDE
    {T height=river_height_interpolation.Clamped_To_Array(river_height_grid,river_height,X);
    T dy=grid.uniform_grid.dX.y,tolerance=(T).01*dy;
    //while(river_implicit_surface(TV(X.x,height+dy,X.y))<=tolerance)height+=dy;
    if(river_implicit_surface(TV(X.x,height,X.y))>=-tolerance) height-=dy;
    return height;}

    T Initial_Phi(const TV& X) const PHYSBAM_OVERRIDE
    {T phi=initial_plane.Signed_Distance(X);
    for(int s=1;s<=sources.m;s++) phi=min(phi,sources(s).Signed_Distance(world_to_source(s).Homogeneous_Times(X)));
    for(int s=1;s<=sources.m;s++) phi=max(phi,sources(s).Signed_Distance(world_to_source(s).Homogeneous_Times(X))); // initial water only from source
    phi=max(phi,-river_implicit_surface(X));
    phi=max(phi,river_height_interpolation.Clamped_To_Array(river_height_grid,river_height,X.Horizontal_Vector())-X.y-grid.uniform_grid.dX.y);
    return phi;}

    T Initial_Phi_Object(const TV& X) const PHYSBAM_OVERRIDE
    {return 1;}

    void Construct_Levelsets_For_Objects(const T time)
    {/*tests->Update_Rigid_Bodies(time);*/SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Construct_Levelsets_For_Objects(time);}

    void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Update_Fluid_Parameters(dt,time);/*tests->Update_Sources(time);*/}

    bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
    {if(time==initial_time) return false; // skip initial call
    for(int s=1;s<=sources.m;s++)Adjust_Phi_With_Source(sources(s),world_to_source(s),true);return false;}

    void Get_Source_Reseed_Mask(ARRAY<bool>*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
    {for(int s=1;s<=sources.m;s++)Get_Source_Reseed_Mask(sources(s),world_to_source(s),cell_centered_mask,s==1);}

    void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
    {for(int s=1;s<=sources.m;s++)Get_Source_Velocities(sources(s),world_to_source(s),source_velocity(s));}

    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
    {/*tests->Limit_Dt(dt,time);*/}

    void Get_Cell_Should_Be_Long(ARRAY<bool>& cell_should_be_long,const T time) const PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Get_Cell_Should_Be_Long(cell_should_be_long,time);
    for(int s=1;s<=sources.m;s++)Get_Source_Cell_Should_Be_Long(sources(s),world_to_source(s),cell_should_be_long,true);
    /*BOX<TV> refine_box(4,5.4,-1.02,4,-5,5);
    const ARRAY<T>& phi=particle_levelset.phi;
    for(CELL_ITERATOR cell(grid,0);cell;cell++)if(cell.Short() && refine_box.Lazy_Inside(cell.X())){int c=cell.Cell();
        if(phi(c)<0) cell_should_be_long(c)=false;}*/}

    void Transfer_Extra_State(const T_GRID& new_grid) PHYSBAM_OVERRIDE
    {SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::Transfer_Extra_State(new_grid);
    // ensure that we don't get water beneath the ground
    ARRAY<T>& phi=particle_levelset.phi;
    for(HORIZONTAL_CELL_ITERATOR iterator(new_grid.horizontal_grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
        int c=new_grid.columns(iterator.Cell_Index())(1).cell;
        phi(c)=phi(c+1)=grid.positive_bandwidth;}
    // delete removed particles that fall beneath the ground
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*particle_levelset.removed_negative_particles_in_long_cells;
    for(int k=particles.array_collection->Size();k>=1;k--)if(particles.X(k).y<river_height_interpolation.Clamped_To_Array(river_height_grid,river_height,particles.X(k).Horizontal_Vector()))
        particles.array_collection->Delete_Element(k);}

//#####################################################################
};
}
#endif
#endif
