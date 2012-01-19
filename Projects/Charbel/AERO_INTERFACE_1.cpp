//#####################################################################
// Copyright 2007-2008, Jon Gretarsson, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_3D.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/POISSON_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/FLOOD_FILL_MPI.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Particles/COMPRESSIBLE_FLUID_PARTICLES.h>
#include "AERO_INTERFACE_1.h"
#include <Deformable_Objects/DEFORMABLE_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> AERO_INTERFACE_1<T>::
AERO_INTERFACE_1(T_GRID domain,const std::string& output_directory_input,int substeps_input)
    :stream_type((RW())),grid(domain.Get_MAC_Grid_At_Regular_Positions()),global_grid(new T_GRID(grid)),euler(domain.Get_MAC_Grid_At_Regular_Positions()),
    rigid_body(0),collision_body_list(0),inaccurate_union(new FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>(grid)),mpi_grid(0),mpi_charbel(0),mpi_boundary(0),
    solid_extrapolation_bandwidth(0),fluid_extrapolation_bandwidth(0),substeps_level(substeps_input),first_frame(1),output_number(0),output_directory(output_directory_input)
{
    euler.Initialize_Domain(grid);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> AERO_INTERFACE_1<T>::
~AERO_INTERFACE_1()
{
    if(rigid_body) delete rigid_body;if(mpi_grid) delete mpi_grid;if(mpi_charbel) delete mpi_charbel;if(mpi_boundary) delete mpi_boundary;
}
//#####################################################################
// Function Phi
//#####################################################################
template<class T> const T AERO_INTERFACE_1<T>::
Phi(const TV& position) const
{
    PHYSBAM_ASSERT(inaccurate_union.get());return inaccurate_union->levelset(levelset_frame.Inverse_Times(position));
}
//#####################################################################
// Function Gradient
//#####################################################################
template<class T> const VECTOR<T,3> AERO_INTERFACE_1<T>::
Gradient(const TV& position) const
{
    PHYSBAM_ASSERT(inaccurate_union.get());return interpolation.Clamped_To_Array((mpi_grid?mpi_grid->local_grid:grid),gradient,position);
}
//#####################################################################
// Function Inside_Solid
//#####################################################################
template<class T> const bool AERO_INTERFACE_1<T>::
Inside_Solid(const TV& position) const
{
    PHYSBAM_ASSERT(inaccurate_union.get());return inaccurate_union->levelset(position)<(T)0;
}
//#####################################################################
// Function Cell_Should_Be_Populated
//#####################################################################
template<class T> const bool AERO_INTERFACE_1<T>::
Cell_Should_Be_Populated(const TV& position, const T& depth) const
{
    PHYSBAM_ASSERT(inaccurate_union.get());
    T phi=inaccurate_union->levelset(position);
    return phi>=0 && phi<=depth;
}
//#####################################################################
// Function Initialize_MPI
//#####################################################################
template<class T> void AERO_INTERFACE_1<T>::
Initialize_MPI(TETRAHEDRALIZED_VOLUME<T>& tet_volume,ARRAY<int>& local_to_global_map,const int& global_particle_count,T depth, TV* dim_min, TV* dim_max){
    LOG::cout << "Original GRID " << grid << std::endl;
    mpi_grid=new MPI_UNIFORM_GRID<T_GRID>(grid,3,false);
    LOG::cout << "MPI GRID " << grid << std::endl;
    euler.Initialize_Domain(grid); // reinitialize with the local grid
    euler.mpi_grid=mpi_grid;
    mpi_boundary=new BOUNDARY_MPI<T_GRID,typename T_GRID::SCALAR,T_GRID::dimension+2>(mpi_grid,*euler.boundary,3);
    euler.Set_Custom_Boundary(*mpi_boundary);
    if(dim_min) *dim_min=grid.Xmin();
    if(dim_max) *dim_max=grid.Xmax();// Report back dimension information

    TV_INT start_index;TV_INT end_index;
    phi_all_solids.Resize(grid);
    T_ARRAYS_VECTOR gradient_global=gradient;gradient.Resize(grid);
    for(int axis=0;axis<T_GRID::dimension;axis++)
        start_index[axis]=mpi_grid->boundaries(axis)(mpi_grid->coordinates[axis]);

    ARRAY<BOX<TV> > local_grids;
    for(int proc=0;proc<mpi_grid->Number_Of_Processors();proc++){
        TV_INT proc_coordinates=mpi_grid->all_coordinates(proc);
        for(int axis=0;axis<T_GRID::dimension;axis++){
            start_index[axis]=mpi_grid->boundaries(axis)(proc_coordinates[axis]);
            end_index[axis]=mpi_grid->boundaries(axis)(proc_coordinates[axis]+1);}
        local_grids.Append(BOX<TV>(mpi_grid->global_grid.X(start_index),mpi_grid->global_grid.X(end_index)).Thickened(tet_volume.Maximum_Edge_Length(0)));}
    ARRAY<ARRAY<int> > tets_to_send(mpi_grid->Number_Of_Processors());
    tet_volume.Update_Tetrahedron_List();
    tet_volume.Initialize_Hierarchy();
    for(int proc=0;proc<mpi_grid->Number_Of_Processors();proc++) tet_volume.hierarchy->Intersection_List(local_grids(proc),tets_to_send(proc));
    delete tet_volume.hierarchy;tet_volume.hierarchy=0;
    delete tet_volume.tetrahedron_list;tet_volume.tetrahedron_list=0;
    // Exchange
    mpi_charbel=new MPI_CHARBEL<T>();
    mpi_charbel->Setup_AEROF_PhysBAM_Mapping(tet_volume,tets_to_send,local_to_global_map,global_particle_count,grid.Domain());

    // Resize local data types to use the local grid
    inaccurate_union.reset(new FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>(mpi_grid->local_grid));
    
    LOG::Time("Initializing the acceleration");
    mpi_charbel->local_tet_volume->Update_Tetrahedron_List();
    mpi_charbel->local_tet_volume->Initialize_Hierarchy();
}
//#####################################################################
// Function Compute_Level_Set
//#####################################################################
template<class T> void AERO_INTERFACE_1<T>::
Compute_Level_Set(TRIANGULATED_SURFACE<T>& triangulated_surface,COMPRESSIBLE_FLUID_PARTICLES<TV>* particles_aerof)
{
    triangulated_surface.Initialize_Hierarchy();
    triangulated_surface.Update_Bounding_Box();
    triangulated_surface.Update_Triangle_List();
    triangulated_surface.Update_Triangle_List();
    triangulated_surface.mesh.Initialize_Adjacent_Elements();
    triangulated_surface.mesh.Initialize_Incident_Elements();
    rigid_body=RIGID_BODY<TV>::Create_From_Triangulated_Surface(triangulated_surface,rigid_body_particles,(T)1);
    collision_body_list.reset(new FLUID_COLLISION_BODY_LIST_UNIFORM<T_GRID>(mpi_grid?mpi_grid->local_grid:grid));
    collision_body_list->Add_Body(rigid_body);

    if(mpi_grid){
        PHYSBAM_ASSERT(inaccurate_union.get());

        int ghost_cells=7;
        ARRAY<VECTOR<int,3> > initialized_indices;
        T_FACE_ARRAYS_BOOL psi_N(mpi_grid->local_grid,ghost_cells);
        T_ARRAYS_BOOL done(mpi_grid->local_grid,ghost_cells);done.Fill(false);T_ARRAYS_SCALAR phi(mpi_grid->local_grid,ghost_cells);
        { // Temporary variables that only flood-fill cares about
            collision_body_list->Rasterize_Objects();
            T_FACE_ARRAYS_SCALAR face_velocities(mpi_grid->local_grid,ghost_cells); // TODO(jontg): We don't need face velocities for the levelset...
            collision_body_list->Compute_Psi_N_No_Velocity(psi_N,face_velocities);

            T_ARRAYS_INT filled_region_colors(mpi_grid->local_grid);filled_region_colors.Fill(0);
            FLOOD_FILL_3D local_flood_fill;local_flood_fill.Optimize_Fill_For_Single_Cell_Regions(true);
            int number_of_regions=local_flood_fill.Flood_Fill(filled_region_colors,psi_N);
            ARRAY<ARRAY<int> > filled_region_ranks;
            FLOOD_FILL_MPI<T_GRID> flood_fill(*mpi_grid,mpi_grid->local_grid,psi_N,number_of_regions,filled_region_colors,filled_region_ranks,NULL);
            int number_of_global_regions=flood_fill.Synchronize_Colors();
            LOG::cout<<"Colored regions for local flood-fill: "<<number_of_regions<<std::endl;
            LOG::cout<<"Colored regions for global flood-fill: "<<number_of_global_regions<<std::endl;

            ARRAY<bool,VECTOR<int,1> > have_checked_color(1,number_of_regions);have_checked_color.Fill(false);
            for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(mpi_grid->local_grid);iterator.Valid();iterator.Next()){bool touches_solid=false;TV_INT cell_index=iterator.Cell_Index();
                for(int axis=0;axis<TV::dimension;axis++) touches_solid=touches_solid || psi_N(axis,iterator.First_Face_Index(axis)) || psi_N(axis,iterator.Second_Face_Index(axis));
                if(touches_solid){done(cell_index)=true;phi(cell_index)=triangulated_surface.Calculate_Signed_Distance(iterator.Location());initialized_indices.Append(cell_index);}}
        }

        int solid_cell_count_to_extrapolate=(int)ceil(solid_extrapolation_bandwidth/grid.DX().Min());
        int fluid_cell_count_to_extrapolate=(int)ceil(fluid_extrapolation_bandwidth/grid.DX().Min());
        int solid_and_fluid_cell_count_to_extrapolate=solid_cell_count_to_extrapolate+fluid_cell_count_to_extrapolate;
        mpi_charbel->Reduce_Max(solid_cell_count_to_extrapolate);
        mpi_charbel->Reduce_Max(fluid_cell_count_to_extrapolate);
        mpi_charbel->Reduce_Max(solid_and_fluid_cell_count_to_extrapolate);

        T_ARRAYS_BOOL done_ghost(mpi_grid->local_grid,solid_and_fluid_cell_count_to_extrapolate+1);
        // Create a boundary that extrapolates by default, and fills ghost_cells otherwise
        T_ARRAYS_BOOL::Put(done,done_ghost);
        BOUNDARY_UNIFORM<T_GRID,bool> done_boundary;
        VECTOR<bool,2*T_GRID::dimension> valid_wall;
        for(int axis=0;axis<T_GRID::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++)
            valid_wall(2*axis+axis_side-2)=!mpi_grid->Neighbor(axis,axis_side);
        done_boundary.Set_Constant_Extrapolation(valid_wall(1),valid_wall(2),valid_wall(3),valid_wall(4),valid_wall(5),valid_wall(6));
        BOUNDARY_MPI<T_GRID,bool> mpi_done_boundary(mpi_grid,done_boundary,solid_and_fluid_cell_count_to_extrapolate);

        T distance_travelled(0);T distance_to_travel=max(solid_extrapolation_bandwidth,fluid_extrapolation_bandwidth);
        int count=0;
        while(distance_travelled<=distance_to_travel){++count;
            mpi_done_boundary.Fill_Ghost_Cells(grid,done,done_ghost,0,(T)0,solid_and_fluid_cell_count_to_extrapolate);
            initialized_indices.Remove_All();
            LEVELSET_3D<T_GRID>& levelset=inaccurate_union->levelset.levelset;
            levelset.Fast_Marching_Method(done_ghost,0,grid.DX().Min()*ghost_cells);
            distance_travelled+=grid.DX().Min()*ghost_cells;
            if(distance_travelled>distance_to_travel){
                levelset.Compute_Gradient(gradient);
                phi_all_solids=levelset.phi;
                phi_all_solids*=(T)-1;}}
        LOG::cout<<"Performed "<<count<<" MPI boundary syncs to populate the levelset"<<std::endl;
        
        COMPRESSIBLE_FLUID_PARTICLES<TV>* particles;
        mpi_charbel->Exchange_Compressible_Data(*particles_aerof);
        particles=&mpi_charbel->physbam_particles;

        for(int i=0;i<particles->number;i++){
            particles->phi(i)=inaccurate_union->levelset(particles->X(i));
            particles->grad_phi(i)=interpolation.Clamped_To_Array(mpi_grid->local_grid,gradient,particles->X(i));}

        LOG::Time("Exchanging Compressible Data Back");
        mpi_charbel->Exchange_Back_Compressible_Data(*particles_aerof);}
    else{
        LEVELSET_UNIFORM<T_GRID>& levelset=inaccurate_union->levelset.levelset;
        LEVELSET_MAKER_UNIFORM<T> levelset_maker;
        levelset_maker.Verbose_Mode();
        levelset_maker.Compute_Level_Set(triangulated_surface,*global_grid,levelset.phi);
        levelset.Compute_Gradient(gradient);
        phi_all_solids=levelset.phi;
        phi_all_solids*=(T)-1;}

    rigid_body->Add_Structure(inaccurate_union->levelset);  // TODO(jontg): At the moment this is expecting an object-space implicit object and receiving a world-space implicit object
}
//#####################################################################
// Function Closest_Boundary_Point
//#####################################################################
template<class T> int AERO_INTERFACE_1<T>::
Closest_Boundary_Point(const TV& location,TV& closest_bondary_point,T& distance)
{
    COLLISION_GEOMETRY_ID body_id;int simplex_id;
    closest_bondary_point=collision_body_list->Closest_Boundary_Point(location,fabs(2*Phi(location)),distance,body_id,simplex_id);
    return simplex_id;
}
//#####################################################################
// Function Initialize_Acceleration_Structures
//#####################################################################
template<class T> void AERO_INTERFACE_1<T>::
Initialize_Acceleration_Structures(TETRAHEDRALIZED_VOLUME<T>& tet_volume,T depth)
{
    LOG::Time("Initializing the acceleration");
    tet_volume.Update_Tetrahedron_List();
    tet_volume.Initialize_Hierarchy();

    LOG::Time("Updating bounding tets");
    ARRAY<bool> particle_not_inside_solid(tet_volume.particles.array_collection->Size());
    for(int i=0;i<particle_not_inside_solid.m;i++)
        particle_not_inside_solid(i)=!Inside_Solid(tet_volume.particles.X(i));

    T_GRID& current_grid(mpi_grid?mpi_grid->local_grid:grid);
    bounding_tet_id.Resize(current_grid);
    bounding_tet_id.Fill(0);
    T max_fluid_bandwidth=(T)0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(current_grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();TV position=current_grid.X(cell_index);
        if(Cell_Should_Be_Populated(position,depth+grid.DX().Max())){
            ARRAY<int> intersection_list;tet_volume.hierarchy->Intersection_List(position,intersection_list);
            for(int i=0;i<intersection_list.m;i++){
                TETRAHEDRON<T>& current_tet=(*tet_volume.tetrahedron_list)(intersection_list(i));
                VECTOR<int,4>& current_tet_nodes=tet_volume.mesh.elements(intersection_list(i));
                bool tet_has_no_ghost_nodes=true;for(int j=0;j<4;j++) tet_has_no_ghost_nodes&=particle_not_inside_solid(current_tet_nodes(j));
                if(tet_has_no_ghost_nodes && current_tet.Inside(position)){
                    bounding_tet_id(cell_index)=intersection_list(i);break;}}
            if(!bounding_tet_id(cell_index))
                max_fluid_bandwidth=max(max_fluid_bandwidth,abs(phi_all_solids(cell_index)));}}

    fluid_extrapolation_bandwidth=max_fluid_bandwidth+grid.DX().Max();
    solid_extrapolation_bandwidth=depth;
    phi_all_solids+=fluid_extrapolation_bandwidth;

    LOG::Stop_Time();
}
//#####################################################################
// Function Compute_Ghost_Cells
//#####################################################################
template<class T> T Point_From_Barycentric_Coordinates(const VECTOR<T,3>& weights,const T& x1,const T& x2,const T& x3,const T& x4)
{return weights.x*x1+weights.y*x2+weights.z*x3+((T)1-weights.x-weights.y-weights.z)*x4;}

template<class T> void AERO_INTERFACE_1<T>::
Compute_Ghost_Cells(TETRAHEDRALIZED_VOLUME<T>& tet_volume_aerof,COMPRESSIBLE_FLUID_PARTICLES<TV>& particles_aerof)
{
    PHYSBAM_ASSERT(inaccurate_union.get());
    T_ARRAYS_BOOL done(mpi_grid?mpi_grid->local_grid:grid);
    done.Fill(false);

    LOG::Time("Populating the Eulerian Grid");
    euler.fluid_affects_solid=false;
    euler.Initialize_Solid_Fluid_Coupling(0,collision_body_list.get());

    int solid_cell_count_to_extrapolate=(int)ceil(solid_extrapolation_bandwidth/grid.DX().Min());
    int fluid_cell_count_to_extrapolate=(int)ceil(fluid_extrapolation_bandwidth/grid.DX().Min());
    int solid_and_fluid_cell_count_to_extrapolate=solid_cell_count_to_extrapolate+fluid_cell_count_to_extrapolate;
    if(mpi_grid){
      mpi_charbel->Reduce_Max(solid_cell_count_to_extrapolate);
      mpi_charbel->Reduce_Max(fluid_cell_count_to_extrapolate);
      mpi_charbel->Reduce_Max(solid_and_fluid_cell_count_to_extrapolate);}

    TETRAHEDRALIZED_VOLUME<T>* tet_volume;COMPRESSIBLE_FLUID_PARTICLES<TV>* particles;
    if(mpi_grid){
        mpi_charbel->Exchange_Compressible_Data(particles_aerof);
        tet_volume=mpi_charbel->local_tet_volume;
        particles=&mpi_charbel->physbam_particles;}
    else{
        tet_volume=&tet_volume_aerof;
        particles=&particles_aerof;}

    // TODO Copy values over into a local thing, then later interpolate back
    LOG::Time("Copying particle data to local storage");
    T_ARRAYS_STORED_DIMENSION stored_values(grid);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(bounding_tet_id(cell_index)){
            VECTOR<int,4>& corners=tet_volume->mesh.elements(bounding_tet_id(cell_index));
            VECTOR<VECTOR<T,T_GRID::dimension+2>,4>& current_cell_stored_values=stored_values(cell_index);
            for(int i=0;i<4;i++){
                VECTOR<T,T_GRID::dimension+2>& current_node_stored_value=current_cell_stored_values(i);
                current_node_stored_value(1)=particles->rho(corners(i));
                for(int j=0;j<T_GRID::dimension;j++) current_node_stored_value(j+1)=particles->V(corners(i))(j);
                current_node_stored_value(T_GRID::dimension+2)=particles->E(corners(i));}}}

    LOG::Time("Interpolating particle data to the Eulerian grid");
    T min_phi=fluid_extrapolation_bandwidth;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        TV position=grid.X(cell_index);
        if(bounding_tet_id(cell_index)){
            TETRAHEDRON<T>& current_tet=(*tet_volume->tetrahedron_list)(bounding_tet_id(cell_index));
            TV weights=current_tet.Barycentric_Coordinates(position);
            VECTOR<VECTOR<T,T_GRID::dimension+2>,4>& current_tet_values=stored_values(cell_index);
            for(int j=0;j<T_GRID::dimension+2;j++)
                euler.U(cell_index)(j)=Point_From_Barycentric_Coordinates(weights,current_tet_values(1)(j),current_tet_values(2)(j),current_tet_values(3)(j),current_tet_values(4)(j));
            done(cell_index)=true;}
        else if(phi_all_solids(cell_index)<min_phi) done(cell_index)=true;}

    LOG::Time("Extrapolating state into solid");
    
    EXTRAPOLATION_UNIFORM<T_GRID,VECTOR<T,T_GRID::dimension+2> > extrapolation(grid,phi_all_solids,euler.U,
        (mpi_grid?solid_and_fluid_cell_count_to_extrapolate:fluid_cell_count_to_extrapolate));
    extrapolation.Set_Band_Width((T)(solid_and_fluid_cell_count_to_extrapolate));
    
    if(mpi_grid){
      BOUNDARY_MPI<T_GRID,VECTOR<T,T_GRID::dimension+2> > mpi_boundary(mpi_grid,*extrapolation.boundary,solid_and_fluid_cell_count_to_extrapolate);
      extrapolation.Set_Custom_Boundary(mpi_boundary);
      BOUNDARY_MPI<T_GRID,T> mpi_levelset_boundary(mpi_grid,*extrapolation.levelset.boundary,solid_and_fluid_cell_count_to_extrapolate);
      extrapolation.levelset.Set_Custom_Boundary(mpi_levelset_boundary);
      T_ARRAYS_BOOL done_ghost(grid,solid_and_fluid_cell_count_to_extrapolate);
      // Create a boundary that extrapolates by default, and fills ghost_cells otherwise
      T_ARRAYS_BOOL::Put(done,done_ghost);
      BOUNDARY_UNIFORM<T_GRID,bool> done_boundary;
      VECTOR<bool,2*T_GRID::dimension> valid_wall;
      for(int axis=0;axis<T_GRID::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++)
          valid_wall(2*axis+axis_side-2)=!mpi_grid->Neighbor(axis,axis_side);
      done_boundary.Set_Constant_Extrapolation(valid_wall(1),valid_wall(2),valid_wall(3),valid_wall(4),valid_wall(5),valid_wall(6));
      BOUNDARY_MPI<T_GRID,bool> mpi_done_boundary(mpi_grid,done_boundary,solid_and_fluid_cell_count_to_extrapolate);
      mpi_done_boundary.Fill_Ghost_Cells(grid,done,done_ghost,0,(T)0,solid_and_fluid_cell_count_to_extrapolate);
      extrapolation.Set_Custom_Seed_Done(&done_ghost);
      extrapolation.Extrapolate((T)0,true);}
    else{
      extrapolation.Set_Custom_Seed_Done(&done);
      extrapolation.Extrapolate((T)0,false);}

    LOG::Time("Interpolating to ghost nodes");

    T_LINEAR_INTERPOLATION_DIMENSION interpolation_dimension;
    T max_distance,object_velocity_normal_component;TV location,normal_direction,reflected_point;
    BOX<TV> domain=grid.Domain();
    for(int i=0;i<particles->number;i++){
        TV& location=particles->X(i);
        if(domain.Inside(location,(T)0) && Inside_Solid(location) && abs(inaccurate_union->levelset(location)) <= solid_extrapolation_bandwidth){
            T_ARRAYS_ELEMENT values=interpolation_dimension.Clamped_To_Array(grid,euler.U,location);
            max_distance=(T)-2*inaccurate_union->levelset(location);
            euler.Get_Neumann_Data(location,max_distance,normal_direction,object_velocity_normal_component,reflected_point);

            // High order interpolation
            // TODO inaccurate_union->Inside_Any_Simplex should take simplex_id as a ptr that defaults to null...
            int simplex_id;
            BOX<TV> domain=grid.Domain();
            if(euler.use_higher_order_solid_extrapolation && !euler.inaccurate_union->Inside_Any_Simplex(reflected_point,simplex_id) && domain.Lazy_Inside(reflected_point))
                values=interpolation_dimension.Clamped_To_Array(grid,euler.U,reflected_point);

            euler.conservation->object_boundary->Apply_Neumann_Boundary_Condition(values,normal_direction,object_velocity_normal_component);

            particles->rho(i)=values(1);
            TV& velocity=particles->V(i);for(int j=0;j<TV::dimension;j++) velocity(j)=values(j+1);
            particles->E(i)=values(2+T_GRID::dimension);}}

    if(mpi_grid){
      LOG::Time("Exchanging Compressible Data Back");
      mpi_charbel->Exchange_Back_Compressible_Data(particles_aerof);    }

    LOG::Stop_Time();
}
//#####################################################################
// Function Compute_Lift
//#####################################################################
template<class T> VECTOR<T,3> AERO_INTERFACE_1<T>::
Compute_Lift()
{
    PHYSBAM_NOT_IMPLEMENTED();  // Should calculate Jp as is done by the SEARCH_CONTROLLER
    return TV();
}
//#####################################################################
// Function Write_Fluid_Output_Files
//#####################################################################
template<class T> void AERO_INTERFACE_1<T>::
Write_Fluid_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/grid",grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/grid."+f,grid);

    EOS_GAMMA<T> *tmp_eos = dynamic_cast<EOS_GAMMA<T>*>(euler.eos);
    T_ARRAYS_DIMENSION_SCALAR U_ghost(grid,3);
    euler.boundary->Fill_Ghost_Cells(grid,euler.U,U_ghost,(T)0.01,0); // TODO: use real dt/time
    T_ARRAYS_VECTOR velocity(grid,3);
    T_ARRAYS_SCALAR density(grid,3),energy(grid,3),internal_energy(grid,3),pressure(grid,3),
        entropy(grid,3),speedofsound(grid,3),machnumber(grid,3);
    T_ARRAYS_SCALAR density_gradient(grid),pressure_gradient(grid);
    T_ARRAYS_VECTOR velocity_plus_c(grid,3),velocity_minus_c(grid,3),momentum(grid,3);
    T_FACE_ARRAYS_SCALAR face_velocities(grid);
    EULER_PROJECTION_UNIFORM<T_GRID>::Compute_Density_Weighted_Face_Velocities(grid,face_velocities,U_ghost,euler.euler_projection.elliptic_solver->psi_N);
    for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();TV position=grid.X(cell);
        //machnumber(cell)=velocity(cell).Magnitude()/speedofsound(cell);
        if(!Inside_Solid(position)){
            density(cell)=U_ghost(1,cell);
            velocity(cell)=euler.Get_Velocity(U_ghost,cell);
            energy(cell)=U_ghost(cell)(TV::dimension+2);
            internal_energy(cell)=euler.e(U_ghost(cell));

            pressure(cell)=((tmp_eos->gamma-1)*(energy(cell)-((T).5*density(cell)*velocity(cell).Magnitude_Squared())));
            entropy(cell)=tmp_eos->S(density(cell),tmp_eos->e_From_p_And_rho(pressure(cell),density(cell)));
            speedofsound(cell)=tmp_eos->c(density(cell),tmp_eos->e_From_p_And_rho(pressure(cell),density(cell)));

            momentum(cell)=velocity(cell)*density(cell);
            velocity_plus_c(cell)=velocity(cell)+speedofsound(cell)*TV::All_Ones_Vector();
            velocity_minus_c(cell)=velocity(cell)-speedofsound(cell)*TV::All_Ones_Vector();}
        else{
            density(cell)=(T)0;momentum(cell)=TV();
            velocity(cell)=TV();entropy(cell)=(T)0;
            speedofsound(cell)=(T)0;internal_energy(cell)=(T)0;}}
    euler.Get_Variable_Gradient(density_gradient,density);euler.Get_Variable_Gradient(pressure_gradient,pressure);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/density."+f,density);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/density_gradient."+f,density_gradient);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/centered_velocities."+f,velocity);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/pressure."+f,pressure);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/pressure_gradient."+f,pressure_gradient);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/entropy."+f,entropy);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/speedofsound."+f,speedofsound);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/machnumber."+f,machnumber);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/momentum."+f,momentum);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/energy."+f,energy);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/internal_energy."+f,internal_energy);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/velocity_plus_c."+f,velocity_plus_c);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/velocity_minus_c."+f,velocity_minus_c);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/mac_velocities."+f,face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/psi_D."+f,euler.euler_projection.elliptic_solver->psi_D);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/psi_N."+f,euler.euler_projection.elliptic_solver->psi_N);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void AERO_INTERFACE_1<T>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(output_directory);
    Write_Levelset("levelset.phi");
    Write_Triangulated_Surface("triangulated_surface.tri",*rigid_body->simplicial_object);
    Write_First_Frame(frame);
    Write_Solid_Output_Files(frame);
    Write_Fluid_Output_Files(frame);
    Write_Last_Frame(frame);
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class T> void AERO_INTERFACE_1<T>::
Write_Substep(const std::string& title,const int level)
{
    LOG::cout<<"Output directory: "<<output_directory<<std::endl;
    if(level<=substeps_level){
        LOG::cout<<"Writing substep ["<<title<<"]: output_number="<<output_number+1<<std::endl;
        Write_Output_Files(++output_number);}
}
//#####################################################################
template class AERO_INTERFACE_1<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class AERO_INTERFACE_1<double>;
#endif
