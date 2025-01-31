//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX_4X4.h>
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Deformable_Object
//#####################################################################
template<class T> int THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Add_Deformable_Object(DEFORMABLE_BODY_COLLECTION<VECTOR<T,2> >& deformable_body_collection,const int number_of_vertices,const VECTOR<T,2>& start_position,const VECTOR<T,2>& end_position)
{
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(deformable_body_collection.particles);
    int index=deformable_body_collection.Add_Structure(&segmented_curve);
    SEGMENT_MESH& mesh=segmented_curve.mesh;
    GEOMETRY_PARTICLES<VECTOR<T,2> >& particles=segmented_curve.particles;
    mesh.Initialize_Straight_Mesh(number_of_vertices);
    for(int i=0;i<number_of_vertices;i++){
        particles.Add_Element();assert(particles.Size()==i+1);
        particles.X(i)=start_position+((T)i/(number_of_vertices-1))*(end_position-start_position);
        particles.V(i)=VECTOR<T,2>(0,0);}
    return index;
}
//#####################################################################
// Function Add_Circle_Deformable_Object
//#####################################################################
template<class T> int THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Add_Circle_Deformable_Object(DEFORMABLE_BODY_COLLECTION<VECTOR<T,2> >& deformable_body_collection,const int number_of_vertices,const VECTOR<T,2>& center,const T radius)
{
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(deformable_body_collection.particles);
    int index=deformable_body_collection.Add_Structure(&segmented_curve);
    SEGMENT_MESH& mesh=segmented_curve.mesh;
    GEOMETRY_PARTICLES<VECTOR<T,2> >& particles=segmented_curve.particles;
    mesh.Initialize_Straight_Mesh(number_of_vertices,true);
    for(int i=0;i<number_of_vertices;i++){
        particles.Add_Element();assert(particles.Size()==i+1);
        T angle=i*(T)2*(T)pi/number_of_vertices;
        particles.X(i)=center+radius*VECTOR<T,2>(cos(angle+(T).5*(T)pi),sin(angle+(T).5*(T)pi));
        particles.V(i)=VECTOR<T,2>(0,0);}
    return index;
}
//#####################################################################
// Function Add_Grid_Deformable_Object
//#####################################################################
template<class T> int THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Add_Grid_Deformable_Object(DEFORMABLE_BODY_COLLECTION<VECTOR<T,2> >& deformable_body_collection,const GRID<VECTOR<T,2> >& grid,const int edge_subdivision)
{
    SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(deformable_body_collection.particles);
    int index=deformable_body_collection.Add_Structure(&segmented_curve);
    SEGMENT_MESH& mesh=segmented_curve.mesh;
    GEOMETRY_PARTICLES<VECTOR<T,2> >& particles=segmented_curve.particles;
    mesh.number_nodes=grid.counts.x*grid.counts.y+(edge_subdivision-1)*((grid.counts.x-1)*grid.counts.y+(grid.counts.y-1)*grid.counts.x);
    mesh.elements.Exact_Resize(edge_subdivision*((grid.counts.x-1)*grid.counts.y+(grid.counts.y-1)*grid.counts.x));
    int segment_index=0;
    ARRAY<int,VECTOR<int,2> > base_index(grid.Domain_Indices());
    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++){
        int index=base_index(i,j)=particles.Add_Element();
        particles.X(index)=grid.X(VECTOR<int,2>(i,j));particles.V(index)=VECTOR<T,2>(0,0);}
    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++){
        if(i<grid.counts.x){
            int last_index=base_index(i,j);
            for(int k=0;k<edge_subdivision-1;k++){
                int index=particles.Add_Element();
                particles.X(index)=((T)k/edge_subdivision)*VECTOR<T,2>(grid.dX.x,0)+grid.X(VECTOR<int,2>(i,j));particles.V(index)=VECTOR<T,2>(0,0);
                mesh.elements(segment_index++).Set(last_index,index);
                last_index=index;}
            mesh.elements(segment_index++).Set(last_index,base_index(i+1,j));}
        if(j<grid.counts.y){
            int last_index=base_index(i,j);
            for(int k=0;k<edge_subdivision-1;k++){
                int index=particles.Add_Element();
                particles.X(index)=((T)k/edge_subdivision)*VECTOR<T,2>(0,grid.dX.y)+grid.X(VECTOR<int,2>(i,j));particles.V(index)=VECTOR<T,2>(0,0);
                mesh.elements(segment_index++).Set(last_index,index);
                last_index=index;}
            mesh.elements(segment_index++).Set(last_index,base_index(i,j+1));}}
    return index;
}
//#####################################################################
// Function Add_Deformable_Object
//#####################################################################
template<class T> int THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Add_Deformable_Object(DEFORMABLE_BODY_COLLECTION<VECTOR<T,3> >& deformable_body_collection,ARRAY<int>& deformable_body_collection_enslaved_nodes,const GRID<VECTOR<T,2> >& cloth_grid,const MATRIX<T,4>& transform,
    const int constraint_mode)
{
    TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create(deformable_body_collection.particles);
    int index=deformable_body_collection.Add_Structure(&triangulated_surface);
    TRIANGLE_MESH& mesh=triangulated_surface.mesh;
    GEOMETRY_PARTICLES<VECTOR<T,3> >& particles=triangulated_surface.particles;

    mesh.Initialize_Herring_Bone_Mesh(cloth_grid.counts.x,cloth_grid.counts.y);
    particles.Add_Elements(mesh.number_nodes);
    for(int i=0;i<cloth_grid.counts.x;i++) for(int j=0;j<cloth_grid.counts.y;j++){
        int node=i+cloth_grid.counts.x*(j-1);particles.X(node)=transform.Homogeneous_Times(VECTOR<T,3>(cloth_grid.X(VECTOR<int,2>(i,j))));particles.V(node)=VECTOR<T,3>();}

    if(constraint_mode==1){
        deformable_body_collection_enslaved_nodes.Append(1+cloth_grid.counts.x*(cloth_grid.counts.y-1));
        deformable_body_collection_enslaved_nodes.Append(cloth_grid.counts.x*cloth_grid.counts.y);}
    else if(constraint_mode==2){
        deformable_body_collection_enslaved_nodes.Append(1+cloth_grid.counts.x*(cloth_grid.counts.y-1));
        deformable_body_collection_enslaved_nodes.Append(cloth_grid.counts.x*cloth_grid.counts.y);
        deformable_body_collection_enslaved_nodes.Append(cloth_grid.counts.x);
        deformable_body_collection_enslaved_nodes.Append(1);
    }

    return index;
}
//#####################################################################
// Function Add_Deformable_From_File
//#####################################################################
template<class T> int THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Add_Deformable_Object_From_File(const STREAM_TYPE stream_type,DEFORMABLE_BODY_COLLECTION<VECTOR<T,3> >& deformable_body_collection,ARRAY<int>& deformable_body_collection_enslaved_nodes,const std::string& filename,
    const MATRIX<T,4>& transform,PLANE<T>* enslaved_halfplane)
{
    TRIANGULATED_SURFACE<T>& triangulated_surface=*TRIANGULATED_SURFACE<T>::Create(deformable_body_collection.particles);
    int index=deformable_body_collection.Add_Structure(&triangulated_surface);
    DEFORMABLE_PARTICLES<VECTOR<T,3> >& particles=dynamic_cast<DEFORMABLE_PARTICLES<VECTOR<T,3> >&>(triangulated_surface.particles);

    particles.Store_Velocity(false);particles.Store_Mass(false); // need to do this before reading it in
    Read_From_File(filename,triangulated_surface);
    particles.Store_Velocity(true);particles.Store_Mass(true);

    if(enslaved_halfplane){for(int i=0;i<particles.Size();i++) if(enslaved_halfplane->Lazy_Inside(particles.X(i))) deformable_body_collection_enslaved_nodes.Append(i);}

    for(int i=0;i<particles.Size();i++){particles.X(i)=transform.Homogeneous_Times(particles.X(i));particles.V(i)=VECTOR<T,3>();}
    triangulated_surface.Refresh_Auxiliary_Structures();

    return index;
}
//#####################################################################
// Function Add_Rigid_Body_Walls
//#####################################################################
template<class T,class T_EXAMPLE> static void
Add_Rigid_Body_Walls(T_EXAMPLE& example,const GRID<VECTOR<T,2> >& fluid_grid,const T coefficient_of_restitution,const T coefficient_of_friction,ARRAY<int>* walls_added)
{
    typedef VECTOR<T,2> TV;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=example.solid_body_collection.rigid_body_collection;
    VECTOR<T,2> center=fluid_grid.domain.Center(),size=fluid_grid.domain.Edge_Lengths();
    int id;

    if(example.fluids_parameters.domain_walls(0)(0)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies_2D/ground",size.y*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id)=FRAME<TV>(TV(fluid_grid.domain.min_corner.x,center.y),ROTATION<TV>::From_Angle(-(T)pi/2));
        rigid_body_collection.Rigid_Body(id).name="left wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(example.fluids_parameters.domain_walls(0)(1)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies_2D/ground",size.y*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id)=FRAME<TV>(TV(fluid_grid.domain.max_corner.x,center.y),ROTATION<TV>::From_Angle((T)pi/2));
        rigid_body_collection.Rigid_Body(id).name="right wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(example.fluids_parameters.domain_walls(1)(0)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies_2D/ground",size.x*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id).t=TV(center.x,fluid_grid.domain.min_corner.y);
        rigid_body_collection.Rigid_Body(id).name="bottom wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(example.fluids_parameters.domain_walls(1)(1)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies_2D/ground",size.x*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id)=FRAME<TV>(TV(center.x,fluid_grid.domain.max_corner.y),ROTATION<TV>::From_Angle((T)pi));
        rigid_body_collection.Rigid_Body(id).name="top wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}
}
//#####################################################################
// Function Add_Rigid_Body_Walls
//#####################################################################
template<class T,class T_EXAMPLE> static void
Add_Rigid_Body_Walls(T_EXAMPLE& example,const GRID<VECTOR<T,3> >& fluid_grid,const T coefficient_of_restitution,const T coefficient_of_friction,ARRAY<int>* walls_added)
{
    typedef VECTOR<T,3> TV;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=example.solid_body_collection.rigid_body_collection;
    TV center=fluid_grid.domain.Center(),size=fluid_grid.domain.Edge_Lengths();
    int id;

    if(example.fluids_parameters.domain_walls(0)(0)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies/ground",max(size.y,size.z)*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id)=FRAME<TV>(TV(fluid_grid.domain.min_corner.x,center.y,center.z),ROTATION<TV>(-(T).5*(T)pi,TV(0,0,1)));
        rigid_body_collection.Rigid_Body(id).name="left wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(example.fluids_parameters.domain_walls(0)(1)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies/ground",max(size.y,size.z)*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id)=FRAME<TV>(TV(fluid_grid.domain.max_corner.x,center.y,center.z),ROTATION<TV>((T).5*(T)pi,TV(0,0,1)));
        rigid_body_collection.Rigid_Body(id).name="right wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(example.fluids_parameters.domain_walls(1)(0)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies/ground",max(size.x,size.z)*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id)=FRAME<TV>(TV(center.x,fluid_grid.domain.min_corner.y,center.z));
        rigid_body_collection.Rigid_Body(id).name="bottom wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(example.fluids_parameters.domain_walls(1)(1)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies/ground",max(size.x,size.z)*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id)=FRAME<TV>(TV(center.x,fluid_grid.domain.max_corner.y,center.z),ROTATION<TV>((T)pi,TV(1,0,0)));
        rigid_body_collection.Rigid_Body(id).name="top wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(example.fluids_parameters.domain_walls(2)(0)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies/ground",max(size.x,size.y)*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id)=FRAME<TV>(TV(center.x,center.y,fluid_grid.domain.min_corner.z),ROTATION<TV>((T).5*(T)pi,TV(1,0,0)));
        rigid_body_collection.Rigid_Body(id).name="front wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}

    if(example.fluids_parameters.domain_walls(2)(1)){
        id=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies/ground",max(size.x,size.y)*(T).00501);
        rigid_body_collection.Rigid_Body(id).Set_Coefficient_Of_Restitution(coefficient_of_restitution);
        rigid_body_collection.Rigid_Body(id).coefficient_of_friction=coefficient_of_friction;
        rigid_body_collection.rigid_body_particles.frame(id)=FRAME<TV>(TV(center.x,center.y,fluid_grid.domain.max_corner.z),ROTATION<TV>((T)-.5*(T)pi,TV(1,0,0)));
        rigid_body_collection.Rigid_Body(id).name="back wall";
        rigid_body_collection.Rigid_Body(id).is_static=true;
        if(walls_added) walls_added->Append(id);}
}
//#####################################################################
// Function Add_Rigid_Body_Walls
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Add_Rigid_Body_Walls(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T,2> >& example,const T coefficient_of_restitution,const T coefficient_of_friction,ARRAY<int>* walls_added)
{
    ::Add_Rigid_Body_Walls(example,(example.fluids_parameters.mpi_grid?example.fluids_parameters.mpi_grid->global_grid:*example.fluids_parameters.grid),coefficient_of_restitution,coefficient_of_friction,walls_added);
}
//#####################################################################
// Function Add_Rigid_Body_Walls
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Add_Rigid_Body_Walls(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T,3> >& example,const T coefficient_of_restitution,const T coefficient_of_friction,ARRAY<int>* walls_added)
{
    ::Add_Rigid_Body_Walls(example,(example.fluids_parameters.mpi_grid?example.fluids_parameters.mpi_grid->global_grid:*example.fluids_parameters.grid),coefficient_of_restitution,coefficient_of_friction,walls_added);
}
//#####################################################################
// Function Set_Deformable_Object_Parameters_2D
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Set_Deformable_Object_Parameters_2D(const int id,T& density,PARAMETER_LIST& parameter_list)
{
    std::string prefix=LOG::sprintf("deformable_object_%d",id);
    parameter_list.Get_Parameter_In_Place(prefix+".density",density);
}
//#####################################################################
// Function Set_Deformable_Object_Parameters_3D
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Set_Deformable_Object_Parameters_3D(const int id,T& edge_stiffness_scaling,T& altitude_stiffness_scaling,T& density,PARAMETER_LIST& parameter_list)
{
    std::string prefix=LOG::sprintf("deformable_object_%d",id);
    parameter_list.Get_Parameter_In_Place(prefix+".edge_stiffness_scaling",edge_stiffness_scaling);
    parameter_list.Get_Parameter_In_Place(prefix+".altitude_stiffness_scaling",altitude_stiffness_scaling);
    parameter_list.Get_Parameter_In_Place(prefix+".density",density);
}
//#####################################################################
// Function Set_Rigid_Body_Parameters_2D
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Set_Rigid_Body_Parameters_2D(const int id,T& density,PARAMETER_LIST& parameter_list)
{
    std::string prefix=LOG::sprintf("rigid_body_%d",id);
    parameter_list.Get_Parameter_In_Place(prefix+".density",density);
}
//#####################################################################
// Function Set_Rigid_Body_Parameters_3D
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Set_Rigid_Body_Parameters_3D(const int id,T& density,PARAMETER_LIST& parameter_list)
{
    std::string prefix=LOG::sprintf("rigid_body_%d",id);
    parameter_list.Get_Parameter_In_Place(prefix+".density",density);
}
//#####################################################################
// Function Set_Mass
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Set_Mass(TRIANGULATED_SURFACE<T>& triangulated_surface,const T mass,const bool use_constant_mass)
{
    T density=mass/triangulated_surface.Total_Area();
    SOLIDS_STANDARD_TESTS<VECTOR<T,3> >::Set_Mass_Of_Particles(triangulated_surface,density,use_constant_mass);
}
//#####################################################################
// Function Set_Density
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Set_Density(TRIANGULATED_SURFACE<T>& triangulated_surface,const T density,const bool use_constant_mass)
{
    SOLIDS_STANDARD_TESTS<VECTOR<T,3> >::Set_Mass_Of_Particles(triangulated_surface,density,use_constant_mass);
}
//#####################################################################
// Function Set_Mass
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Set_Mass(SEGMENTED_CURVE_2D<T>& segmented_curve,const T mass,const bool use_constant_mass)
{
    T density=mass/segmented_curve.Total_Length();
    SOLIDS_STANDARD_TESTS<VECTOR<T,2> >::Set_Mass_Of_Particles(segmented_curve,density,use_constant_mass);
}
//#####################################################################
// Function Set_Density
//#####################################################################
template<class T> void THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::
Set_Density(SEGMENTED_CURVE_2D<T>& segmented_curve,const T density,const bool use_constant_mass)
{
    SOLIDS_STANDARD_TESTS<VECTOR<T,2> >::Set_Mass_Of_Particles(segmented_curve,density,use_constant_mass);
}
//#####################################################################
template<class T_EXAMPLE> static void Set_Example_Common_Parameters_From_Parameter_List(T_EXAMPLE& example,PARAMETER_LIST& parameter_list)
{
    parameter_list.Get_Parameter_In_Place("last_frame",example.last_frame);
    parameter_list.Get_Parameter_In_Place("frame_rate",example.frame_rate);
    parameter_list.Get_Parameter_In_Place("restart",example.restart);
    parameter_list.Get_Parameter_In_Place("verbose_dt",example.solids_parameters.verbose_dt);
    //parameter_list.Get_Parameter_In_Place("use_collision_aware_velocity_extrapolation",example.use_collision_aware_velocity_extrapolation);
    //parameter_list.Get_Parameter_In_Place("use_collision_aware_signed_distance",example.use_collision_aware_signed_distance);
    //parameter_list.Get_Parameter_In_Place("clamp_phi_with_collision_bodies",example.clamp_phi_with_collision_bodies);
}

template<class TV> static void Set_Parameters_From_Parameter_List(FLUIDS_PARAMETERS<TV>& fluids_parameters,PARAMETER_LIST& parameter_list)
{
    parameter_list.Get_Parameter_In_Place("fluids_parameters.cfl",fluids_parameters.cfl);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.use_vorticity_confinement",fluids_parameters.use_vorticity_confinement);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.use_vorticity_confinement_fuel",fluids_parameters.use_vorticity_confinement_fuel);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.confinement_parameter",fluids_parameters.confinement_parameter);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.confinement_parameter_fuel",fluids_parameters.confinement_parameter_fuel);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.temperature_buoyancy_constant",fluids_parameters.temperature_buoyancy_constant);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.density_buoyancy_constant",fluids_parameters.density_buoyancy_constant);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.density_buoyancy_threshold",fluids_parameters.density_buoyancy_threshold);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.normal_flame_speed",fluids_parameters.normal_flame_speed);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.density",fluids_parameters.density);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.density_fuel",fluids_parameters.density_fuel);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.incompressible_iterations",fluids_parameters.incompressible_iterations);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.incompressible_tolerance",fluids_parameters.incompressible_tolerance);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.maximum_tree_depth",fluids_parameters.maximum_tree_depth);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.write_debug_data",fluids_parameters.write_debug_data);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.write_ghost_values",fluids_parameters.write_ghost_values);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.restart_data_write_rate",fluids_parameters.restart_data_write_rate);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.simulate",fluids_parameters.simulate);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.number_particles_per_cell",fluids_parameters.number_particles_per_cell);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.thin_shells_refine_near_objects",fluids_parameters.thin_shells_refine_near_objects);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.delete_fluid_inside_objects",fluids_parameters.delete_fluid_inside_objects);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.object_friction",fluids_parameters.object_friction);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.second_order_cut_cell_method",fluids_parameters.second_order_cut_cell_method);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.solve_neumann_regions",fluids_parameters.solve_neumann_regions);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.particle_half_bandwidth",fluids_parameters.particle_half_bandwidth);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.use_external_velocity",fluids_parameters.use_external_velocity);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.refine_fmm_initialization_with_iterative_solver",fluids_parameters.refine_fmm_initialization_with_iterative_solver);
    parameter_list.Get_Parameter_In_Place("fluids_parameters.modify_wall_tangential_velocities",fluids_parameters.modify_wall_tangential_velocities);

    parameter_list.Get_Parameter_In_Place("fluids_parameters.domain_walls(0)(0)",fluids_parameters.domain_walls(0)(0));
    parameter_list.Get_Parameter_In_Place("fluids_parameters.domain_walls(0)(1)",fluids_parameters.domain_walls(0)(1));
    parameter_list.Get_Parameter_In_Place("fluids_parameters.domain_walls(1)(0)",fluids_parameters.domain_walls(1)(0));
    parameter_list.Get_Parameter_In_Place("fluids_parameters.domain_walls(1)(1)",fluids_parameters.domain_walls(1)(1));
    parameter_list.Get_Parameter_In_Place("fluids_parameters.domain_walls(2)(0)",fluids_parameters.domain_walls(2)(0));
    parameter_list.Get_Parameter_In_Place("fluids_parameters.domain_walls(2)(1)",fluids_parameters.domain_walls(2)(1));

    if(parameter_list.Is_Defined("fluids_parameters.levelset_substeps")) PHYSBAM_FATAL_ERROR("levelset_substeps has been renamed to scalar_substeps -- update your param file!");

    parameter_list.Get_Parameter_In_Place("fluids_parameters.reincorporate_removed_particle_velocity",fluids_parameters.reincorporate_removed_particle_velocity);
}
//#####################################################################
namespace PhysBAM{
template class THIN_SHELLS_FLUID_COUPLING_UTILITIES<float>;
template class THIN_SHELLS_FLUID_COUPLING_UTILITIES<double>;
}
