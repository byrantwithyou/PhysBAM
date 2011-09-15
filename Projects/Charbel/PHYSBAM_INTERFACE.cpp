//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <climits>
#include "PHYSBAM_INTERFACE.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> PhysBAMInterface<T>::
PhysBAMInterface(TRIANGULATED_SURFACE<T> &triangulated_surface_input)
    :triangle_mesh(triangulated_surface_input.mesh),particles(triangulated_surface_input.particles),triangulated_surface(&triangulated_surface_input)
{
    triangle_hierarchy=new TRIANGLE_HIERARCHY<T>(triangle_mesh,particles,true,0);
    triangulated_surface->hierarchy=triangle_hierarchy;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> PhysBAMInterface<T>::
PhysBAMInterface(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :triangle_mesh(triangle_mesh_input),particles(particles_input),triangulated_surface(new TRIANGULATED_SURFACE<T>(triangle_mesh_input,particles_input))
{
    triangulated_surface->Update_Triangle_List();
    triangle_hierarchy=new TRIANGLE_HIERARCHY<T>(triangle_mesh,particles,true,0);
    triangulated_surface->hierarchy=triangle_hierarchy;
    // collision_geometry=new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(triangulated_surface);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> PhysBAMInterface<T>::
~PhysBAMInterface()
{
    delete triangle_hierarchy;
    delete triangulated_surface;
}
//#####################################################################
// Update
//#####################################################################
template<class T> void PhysBAMInterface<T>::
Update(const bool rebuild_hierarchy)
{
    if(triangulated_surface->triangle_list) triangulated_surface->Update_Triangle_List();
    if(rebuild_hierarchy){
        delete triangle_hierarchy;
        triangle_hierarchy=new TRIANGLE_HIERARCHY<T>(triangle_mesh,particles,true,0);
        triangulated_surface->hierarchy=triangle_hierarchy;}
    else triangle_hierarchy->Update_Boxes();
}
//#####################################################################
// HasCloseTriangle
//#####################################################################
template<class T> bool PhysBAMInterface<T>::
HasCloseTriangle(const TV min_corner,const TV max_corner,bool* is_occluded) const
{
    RANGE<TV> bounding_box(min_corner, max_corner);
    ARRAY<int> candidates;

    triangle_hierarchy->Intersection_List(bounding_box,candidates,1e-5);
    return candidates.Size()>0;
}
//#####################################################################
// HasCloseTriangle
//#####################################################################
template<class T> bool PhysBAMInterface<T>::
HasCloseTriangle(const TV position,const TV min_corner,const TV max_corner,const T thickness,bool* is_occluded) const
{
    RANGE<TV> bounding_box(min_corner, max_corner);
    ARRAY<int> candidates;
    triangle_hierarchy->Intersection_List(bounding_box,candidates,thickness);

    if(is_occluded){*is_occluded=false;
        for(int i=1;i<=candidates.Size() && !*is_occluded;++i) *is_occluded=triangulated_surface->Get_Element(candidates(i)).Point_Inside_Triangle(position,thickness*(T).5);}

    return candidates.Size()>0;
}
//#####################################################################
// GetClosestTriangle
//#####################################################################
template<class T> int PhysBAMInterface<T>::
GetClosestTriangle(const TV position,const TV min_corner,const TV max_corner,T* distance) const
{
    RANGE<TV> bounding_box(min_corner, max_corner);
    ARRAY<int> candidates;

    triangle_hierarchy->Intersection_List(bounding_box,candidates,1e-5);
    T dist=(T)FLT_MAX;int closest_triangle=-1;
    for(int i=1;i<=candidates.Size();++i){
        T local_dist=triangulated_surface->Get_Element(candidates(i)).Distance_To_Triangle(position);
        if(local_dist < dist){
            dist=local_dist;closest_triangle=candidates(i);}}

    if(closest_triangle >= 0 && distance) *distance=(static_cast<PLANE<T> >(
        triangulated_surface->Get_Element(closest_triangle)).Inside(position,(T)0)?-dist:dist);
    return closest_triangle;
}
//#####################################################################
// Intersect
//#####################################################################
template<class T> IntersectionResult<T> PhysBAMInterface<T>::
Intersect(const TV& start,const TV& end,const T thickness) const
{
    IntersectionResult<T> result;
    RAY<TV> edge(SEGMENT_3D<T>(start,end));
    TV weights;T old_t_max=edge.t_max;
    if(triangle_hierarchy->Intersection(edge,weights,thickness*(T).5)){
        result.triangleID=edge.aggregate_id;
        result.zeta[0]=weights(1);result.zeta[1]=weights(2);result.zeta[2]=weights(3);
        result.alpha=(T)1-edge.t_max/old_t_max;}
    else
        result.triangleID=-1;
    return result;
}
//#####################################################################
// Intersect
//#####################################################################
template<class T> void PhysBAMInterface<T>::
Intersect(const ARRAY<TV>& node_positions,ARRAY<bool>& occluded_node,ARRAY<TRIPLE<VECTOR<int,3>,IntersectionResult<T>,IntersectionResult<T> > >& edges_and_results,const T thickness) const
{
    // Detect occluded nodes as a separate pass
    int triangle_id; for(int i=1;i<=node_positions.m;++i) occluded_node(i)=triangulated_surface->Inside_Any_Triangle(node_positions(i),triangle_id,thickness*(T).5);

    const int PERMITTED_ATTEMPTS=10;
    int retryAttempts=0;
    for(int i=1;i<=edges_and_results.m;++i){
        VECTOR<int,3>& edge_nodes=edges_and_results(i).x;
        bool left_node_occluded=occluded_node(edge_nodes(1)),right_node_occluded=occluded_node(edge_nodes(2));
        IntersectionResult<T>& current_working_result=edges_and_results(i).y;
        IntersectionResult<T>& reverse_working_result=edges_and_results(i).z;

        RAY<TV> edge(SEGMENT_3D<T>(node_positions(edge_nodes(1)),node_positions(edge_nodes(2)))),
            reverse_edge(SEGMENT_3D<T>(node_positions(edge_nodes(2)),node_positions(edge_nodes(1))));
        TV weights(0,0,0),reverse_weights(0,0,0);
        T old_t_max=edge.t_max,
          reverse_old_t_max=reverse_edge.t_max;
        if(!left_node_occluded){
            triangle_hierarchy->Intersection(edge,weights,thickness*(T).5);
            if(right_node_occluded && edge.intersection_location==RAY<TV>::LOCATION_UNKNOWN){ // This edge connects to an occluded node, and thus MUST report an intersection.
                for(int number_of_attempts=0;number_of_attempts<PERMITTED_ATTEMPTS && edge.intersection_location==RAY<TV>::LOCATION_UNKNOWN;++number_of_attempts){
                    edge.t_max+=thickness;old_t_max+=thickness*(T)(1 << number_of_attempts);triangle_hierarchy->Intersection(edge,weights,thickness*(T).5);++retryAttempts;}
                if(edge.intersection_location==RAY<TV>::LOCATION_UNKNOWN){LOG::cout<<"An occluded node appears to be visible... BUG 1!"<<std::endl;exit(-1);}}
            else if(edge.intersection_location==RAY<TV>::START_POINT){ // This edge should not believe it starts at an occluded point, so shrink the thickness a little.
                T modified_thickness=thickness;
                for(int number_of_attempts=0;number_of_attempts<PERMITTED_ATTEMPTS && edge.intersection_location==RAY<TV>::START_POINT;++number_of_attempts){
                    modified_thickness*=(T).5;triangle_hierarchy->Intersection(edge,weights,modified_thickness*(T).5);++retryAttempts;}
                if(edge.intersection_location==RAY<TV>::START_POINT){LOG::cout<<"A visible node appears to be occluded... BUG!"<<std::endl;exit(-1);}}

            if(edge.intersection_location==RAY<TV>::LOCATION_UNKNOWN) current_working_result.triangleID=-1; // Truely no intersection
            else{
                current_working_result.triangleID=edge.aggregate_id;current_working_result.alpha=(T)1-edge.t_max/old_t_max;
                current_working_result.zeta[0]=weights(1);current_working_result.zeta[1]=weights(2);current_working_result.zeta[2]=weights(3);}}
        else current_working_result.triangleID=INT_MAX; // TODO(jontg): HACK - Report an intersection that should never be used.

        if(!right_node_occluded){
            triangle_hierarchy->Intersection(reverse_edge,reverse_weights,thickness*(T).5);
            if(left_node_occluded && reverse_edge.intersection_location==RAY<TV>::LOCATION_UNKNOWN){ // This edge connects to an occluded node, and thus MUST report an intersection.
                for(int number_of_attempts=0;number_of_attempts<PERMITTED_ATTEMPTS && reverse_edge.intersection_location==RAY<TV>::LOCATION_UNKNOWN;++number_of_attempts){
                    reverse_edge.t_max+=thickness;reverse_old_t_max+=thickness;triangle_hierarchy->Intersection(reverse_edge,reverse_weights,thickness*(T).5);++retryAttempts;}
                if(reverse_edge.intersection_location==RAY<TV>::LOCATION_UNKNOWN){LOG::cout<<"An occluded node appears to be visible... BUG 2!"<<std::endl;exit(-1);}}
            else if(reverse_edge.intersection_location==RAY<TV>::START_POINT){ // This reverse_edge should not believe it starts at an occluded point, so shrink the thickness a little.
                T modified_thickness=thickness;
                for(int number_of_attempts=0;number_of_attempts<PERMITTED_ATTEMPTS && reverse_edge.intersection_location==RAY<TV>::START_POINT;++number_of_attempts){
                    modified_thickness*=(T).5;triangle_hierarchy->Intersection(reverse_edge,reverse_weights,modified_thickness*(T).5);++retryAttempts;}
                if(reverse_edge.intersection_location==RAY<TV>::START_POINT){LOG::cout<<"A visible node appears to be occluded... BUG!"<<std::endl;exit(-1);}}
            
            if(reverse_edge.intersection_location==RAY<TV>::LOCATION_UNKNOWN) reverse_working_result.triangleID=-1; // Truely no intersection
            else{
                reverse_working_result.triangleID=reverse_edge.aggregate_id;reverse_working_result.alpha=(T)1-reverse_edge.t_max/reverse_old_t_max;
                reverse_working_result.zeta[0]=reverse_weights(1);reverse_working_result.zeta[1]=reverse_weights(2);reverse_working_result.zeta[2]=reverse_weights(3);}}
        else reverse_working_result.triangleID=INT_MAX;} // TODO(jontg): HACK - Report an intersection that should never be used.

    if(retryAttempts != 0) LOG::cout<<"Intersect Finished with "<<retryAttempts<<" retry attempts"<<std::endl;
}
//#####################################################################
// computeSweptNodes
//#####################################################################
template<class T> void PhysBAMInterface<T>::
computeSweptNodes(const ARRAY<TV>& node_positions,ARRAY<bool>& swept_node,const T dt,const T thickness)
{
    // TODO(jontg): In order to work, we need to...
    //   1. Initialize and keep around collision_geometry of type DEFORMABLE_OBJECT_FLUID_COLLISIONS
    //   2. Keep the saved states updated properly
    // for(int i=1;i<=node_positions.m;++i)
    //     swept_node(i)=collision_geometry.Any_Simplex_Crossover(node_positions(i),node_positions(i),dt);
}
//#####################################################################
template class PhysBAMInterface<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PhysBAMInterface<double>;
#endif
