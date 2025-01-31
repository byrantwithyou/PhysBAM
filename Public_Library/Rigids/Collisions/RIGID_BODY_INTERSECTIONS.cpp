//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/PAIR.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_PARTITION.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/PARTICLE_LEVELSET_VISITOR.h>
#include <Rigids/Collisions/PARTICLES_IN_IMPLICIT_OBJECT.h>
#include <Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>

using namespace PhysBAM;
//#####################################################################
// Function Destructor
//#####################################################################
template<class TV> RIGID_BODY_INTERSECTIONS<TV>::
~RIGID_BODY_INTERSECTIONS()
{
    typedef typename HASHTABLE<const GEOMETRY_PARTICLES<TV>*,const PARTICLE_HIERARCHY<TV>*>::ITERATOR T_HIERARCHY_ITERATOR;
    for(T_HIERARCHY_ITERATOR i(particle_hierarchies);i.Valid();i.Next()) delete i.Data();
}
//#####################################################################
// Function Intersection_Check
//#####################################################################
template<class TV> bool RIGID_BODY_INTERSECTIONS<TV>::
Intersection_Check(const int id_1,const int id_2,int& particle_body,int& levelset_body,const T thickness) const
{
    if(!Bounding_Boxes_Intersect(id_1,id_2,thickness)) return false;
    else return Find_Any_Intersection(id_1,id_2,particle_body,levelset_body);
}
//#####################################################################
// Function Bounding_Boxes_Intersect
//#####################################################################
template<class TV> bool RIGID_BODY_INTERSECTIONS<TV>::
Bounding_Boxes_Intersect(const int id_1,const int id_2,const T thickness) const
{
    RIGID_BODY<TV> &body_1=rigid_body_collection.Rigid_Body(id_1),&body_2=rigid_body_collection.Rigid_Body(id_2);
    if(!body_1.Is_Simulated() && !body_2.Is_Simulated()) return false; // don't check when neither object is dynamic
    return body_1.Bounding_Boxes_Intersect(body_2,thickness);
}
//#####################################################################
// Function Find_Any_Intersection
//#####################################################################
template<class TV> bool RIGID_BODY_INTERSECTIONS<TV>::
Find_Any_Intersection(const int id_1,const int id_2,int& particle_body,int& levelset_body) const
{
    RIGID_BODY<TV> &body0=rigid_body_collection.Rigid_Body(id_1),&body1=rigid_body_collection.Rigid_Body(id_2);
    Initialize_Transformation_From_Body1_To_Body2_Coordinates(body0,body1);
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> > intersection_list;
    if(!use_triangle_hierarchy){
        if(body0.simplicial_object && body1.implicit_object){
            Particles_In_Levelset(id_1,id_2,intersection_list,0,true);
            if(intersection_list.m){particle_body=id_1;levelset_body=id_2;return true;}}
        if(body0.implicit_object && body1.simplicial_object){
            Flip_Transformation();Particles_In_Levelset(id_2,id_1,intersection_list,0,true);
            if(intersection_list.m){particle_body=id_2;levelset_body=id_1;return true;}}}
    else{ // using triangle hierarchy
        if(!use_edge_intersection){
            ARRAY<int> triangle_list;triangle_list.Preallocate(50);
            if(body0.simplicial_object && body1.implicit_object){
                PARTICLES_IN_IMPLICIT_OBJECT<TV>::Get_Interfering_Simplices(body0,body1,triangle_list,rotation,translation,use_triangle_hierarchy_center_phi_test);
                PARTICLES_IN_IMPLICIT_OBJECT<TV>::Intersections_Using_Hierarchy(body0,body1,triangle_list,intersection_list,0,true,rotation,translation);
                if(intersection_list.m){particle_body=id_1;levelset_body=id_2;return true;}}
            if(body0.implicit_object && body1.simplicial_object){
                Flip_Transformation();PARTICLES_IN_IMPLICIT_OBJECT<TV>::Get_Interfering_Simplices(body1,body0,triangle_list,rotation,translation,use_triangle_hierarchy_center_phi_test);
                PARTICLES_IN_IMPLICIT_OBJECT<TV>::Intersections_Using_Hierarchy(body1,body0,triangle_list,intersection_list,0,true,rotation,translation);
                if(intersection_list.m){particle_body=id_2;levelset_body=id_1;return true;}}}
        else{ // use edge intersection
            ARRAY<int> triangle_list1,triangle_list2;triangle_list1.Preallocate(50);triangle_list2.Preallocate(50);
            if(body0.simplicial_object) PARTICLES_IN_IMPLICIT_OBJECT<TV>::Get_Interfering_Simplices(body0,body1,triangle_list1,rotation,translation,use_triangle_hierarchy_center_phi_test);
            Flip_Transformation();
            if(body1.simplicial_object){PARTICLES_IN_IMPLICIT_OBJECT<TV>::Get_Interfering_Simplices(body1,body0,triangle_list2,rotation,translation,use_triangle_hierarchy_center_phi_test);
                if(body0.implicit_object){PARTICLES_IN_IMPLICIT_OBJECT<TV>::Intersections_Using_Hierarchy_And_Edges(body1,body0,triangle_list2,triangle_list1,intersection_list,0,true,rotation,translation);
                    if(intersection_list.m){particle_body=id_2;levelset_body=id_1;return true;}}}
            if(body0.simplicial_object && body1.implicit_object){
                Flip_Transformation();
                PARTICLES_IN_IMPLICIT_OBJECT<TV>::Intersections_Using_Hierarchy_And_Edges(body0,body1,triangle_list1,triangle_list2,intersection_list,0,true,rotation,translation);
                if(intersection_list.m){particle_body=id_1;levelset_body=id_2;return true;}}}}
    return false;
}
//#####################################################################
// Function Append_All_Intersections
//#####################################################################
template<class TV> void RIGID_BODY_INTERSECTIONS<TV>::
Append_All_Intersections(const int id_1,const int id_2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T contour_value) const
{
    RIGID_BODY<TV> &body0=rigid_body_collection.Rigid_Body(id_1),&body1=rigid_body_collection.Rigid_Body(id_2);
    PARTICLES_IN_IMPLICIT_OBJECT<TV>::Append_All_Intersections(body0,body1,particle_intersections,contour_value,use_triangle_hierarchy,use_edge_intersection,use_triangle_hierarchy_center_phi_test);
}
//#####################################################################
// Function Particle_Partition
//#####################################################################
template<class TV> const PARTICLE_PARTITION<TV>& RIGID_BODY_INTERSECTIONS<TV>::
Particle_Partition(const RIGID_BODY<TV>& rigid_body) const
{
    // TODO: regenerate if particle_partition_size changes
    PHYSBAM_ASSERT(rigid_body.simplicial_object);
    if(!rigid_body.simplicial_object->particle_partition)
        rigid_body.simplicial_object->Initialize_Particle_Partition(particle_partition_size);
    return *rigid_body.simplicial_object->particle_partition;
}
//#####################################################################
// Function Particle_Hierarchy
//#####################################################################
template<class TV> const PARTICLE_HIERARCHY<TV>& RIGID_BODY_INTERSECTIONS<TV>::
Particle_Hierarchy(const RIGID_BODY<TV>& rigid_body) const
{
    PHYSBAM_ASSERT(rigid_body.simplicial_object);
    const GEOMETRY_PARTICLES<TV>& particles=rigid_body.simplicial_object->particles;
    try{
        return *particle_hierarchies.Get(&particles);}
    catch(KEY_ERROR&){
        PARTICLE_HIERARCHY<TV>* particle_hierarchy(new PARTICLE_HIERARCHY<TV>(particles.X));
        particle_hierarchy->Update_Box_Radii();
        particle_hierarchies.Insert(&particles,particle_hierarchy);
        return *particle_hierarchy;}
}
//#####################################################################
// Function Particles_In_Levelset
//#####################################################################
template<class TV> void RIGID_BODY_INTERSECTIONS<TV>::
Particles_In_Levelset(const int particle_body_id,const int levelset_body_id,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,
    const T contour_value,const bool exit_early) const
{
    RIGID_BODY<TV> &body0=rigid_body_collection.Rigid_Body(particle_body_id),&body1=rigid_body_collection.Rigid_Body(levelset_body_id);
    if(use_particle_hierarchy){
        PARTICLES_IN_IMPLICIT_OBJECT<TV>::Particles_In_Implicit_Object_Hierarchy(body0,body1,particle_intersections,contour_value,particle_hierarchies);}
    else if(!use_particle_partition){
        PARTICLES_IN_IMPLICIT_OBJECT<TV>::Particles_In_Implicit_Object(body0,body1,particle_intersections,contour_value,exit_early);}
    else{ // use particle partition
        PARTICLES_IN_IMPLICIT_OBJECT<TV>::Particles_In_Implicit_Object_Partition(body0,body1,particle_intersections,contour_value,use_particle_partition_center_phi_test,particle_partition_size,exit_early);
    }
}
//#####################################################################
// Function Oriented_Box2_In_Body1_Coordinates
//#####################################################################
template<class TV> typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX RIGID_BODY_INTERSECTIONS<TV>::
Oriented_Box2_In_Body1_Coordinates(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1) const
{
    return T_ORIENTED_BOX(body1.Object_Space_Bounding_Box(),body0.Frame().Inverse_Times(body1.Frame()));
}
//#####################################################################
// Function Initialize_Transformation_From_Body1_To_Body2_Coordinates
//#####################################################################
template<class TV> void RIGID_BODY_INTERSECTIONS<TV>::
Initialize_Transformation_From_Body1_To_Body2_Coordinates(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1) const
{
    FRAME<TV> frame=body1.Frame().Inverse_Times(body0.Frame());
    rotation=frame.r.Rotation_Matrix();translation=frame.t;
    rotation_reverse=rotation.Transposed();translation_reverse=-(rotation_reverse*translation); // reverse is from body1 to body0 coordinates
}
//#####################################################################
namespace PhysBAM{
template class RIGID_BODY_INTERSECTIONS<VECTOR<float,1> >;
template class RIGID_BODY_INTERSECTIONS<VECTOR<float,2> >;
template class RIGID_BODY_INTERSECTIONS<VECTOR<float,3> >;
template class RIGID_BODY_INTERSECTIONS<VECTOR<double,1> >;
template class RIGID_BODY_INTERSECTIONS<VECTOR<double,2> >;
template class RIGID_BODY_INTERSECTIONS<VECTOR<double,3> >;
}
