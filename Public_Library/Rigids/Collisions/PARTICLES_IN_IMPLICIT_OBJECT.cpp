//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace PARTICLES_IN_IMPLICIT_OBJECT
//##################################################################### 
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
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
#include <Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>

namespace PhysBAM{

//#####################################################################
// Function Append_All_Intersections
//#####################################################################
template<class TV> void PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Append_All_Intersections(RIGID_BODY<TV>& body0, RIGID_BODY<TV>& body1,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_triangle_hierarchy,const bool use_edge_intersection,const bool use_triangle_hierarchy_center_phi_test)
{
    if(!use_triangle_hierarchy){  
        Append_All_Intersections_Points(body0,body1,particle_intersections,contour_value);}
    else{ // using triangle hierarchy
        if(!use_edge_intersection){ 
            Append_All_Intersections_Triangles(body0,body1,particle_intersections,contour_value,use_triangle_hierarchy_center_phi_test);}
        else{ // use edge intersections
            Append_All_Intersections_Edges(body0,body1,particle_intersections,contour_value,use_triangle_hierarchy_center_phi_test);}
    }
}
template<class TV> void PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Append_All_Intersections_Points(RIGID_BODY<TV>& body0, RIGID_BODY<TV>& body1,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value)
{
    if(body0.simplicial_object && body1.implicit_object){
        Particles_In_Implicit_Object(body0,body1,particle_intersections,contour_value,false);}
    if(body0.implicit_object && body1.simplicial_object){
        Particles_In_Implicit_Object(body1,body0,particle_intersections,contour_value,false);}
}
//#####################################################################
// Function Append_All_Intersections
//#####################################################################
template<class TV> void PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Append_All_Intersections_Triangles(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_triangle_hierarchy_center_phi_test)
{
    typedef typename TV::SCALAR T;
    static const int d=TV::m;

    FRAME<TV> frame=body1.Frame().Inverse_Times(body0.Frame());
    MATRIX<T,d> rotation=frame.r.Rotation_Matrix();
    VECTOR<T,d> translation=frame.t;
    MATRIX<T,d> rotation_reverse=rotation.Transposed();
    VECTOR<T,d> translation_reverse=-(rotation_reverse*translation);

    //int id_1=body0.particle_index,id_2=body1.particle_index;
    ARRAY<int> triangle_list;triangle_list.Preallocate(50);
    if(body0.simplicial_object && body1.implicit_object){
        Get_Interfering_Simplices(body0,body1,triangle_list,rotation,translation,use_triangle_hierarchy_center_phi_test);
        Intersections_Using_Hierarchy(body0,body1,triangle_list,particle_intersections,contour_value,false,rotation,translation);}
    if(body0.implicit_object && body1.simplicial_object){
        Get_Interfering_Simplices(body1,body0,triangle_list,rotation_reverse,translation_reverse,use_triangle_hierarchy_center_phi_test);
        Intersections_Using_Hierarchy(body1,body0,triangle_list,particle_intersections,contour_value,false,rotation_reverse,translation_reverse);}
}
//#####################################################################
// Function Append_All_Intersections
//#####################################################################
template<class TV> void PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Append_All_Intersections_Edges(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_triangle_hierarchy_center_phi_test)
{
    typedef typename TV::SCALAR T;
    static const int d=TV::m;

    FRAME<TV> frame=body1.Frame().Inverse_Times(body0.Frame());
    MATRIX<T,d> rotation=frame.r.Rotation_Matrix();
    VECTOR<T,d> translation=frame.t;
    MATRIX<T,d> rotation_reverse=rotation.Transposed();
    VECTOR<T,d> translation_reverse=-(rotation_reverse*translation);

    //int id_1=body0.particle_index,id_2=body1.particle_index;
    ARRAY<int> triangle_list1,triangle_list2;triangle_list1.Preallocate(50);triangle_list2.Preallocate(50);
    if(body0.simplicial_object)
        Get_Interfering_Simplices(body0,body1,triangle_list1,rotation,translation,use_triangle_hierarchy_center_phi_test);
    if(body1.simplicial_object){
        Get_Interfering_Simplices(body1,body0,triangle_list2,rotation_reverse,translation_reverse,use_triangle_hierarchy_center_phi_test);
        if(body0.implicit_object) Intersections_Using_Hierarchy_And_Edges(body1,body0,triangle_list2,triangle_list1,particle_intersections,contour_value,false,rotation_reverse,translation_reverse);}
    if(body0.simplicial_object && body1.implicit_object){            
        Intersections_Using_Hierarchy_And_Edges(body0,body1,triangle_list1,triangle_list2,particle_intersections,contour_value,false,rotation,translation);}
}
//#####################################################################
// Function Simplex_Hierarchy
//#####################################################################
template<class TV> const typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::HIERARCHY& PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Simplex_Hierarchy(const RIGID_BODY<TV>& rigid_body)
{
    PHYSBAM_ASSERT(rigid_body.simplicial_object);
    if(!rigid_body.simplicial_object->hierarchy){
        rigid_body.simplicial_object->Initialize_Hierarchy();
        rigid_body.simplicial_object->hierarchy->Update_Box_Radii();}
    return *rigid_body.simplicial_object->hierarchy;
}
//#####################################################################
// Function Get_Interfering_Simplices
//#####################################################################
template<class TV> void PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Get_Interfering_Simplices(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1,ARRAY<int>& simplex_list,MATRIX<typename TV::SCALAR,TV::m,TV::m>& rotation,TV& translation,const bool use_triangle_hierarchy_center_phi_test)
{
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;

    simplex_list.Remove_All();
    const typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::HIERARCHY& hierarchy=Simplex_Hierarchy(body0);
    if(use_triangle_hierarchy_center_phi_test && body1.implicit_object)
        hierarchy.Intersection_List(*body1.implicit_object->object_space_implicit_object,rotation,translation,simplex_list);
    else
    {
        FRAME<TV> frame=body0.Frame().Inverse_Times(body0.Frame());
        hierarchy.Intersection_List(T_ORIENTED_BOX(body0.Object_Space_Bounding_Box(),body0.Frame().Inverse_Times(body0.Frame())),simplex_list);
    }
}
//#####################################################################
// Function Transform_From_Body1_To_Body2_Coordinates
//#####################################################################
template<class TV> TV
Transform_From_Body1_To_Body2_Coordinates(TV& v,MATRIX<typename TV::SCALAR,TV::m,TV::m>& rotation,TV& translation)
{return rotation*v+translation;}
//#####################################################################
// Function Intersections_Using_Hierarchy
//#####################################################################
template<class TV> void  PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Intersections_Using_Hierarchy(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<int>& simplex_list,
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool exit_early,MATRIX<typename TV::SCALAR,TV::m,TV::m>& rotation,TV& translation)
{
    ARRAY<bool> checked(particle_body.simplicial_object->particles.Size());
    ARRAY_VIEW<bool>* collidable=particle_body.simplicial_object->particles.template Get_Array<bool>("collidable");
    for(int t=0;t<simplex_list.m;t++){const VECTOR<int,TV::m>& nodes=particle_body.simplicial_object->mesh.elements(simplex_list(t));
        for(int i=0;i<nodes.m;i++){
            if(!checked(nodes[i])){checked(nodes[i])=true;
                if((!collidable || (*collidable)(nodes[i])) && levelset_body.implicit_object->object_space_implicit_object->Lazy_Inside(Transform_From_Body1_To_Body2_Coordinates(particle_body.simplicial_object->particles.X(nodes[i]),rotation,translation))){
                    particle_intersections.Append(RIGID_BODY_PARTICLE_INTERSECTION<TV>(particle_body.simplicial_object->particles.X(nodes[i]),nodes[i],particle_body.particle_index,levelset_body.particle_index));
                    if(exit_early) return;}}}}
}
//#####################################################################
// Function Intersections_Using_Hierarchy_And_Edges
//#####################################################################
// body1 doesn't require a triangulated surface, but if it has one it will be used for edge-face intersections
template<class T> void
Intersections_Using_Hierarchy_And_Edges_Helper(RIGID_BODY<VECTOR<T,1> >& body0,RIGID_BODY<VECTOR<T,1> >& body1,ARRAY<int>& simplex_list1,ARRAY<int>& simplex_list2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<VECTOR<T,1> > >& particle_intersections,const T contour_value,const bool exit_early,MATRIX<T,1>& rotation,VECTOR<T,1>& translation)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Intersections_Using_Hierarchy_And_Edges_Helper
//#####################################################################
template<class T> void
Intersections_Using_Hierarchy_And_Edges_Helper(RIGID_BODY<VECTOR<T,2> >& body0,RIGID_BODY<VECTOR<T,2> >& body1,ARRAY<int>& simplex_list1,ARRAY<int>& simplex_list2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<VECTOR<T,2> > >& particle_intersections,const T contour_value,const bool exit_early,MATRIX<T,2>& rotation,VECTOR<T,2>& translation)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Intersections_Using_Hierarchy_And_Edges_Helper
//#####################################################################
template<class T> void
Intersections_Using_Hierarchy_And_Edges_Helper(RIGID_BODY<VECTOR<T,3> >& body0,RIGID_BODY<VECTOR<T,3> >& body1,ARRAY<int>& triangle_list1,ARRAY<int>& triangle_list2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<VECTOR<T,3> > >& particle_intersections,const T contour_value,const bool exit_early,MATRIX<T,3>& rotation,VECTOR<T,3>& translation)
{
    typedef VECTOR<T,3> TV;

    int id1=body0.particle_index,id2=body1.particle_index;

    ARRAY<bool> checked(body0.simplicial_object->particles.Size()),segment_checked(body0.simplicial_object->mesh.segment_mesh->elements.m);
    ARRAY<T> phi_value(body0.simplicial_object->particles.Size()); // only meaningful if node is outside
    ARRAY_VIEW<bool>* collidable=body0.simplicial_object->particles.template Get_Array<bool>("collidable");
    T value;TRIANGLE_MESH& mesh=body0.simplicial_object->mesh;
    if(!mesh.segment_mesh) mesh.Initialize_Segment_Mesh();
    if(!mesh.segment_mesh->incident_elements) mesh.segment_mesh->Initialize_Incident_Elements();
    if(!mesh.element_edges) mesh.Initialize_Element_Edges();
    if(!body0.simplicial_object->segment_lengths) body0.simplicial_object->Initialize_Segment_Lengths();
    if(!body0.simplicial_object->triangle_list) body0.simplicial_object->Update_Triangle_List();
    for(int t=0;t<triangle_list1.m;t++) for(int i=0;i<3;i++){
        int edge=(*mesh.element_edges)(triangle_list1(t))(i);
        if(!segment_checked(edge)){
            int node1=mesh.segment_mesh->elements(edge)(0);
            if(!checked(node1) && (!collidable || (*collidable)(node1))){
                if(!body1.implicit_object->object_space_implicit_object->Lazy_Outside_Extended_Levelset_And_Value(
                       Transform_From_Body1_To_Body2_Coordinates(body0.simplicial_object->particles.X(node1),rotation,translation),value,contour_value)){
                    particle_intersections.Append(RIGID_BODY_PARTICLE_INTERSECTION<TV>(body0.simplicial_object->particles.X(node1),node1,id1,id2));
                    if(exit_early) return;
                    for(int j=0;j<(*mesh.segment_mesh->incident_elements)(node1).m;j++) // mark incident edges as checked
                        segment_checked((*mesh.segment_mesh->incident_elements)(node1)(j))=true;}
                else phi_value(node1)=value;}
            int node2=mesh.segment_mesh->elements(edge)(1);
            if(!checked(node2) && (!collidable || (*collidable)(node1))){
                if(!body1.implicit_object->object_space_implicit_object->Lazy_Outside_Extended_Levelset_And_Value(
                       Transform_From_Body1_To_Body2_Coordinates(body0.simplicial_object->particles.X(node2),rotation,translation),value,contour_value)){
                    particle_intersections.Append(RIGID_BODY_PARTICLE_INTERSECTION<TV>(body0.simplicial_object->particles.X(node2),node2,id1,id2));
                    if(exit_early) return;
                    for(int j=0;j<(*mesh.segment_mesh->incident_elements)(node2).m;j++) // mark incident edges as checked
                        segment_checked((*mesh.segment_mesh->incident_elements)(node2)(j))=true;}
                else phi_value(node2)=value;}
            if(collidable && !((*collidable)(node1) && (*collidable)(node2))) {checked(node1)=true;checked(node2)=true;segment_checked(edge)=true;continue;}
            if(phi_value(node1) > 0 && phi_value(node2) > 0 && phi_value(node1)+phi_value(node2) <= (*body0.simplicial_object->segment_lengths)(edge) && body1.simplicial_object){
                TV p1=body0.simplicial_object->particles.X(node1),p2=body0.simplicial_object->particles.X(node2),
                    x0=Transform_From_Body1_To_Body2_Coordinates(p1,rotation,translation),x1=Transform_From_Body1_To_Body2_Coordinates(p2,rotation,translation);
                RAY<VECTOR<T,3> > ray(SEGMENT_3D<T>(x0,x1));T t_max=ray.t_max,one_over_t_max=1/t_max;
                ARRAY<PAIR<T,bool> > intersections;
                intersections.Preallocate(10); // pair is (t_intersect, going_in)
                for(int k=0;k<triangle_list2.m;k++){ray.t_max=t_max;
                    if(INTERSECTION::Lazy_Intersects(ray,(*body1.simplicial_object->triangle_list)(triangle_list2(k)))){
                        bool going_in=TV::Dot_Product((*body1.simplicial_object->triangle_list)(triangle_list2(k)).Normal(),ray.direction) < 0;
                        intersections.Append(PAIR<T,bool>(ray.t_max,going_in));}}
                intersections.Sort([](const auto& a,const auto& b){return a.x<b.x;});
                for(int kk=1;kk<intersections.m;kk++) // find midpoints of intersecting edges
                    if(intersections(kk).y && !intersections(kk+1).y){              
                        T s=(T).5*(intersections(kk).x+intersections(kk+1).x)*one_over_t_max;
                        // TODO: FIX THIS -- improve interface so we can return an interpolation point on segment instead of just node1
                        particle_intersections.Append(RIGID_BODY_PARTICLE_INTERSECTION<TV>((1-s)*p1+s*p2,node1,id1,id2));
                        if(exit_early) return;
                        kk++;}} // skip an extra one here
            checked(node1)=true;checked(node2)=true;segment_checked(edge)=true;}}
}
//#####################################################################
// Function Intersections_Using_Hierarchy_And_Edges
//#####################################################################
template<class TV> void PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Intersections_Using_Hierarchy_And_Edges(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,ARRAY<int>& simplex_list1,ARRAY<int>& simplex_list2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool exit_early,MATRIX<typename TV::SCALAR,TV::m,TV::m>& rotation,TV& translation)
{
    Intersections_Using_Hierarchy_And_Edges_Helper(body0,body1,simplex_list1,simplex_list2,particle_intersections,contour_value,exit_early,rotation,translation);
}
//#####################################################################
// Function Particles_In_Implicit_Object
//#####################################################################
template<class TV> void PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Particles_In_Implicit_Object(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool exit_early){
    typedef typename TV::SCALAR T;
    static const int d=TV::m;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;

    ARRAY_VIEW<bool>* collidable=particle_body.simplicial_object->particles.template Get_Array<bool>("collidable");
    FRAME<TV> frame=levelset_body.Frame().Inverse_Times(particle_body.Frame());
    MATRIX<T,d> rotation=frame.r.Rotation_Matrix();
    VECTOR<T,d> translation=frame.t;
    RANGE<TV> bounding_box2_in_body1_coordinates=T_ORIENTED_BOX(levelset_body.Object_Space_Bounding_Box(),particle_body.Frame().Inverse_Times(levelset_body.Frame())).Axis_Aligned_Bounding_Box();

    ARRAY_VIEW<TV>& particles_X=particle_body.simplicial_object->particles.X;
    IMPLICIT_OBJECT<VECTOR<T,d> >& object_space_implicit_object=*levelset_body.implicit_object->object_space_implicit_object;
    for(int p=0;p<particles_X.Size();p++)
        if((!collidable || (*collidable)(p)) && bounding_box2_in_body1_coordinates.Lazy_Inside(particles_X(p)) && 
            object_space_implicit_object.Lazy_Inside(Transform_From_Body1_To_Body2_Coordinates(particles_X(p),rotation,translation),contour_value)){
            particle_intersections.Append(RIGID_BODY_PARTICLE_INTERSECTION<TV>(particles_X(p),p,particle_body.particle_index,levelset_body.particle_index));
            if(exit_early) return;}
}
//#####################################################################
// Function Particle_Hierarchy
//#####################################################################
template<class TV> const PARTICLE_HIERARCHY<TV>& 
Particle_Hierarchy(const RIGID_BODY<TV>& rigid_body,HASHTABLE<const GEOMETRY_PARTICLES<TV>*,const PARTICLE_HIERARCHY<TV>*>& particle_hierarchies)
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
// Function Particles_In_Implicit_Object_Hierarchy
//#####################################################################
template<class TV> void PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Particles_In_Implicit_Object_Hierarchy(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,HASHTABLE<const GEOMETRY_PARTICLES<TV>*,const PARTICLE_HIERARCHY<TV>*>& particle_hierarchies)
{
    typedef typename TV::SCALAR T;
    static const int d=TV::m;

    FRAME<TV> frame=levelset_body.Frame().Inverse_Times(particle_body.Frame());
    MATRIX<T,d> rotation=frame.r.Rotation_Matrix();
    VECTOR<T,d> translation=frame.t;

    IMPLICIT_OBJECT<TV>& object_space_implicit_object=*levelset_body.implicit_object->object_space_implicit_object;
    const PARTICLE_HIERARCHY<TV>& particle_hierarchy=Particle_Hierarchy(particle_body,particle_hierarchies);
    PARTICLE_LEVELSET_VISITOR<TV> visitor(particle_hierarchy,object_space_implicit_object,rotation,translation,contour_value,
        particle_intersections,particle_body.particle_index,levelset_body.particle_index);
    particle_hierarchy.Intersection_List(visitor);
}
//#####################################################################
// Function Particle_Partition
//#####################################################################
template<class TV> const PARTICLE_PARTITION<TV>&
Particle_Partition(const RIGID_BODY<TV>& rigid_body,const VECTOR<int,TV::m>& particle_partition_size)
{
    // TODO: regenerate if particle_partition_size changes
    PHYSBAM_ASSERT(rigid_body.simplicial_object);
    if(!rigid_body.simplicial_object->particle_partition)
        rigid_body.simplicial_object->Initialize_Particle_Partition(particle_partition_size);
    return *rigid_body.simplicial_object->particle_partition;
}
//#####################################################################
// Function Particles_In_Implicit_Object_Partition
//#####################################################################
template<class TV> void PARTICLES_IN_IMPLICIT_OBJECT<TV>::
Particles_In_Implicit_Object_Partition(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_particle_partition_center_phi_test,const VECTOR<int,TV::m>& particle_partition_size,const bool exit_early)
{
    typedef typename TV::SCALAR T;
    const int d=TV::m;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
    
    FRAME<TV> frame=levelset_body.Frame().Inverse_Times(particle_body.Frame());
    MATRIX<T,d> rotation=frame.r.Rotation_Matrix();
    VECTOR<T,d> translation=frame.t;
    RANGE<TV> bounding_box2_in_body1_coordinates=T_ORIENTED_BOX(levelset_body.Object_Space_Bounding_Box(),particle_body.Frame().Inverse_Times(levelset_body.Frame())).Axis_Aligned_Bounding_Box();

    IMPLICIT_OBJECT<VECTOR<T,d> >& object_space_implicit_object=*levelset_body.implicit_object->object_space_implicit_object;
    ARRAY_VIEW<TV>& particles_X=particle_body.simplicial_object->particles.X;
    const PARTICLE_PARTITION<TV>& particle_partition=Particle_Partition(particle_body,particle_partition_size);
    if(!use_particle_partition_center_phi_test){
        RANGE<VECTOR<int,d> > range=particle_partition.Range(bounding_box2_in_body1_coordinates);      
        for(CELL_ITERATOR<TV> iterator(particle_partition.grid,range);iterator.Valid();iterator.Next()){
            const ARRAY<int>& particles_in_cell=particle_partition.partition(iterator.Cell_Index());
            for(int t=0;t<particles_in_cell.m;t++){int p=particles_in_cell(t);
                if(bounding_box2_in_body1_coordinates.Lazy_Inside(particles_X(p)) &&
                    object_space_implicit_object.Lazy_Inside(Transform_From_Body1_To_Body2_Coordinates(particles_X(p),rotation,translation),contour_value)){
                    particle_intersections.Append(RIGID_BODY_PARTICLE_INTERSECTION<TV>(particles_X(p),p,particle_body.particle_index,levelset_body.particle_index));
                    if(exit_early) return;}}}}
    else{ // use particle partition center phi test
        ARRAY<VECTOR<int,d> > intersection_list;
        particle_partition.Intersection_List(object_space_implicit_object,rotation,translation,intersection_list,contour_value);
        for(int i=0;i<intersection_list.m;i++){
            const ARRAY<int>& particles_in_cell=particle_partition.partition(intersection_list(i));
            for(int t=0;t<particles_in_cell.m;t++){int p=particles_in_cell(t);
                if(bounding_box2_in_body1_coordinates.Lazy_Inside(particles_X(p)) && 
                    object_space_implicit_object.Lazy_Inside(Transform_From_Body1_To_Body2_Coordinates(particles_X(p),rotation,translation),contour_value)){
                    particle_intersections.Append(RIGID_BODY_PARTICLE_INTERSECTION<TV>(particles_X(p),p,particle_body.particle_index,levelset_body.particle_index));
                    if(exit_early) return;}}}}
}
//#####################################################################
template struct PARTICLES_IN_IMPLICIT_OBJECT<VECTOR<float,1> >;
template struct PARTICLES_IN_IMPLICIT_OBJECT<VECTOR<float,2> >;
template struct PARTICLES_IN_IMPLICIT_OBJECT<VECTOR<float,3> >;
template struct PARTICLES_IN_IMPLICIT_OBJECT<VECTOR<double,1> >;
template struct PARTICLES_IN_IMPLICIT_OBJECT<VECTOR<double,2> >;
template struct PARTICLES_IN_IMPLICIT_OBJECT<VECTOR<double,3> >;
}
