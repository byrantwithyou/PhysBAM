//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/KD_TREE.h>
#include <Core/Math_Tools/RANGE.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T_ARRAY> PARTICLE_HIERARCHY<TV,T_ARRAY>::
    PARTICLE_HIERARCHY(const T_ARRAY& X_input,const bool update_boxes,const int particles_per_group_input)
    :X(X_input),particles_per_group(particles_per_group_input)
{
    if(X.Size()){Initialize_Hierarchy_Using_KD_Tree();if(update_boxes) Update_Boxes();}else{leaves=0;root=0;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T_ARRAY> PARTICLE_HIERARCHY<TV,T_ARRAY>::
    ~PARTICLE_HIERARCHY()
{
}
//#####################################################################
// Function Update_Boxes
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Update_Boxes(const T extra_thickness)
{
    Update_Leaf_Boxes(extra_thickness);Update_Nonleaf_Boxes();
}
//#####################################################################
// Function Update_Boxes
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Update_Boxes(const ARRAY<TV>& X,const T extra_thickness)
{
    Update_Leaf_Boxes(X,extra_thickness);Update_Nonleaf_Boxes();
}
//#####################################################################
// Function Update_Leaf_Boxes
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Update_Leaf_Boxes(const T extra_thickness)
{
    Calculate_Bounding_Boxes(box_hierarchy);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Intersection_List(const TV& point,ARRAY<int>& intersection_list,const T thickness_over_two) const
{
    if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,point,group_list,thickness_over_two);
        for(int i=0;i<group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,point,intersection_list,thickness_over_two);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Intersection_List(const RANGE<TV>& test_box,ARRAY<int>& intersection_list,const T thickness_over_two) const
{
    if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,test_box,group_list,thickness_over_two);
        for(int i=0;i<group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,test_box,intersection_list,thickness_over_two);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Intersection_List(const ORIENTED_BOX<TV>& test_box,ARRAY<int>& intersection_list) const
{
    if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,test_box,group_list);
        for(int i=0;i<group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,test_box,intersection_list);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Intersection_List(const T_HYPERPLANE& test_plane,ARRAY<int>& intersection_list,const T thickness_over_two) const
{
    if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,test_plane,group_list,thickness_over_two);
        for(int i=0;i<group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,test_plane,intersection_list,thickness_over_two);
}
//#####################################################################
// Function Intersection_List
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Intersection_List(const IMPLICIT_OBJECT<TV>& implicit_object,const MATRIX<T,TV::dimension>& rotation,const TV& translation,ARRAY<int>& intersection_list,const T contour_value) const
{
    if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,implicit_object,rotation,translation,group_list,contour_value);
        for(int i=0;i<group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,implicit_object,rotation,translation,intersection_list,contour_value);
}
//#####################################################################
// Function Initialize_Hierarchy_Using_KD_Tree
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Initialize_Hierarchy_Using_KD_Tree()
{
    particles_in_group.Clean_Memory();
    KD_TREE<TV> kd_tree(false);
    ARRAY<TV> X_copy(X);
    if(particles_per_group){
        kd_tree.Create_Left_Balanced_KD_Tree_With_Grouping(X_copy,particles_in_group,particles_per_group);
        leaves=particles_in_group.m;}
    else{
        kd_tree.Create_Left_Balanced_KD_Tree(X_copy);
        leaves=X.Size();}
    parents.Resize(leaves);
    children.Remove_All();root=Initialize_Hierarchy_Using_KD_Tree_Helper(kd_tree.root_node);
    assert(root==2*leaves-2);box_hierarchy.Resize(root+1);box_radius.Resize(root+1);
}
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes) const
{
    if(particles_per_group) for(int k=0;k<leaves;k++){
        bounding_boxes(k)=RANGE<TV>::Bounding_Box(X.Subset(particles_in_group(k)));}
    else for(int k=0;k<leaves;k++) bounding_boxes(k)=RANGE<TV>(X(k),X(k));
}
//#####################################################################
// Function Calculate_Bounding_Box_Radii
//#####################################################################
// at the boxes center, but may be a tighter bound on the triangle than the box
template<class TV,class T_ARRAY> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius)
{
    if(particles_per_group) for(int k=0;k<leaves;k++){
        TV center=bounding_boxes(k).Center();T max_radius_squared=0;
        for(int i=0;i<particles_in_group(k).m;i++) max_radius_squared=max(max_radius_squared,(X(particles_in_group(k)(i))-center).Magnitude_Squared());
        radius(k)=sqrt(max_radius_squared);}
    else for(int k=0;k<leaves;k++){TV center=bounding_boxes(k).Center();radius(k)=sqrt((X(k)-center).Magnitude_Squared());}
}
//#####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(T,d) \
    template class PARTICLE_HIERARCHY<VECTOR<T,d> >; \
    template class PARTICLE_HIERARCHY<VECTOR<T,d>,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<T,d>,int> > >;
INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
}
