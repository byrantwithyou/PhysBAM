//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_DEFORMABLE_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_POINT_FACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_POINT_FACE<TV>::
COMBINED_COLLISIONS_POINT_FACE(TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,bool prune_pairs_input)
    :BASE(triangle_repulsions_input,prune_pairs_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_POINT_FACE<TV>::
~COMBINED_COLLISIONS_POINT_FACE()
{
}
//#####################################################################
// Function Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_POINT_FACE<TV>::
Discover(const T dt,const T time)
{
    collisions.Remove_All();
    if(!prune_pairs) Discover_Saved_Pairs();
    else Discover_Pruned();
}
template<class T,class TV>
bool Point_Face_Interaction(POINT_SIMPLEX_1D<T>& face,const TV& X,T& thickness,T& distance)
{
    return false;
}
template<class T_FACE,class T,class TV>
bool Point_Face_Interaction(T_FACE& face,const TV& X,T& thickness,T& distance)
{
    return face.Point_Face_Interaction(X,thickness,false,distance);
}
template<class T,class TV>
void Point_Face_Interaction_Data(POINT_SIMPLEX_1D<T>& face,const TV& X,T& distance,TV& normal,TV& weights,bool exit_early)
{
}
template<class T_FACE,class T,class TV>
void Point_Face_Interaction_Data(T_FACE& face,const TV& X,T& distance,TV& normal,TV& weights,bool exit_early)
{
    face.Point_Face_Interaction_Data(X,distance,normal,weights,exit_early);
}
//#####################################################################
// Function Discover_Pruned
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_POINT_FACE<TV>::
Discover_Pruned()
{
    ARRAY_VIEW<const TV> X=triangle_repulsions.geometry.deformable_body_collection.particles.X;
    ARRAY_VIEW<const TV> V=triangle_repulsions.geometry.deformable_body_collection.particles.V;
    for(int i=1;i<=triangle_repulsions.point_face_interaction_pairs.m;i++){
        POINT_FACE_REPULSION_PAIR<TV>& pf=triangle_repulsions.point_face_interaction_pairs(i);

        typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX face(X.Subset(pf.nodes.Remove_Index(1)));
        T thickness=triangle_repulsions.repulsion_thickness_detection_multiplier*pf.Total_Repulsion_Thickness(triangle_repulsions.repulsion_thickness);
        if(!Point_Face_Interaction(face,X(pf.nodes(1)),thickness,pf.distance)) continue;

        typename BASE::COLLISION col;
        col.nodes=pf.nodes;
        Point_Face_Interaction_Data(face,X(pf.nodes[1]),pf.distance,col.normal,pf.weights,false);
        col.weights=-pf.weights.Insert(-1,1);
        col.relative_velocity=ARRAYS_COMPUTATIONS::Weighted_Sum(col.weights,V.Subset(col.nodes));

        if(TV::Dot_Product(col.relative_velocity,col.normal)>=0) continue;

        collisions.Append(col);}
}
//#####################################################################
// Function Discover_Saved_Pairs
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_POINT_FACE<TV>::
Discover_Saved_Pairs()
{
    ARRAY_VIEW<const TV> V=triangle_repulsions.geometry.deformable_body_collection.particles.V;
    for(int i=1;i<=triangle_repulsions.internal_point_face_precomputed.m;i++){
        const PRECOMPUTE_PROJECT_POINT_FACE<TV>& pf=triangle_repulsions.internal_point_face_precomputed(i);
        typename BASE::COLLISION col;
        col.weights=-pf.weights.Insert(-1,1);
        col.normal=pf.normal;
        col.nodes=pf.nodes;
        col.relative_velocity=ARRAYS_COMPUTATIONS::Weighted_Sum(col.weights,V.Subset(col.nodes));
        collisions.Append(col);}
}
template class COMBINED_COLLISIONS_POINT_FACE<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_POINT_FACE<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_POINT_FACE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_POINT_FACE<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_POINT_FACE<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_POINT_FACE<VECTOR<double,3> >;
#endif
