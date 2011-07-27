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
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_EDGE_EDGE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_EDGE_EDGE<TV>::
COMBINED_COLLISIONS_EDGE_EDGE(TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,bool prune_pairs_input)
    :BASE(triangle_repulsions_input,prune_pairs_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_EDGE_EDGE<TV>::
~COMBINED_COLLISIONS_EDGE_EDGE()
{
}
//#####################################################################
// Function Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_EDGE_EDGE<TV>::
Discover(const T dt,const T time)
{
    collisions.Remove_All();
    if(!prune_pairs) Discover_Saved_Pairs();
    else Discover_Pruned();
}
template<class T> static void Edge_Edge_Weight(const VECTOR<T,1>& weights,typename COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<T,1> >::COLLISION& col){}
template<class T> static void Edge_Edge_Weight(const VECTOR<T,2>& weights,typename COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<T,2> >::COLLISION& col){col.weights=VECTOR<T,2>(1,-1);}
template<class T> static void Edge_Edge_Weight(const VECTOR<T,2>& weights,typename COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<T,3> >::COLLISION& col)
{col.weights=VECTOR<T,4>(1-weights(1),weights(1),-(1-weights(2)),-weights(2));}
//#####################################################################
// Function Update_Repulsion_Pairs_Using_History
//#####################################################################
template<class T,class T_WEIGHT> static bool Edge_Edge_Interaction(ARRAY_VIEW<const VECTOR<T,1> > X,ARRAY_VIEW<const VECTOR<T,1> > V,T_WEIGHT& weights,
    typename COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<T,1> >::COLLISION& col,T total_repulsion_thickness,T small_number)
{
    return false;
}
template<class T,class T_WEIGHT> static bool Edge_Edge_Interaction(ARRAY_VIEW<const VECTOR<T,2> > X,ARRAY_VIEW<const VECTOR<T,2> > V,T_WEIGHT& weights,
    typename COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<T,2> >::COLLISION& col,T total_repulsion_thickness,T small_number)
{
    col.normal=X(col.nodes(1))-X(col.nodes(2));
    return col.normal.Normalize()<=total_repulsion_thickness;
}
template<class T,class T_WEIGHT> static bool Edge_Edge_Interaction(ARRAY_VIEW<const VECTOR<T,3> > X,ARRAY_VIEW<const VECTOR<T,3> > V,T_WEIGHT& weights,
    typename COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<T,3> >::COLLISION& col,T total_repulsion_thickness,T small_number)
{
    SEGMENT_3D<T> segment1(X.Subset(VECTOR<int,2>(col.nodes(1),col.nodes(2)))),segment2(X.Subset(VECTOR<int,2>(col.nodes(3),col.nodes(4))));
    T distance;
    bool ret=segment1.Edge_Edge_Interaction(segment2,total_repulsion_thickness,distance,col.normal,weights,false);
    segment1.Edge_Edge_Interaction_Data(segment2,V.Subset(col.nodes),distance,col.normal,weights,small_number);
    return ret;
}
//#####################################################################
// Function Discover_Pruned
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_EDGE_EDGE<TV>::
Discover_Pruned()
{
    ARRAY_VIEW<const TV> X=triangle_repulsions.geometry.deformable_body_collection.particles.X;
    ARRAY_VIEW<const TV> V=triangle_repulsions.geometry.deformable_body_collection.particles.V;
    for(int i=1;i<=triangle_repulsions.edge_edge_interaction_pairs.m;i++){
        EDGE_EDGE_REPULSION_PAIR<TV>& ee=triangle_repulsions.edge_edge_interaction_pairs(i);
        typename BASE::COLLISION col;
        col.nodes=ee.nodes;

        T thickness=triangle_repulsions.repulsion_thickness_detection_multiplier*ee.Total_Repulsion_Thickness(triangle_repulsions.repulsion_thickness);
        if(!Edge_Edge_Interaction(X,V,ee.weights,col,thickness,triangle_repulsions.geometry.small_number)) continue;
        Edge_Edge_Weight(ee.weights,col);

        col.relative_velocity=ARRAYS_COMPUTATIONS::Weighted_Sum(col.weights,V.Subset(col.nodes));

        if(TV::Dot_Product(col.relative_velocity,col.normal)>=0) continue;

        collisions.Append(col);}
}
//#####################################################################
// Function Discover_Saved_Pairs
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_EDGE_EDGE<TV>::
Discover_Saved_Pairs()
{
    ARRAY_VIEW<const TV> X=triangle_repulsions.geometry.deformable_body_collection.particles.X;
    ARRAY_VIEW<const TV> V=triangle_repulsions.geometry.deformable_body_collection.particles.V;
    for(int i=1;i<=triangle_repulsions.internal_edge_edge_precomputed.m;i++){
        PRECOMPUTE_PROJECT_EDGE_EDGE<TV>& ee=triangle_repulsions.internal_edge_edge_precomputed(i);
        typename BASE::COLLISION col;
        col.nodes=ee.nodes;
        col.normal=ee.normal;
        Edge_Edge_Weight(ee.weights,col);
        col.relative_velocity=ARRAYS_COMPUTATIONS::Weighted_Sum(col.weights,V.Subset(col.nodes));
        collisions.Append(col);}
}
template class COMBINED_COLLISIONS_EDGE_EDGE<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_EDGE_EDGE<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_EDGE_EDGE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_EDGE_EDGE<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_EDGE_EDGE<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_EDGE_EDGE<VECTOR<double,3> >;
#endif
