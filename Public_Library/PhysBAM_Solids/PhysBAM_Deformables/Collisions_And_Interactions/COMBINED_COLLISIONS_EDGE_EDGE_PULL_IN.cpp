//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_DEFORMABLE_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<TV>::
COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN(TRIANGLE_COLLISIONS<TV>& triangle_collisions_input,TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,const T dt_input,const bool update_swept_hierarchies_input)
    :BASE(triangle_collisions_input,triangle_repulsions_input,update_swept_hierarchies_input),dt(dt_input)
{
    edge_edge_pairs.Remove_All();
    // wipe out edge_edge_interaction_pairs, refill in Discover
    triangle_repulsions.edge_edge_interaction_pairs.Remove_All();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<TV>::
~COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN()
{
}
//#####################################################################
// Function Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<TV>::
Discover(const T dt,const T time)
{
    Discover_Pruned();
}
template<class T> static void Edge_Edge_Weight(const VECTOR<T,2>& weights,typename COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<T,1>,0>::COLLISION& col){}
template<class T> static void Edge_Edge_Weight(const VECTOR<T,2>& weights,typename COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<T,2>,2>::COLLISION& col){col.weights=VECTOR<T,2>(1,-1);}
template<class T> static void Edge_Edge_Weight(const VECTOR<T,2>& weights,typename COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<T,3>,4>::COLLISION& col)
{col.weights=VECTOR<T,4>(1-weights(1),weights(1),-(1-weights(2)),-weights(2));}
//#####################################################################
// Function Update_Repulsion_Pairs_Using_History
//#####################################################################
template<class T> static bool Edge_Edge_Interaction(ARRAY_VIEW<const VECTOR<T,1> > X,ARRAY_VIEW<const VECTOR<T,1> > V,VECTOR<T,2>& weights,
    T& distance,typename COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<T,1>,0>::COLLISION& col,T total_repulsion_thickness,T small_number)
{
    return false;
}
template<class T> static bool Edge_Edge_Interaction(ARRAY_VIEW<const VECTOR<T,2> > X,ARRAY_VIEW<const VECTOR<T,2> > V,VECTOR<T,2>& weights,
    T& distance,typename COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<T,2>,2>::COLLISION& col,T total_repulsion_thickness,T small_number)
{
    col.normal=X(col.nodes(1))-X(col.nodes(2));
    distance=col.normal.Normalize();
    return distance<=total_repulsion_thickness;
}
template<class T> static bool Edge_Edge_Interaction(ARRAY_VIEW<const VECTOR<T,3> > X,ARRAY_VIEW<const VECTOR<T,3> > V,VECTOR<T,2>& weights,
    T& distance,typename COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<T,3>,4>::COLLISION& col,T total_repulsion_thickness,T small_number)
{
    SEGMENT_3D<T> segment1(X.Subset(VECTOR<int,2>(col.nodes(1),col.nodes(2)))),segment2(X.Subset(VECTOR<int,2>(col.nodes(3),col.nodes(4))));
    bool ret=segment1.Edge_Edge_Interaction(segment2,total_repulsion_thickness,distance,col.normal,weights,false);
    segment1.Edge_Edge_Interaction_Data(segment2,V.Subset(col.nodes),distance,col.normal,weights,small_number);
    return ret;
}
VECTOR<int,2> Sort_Edge_Edge_Nodes(const VECTOR<int,2>& nodes)
{
    return nodes.Sorted();
}
VECTOR<int,4> Sort_Edge_Edge_Nodes(const VECTOR<int,4>& nodes)
{
    VECTOR<int,2> pair1(nodes(1),nodes(2)),pair2(nodes(3),nodes(4));
    pair1=pair1.Sorted();pair2=pair2.Sorted();
    if(pair1(1)<pair2(1)) return VECTOR<int,4>(pair1,pair2);
    return VECTOR<int,4>(pair2,pair1);
}
//#####################################################################
// Function Discover_Pruned
//#####################################################################
template<class T,class T_ARRAY1,class T_ARRAY2> T Create_Edges(const T_ARRAY1& X,const T_ARRAY2& V,const T dt,const VECTOR<int,2>& nodes,const ARRAY<T>& repulsion_thickness,POINT_2D<T>& point1,POINT_2D<T>& point2)
{
    point1=POINT_2D<T>(X(nodes[1]))-dt*V(nodes[1]);
    point2=POINT_2D<T>(X(nodes[2]))-dt*V(nodes[2]);
    return EDGE_EDGE_REPULSION_PAIR<typename T_ARRAY1::ELEMENT>::Total_Repulsion_Thickness(repulsion_thickness,nodes);
}
template<class T,class T_ARRAY1,class T_ARRAY2> T Create_Edges(const T_ARRAY1& X,const T_ARRAY2& V,const T dt,const VECTOR<int,4>& nodes,const ARRAY<T>& repulsion_thickness,SEGMENT_3D<T>& segment1,SEGMENT_3D<T>& segment2)
{
    segment1=SEGMENT_3D<T>(X(nodes[1])-dt*V(nodes[1]),X(nodes[2])-dt*V(nodes[2]));
    segment2=SEGMENT_3D<T>(X(nodes[3])-dt*V(nodes[3]),X(nodes[4])-dt*V(nodes[4]));
    return EDGE_EDGE_REPULSION_PAIR<typename T_ARRAY1::ELEMENT>::Total_Repulsion_Thickness(repulsion_thickness,nodes);
}
template<class TV> void COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<TV>::
Discover_Pruned()
{
    for(int i=flagged_for_removal.m;i>0;i--){
        int index=flagged_for_removal(i);
        edge_edge_pairs.Set(collisions(index).nodes,false);
        collisions.Remove_Index_Lazy(index);
        triangle_repulsions.edge_edge_interaction_pairs.Remove_Index_Lazy(index); // lockstep with collisions
    }
    flagged_for_removal.Remove_All();
    /*for(int i=1;i<=collisions.m;i++)
      collisions(i).total_impulse=TV();*/

    ARRAY_VIEW<const TV> X=triangle_collisions.geometry.deformable_body_collection.particles.X;
    ARRAY_VIEW<const TV> V=triangle_collisions.geometry.deformable_body_collection.particles.V;

    ARRAY<TV> V_null;for(int i=1;i<=2*d-2;i++) V_null.Append(TV());

    for(int i=1;i<=triangle_collisions.edge_edge_pairs_internal.m;i++){
        const VECTOR<int,2*d-2>& nodes=Sort_Edge_Edge_Nodes(triangle_collisions.edge_edge_pairs_internal(i));
        if(edge_edge_pairs.Get_Default(nodes,false)) continue;

        typename BASE::COLLISION col;
        col.nodes=nodes;
        VECTOR<T,2> weights;

        //T repulsion_thickness=EDGE_EDGE_REPULSION_PAIR<TV>::Total_Repulsion_Thickness(triangle_collisions.repulsion_thickness,nodes);
        //if(!Edge_Edge_Interaction(X,V,weights,distance,col,repulsion_thickness,triangle_collisions.geometry.small_number)) continue;

        T_EDGE edge1,edge2;T repulsion_thickness=Create_Edges(X,V,dt,nodes,triangle_collisions.repulsion_thickness,edge1,edge2);
        T detection_thickness=1.1*repulsion_thickness;
        T relative_speed,collision_time;
        T distance=0;
        TV normal;
        INDIRECT_ARRAY<ARRAY_VIEW<const TV>,VECTOR<int,2*d-2>&> V_edges(V,nodes);
        if(!edge1.Edge_Edge_Collision(edge2,V_edges,dt,detection_thickness,collision_time,normal,weights,relative_speed,triangle_collisions.geometry.small_number,false)){
            // move edges to final positions
            //Create_Edges(X,V_null,dt,nodes,triangle_collisions.repulsion_thickness,edge1,edge2);
            collision_time=dt;
            if(!Edge_Edge_Interaction(X,V,weights,distance,col,detection_thickness,triangle_collisions.geometry.small_number)) continue;}

        Edge_Edge_Weight(weights,col);
        col.relative_velocity=ARRAYS_COMPUTATIONS::Weighted_Sum(col.weights,V.Subset(col.nodes));
        col.target_relative_velocity=TV::Dot_Product(col.relative_velocity,col.normal)+(repulsion_thickness-distance)/dt;
        col.total_impulse=TV();

        //if(TV::Dot_Product(col.relative_velocity,col.normal)>=0) continue;

        // set up contact
        EDGE_EDGE_REPULSION_PAIR<TV> ee;
        ee.distance=0;
        ee.normal=col.normal;
        ee.weights=weights;
        ee.nodes=nodes;
        triangle_repulsions.edge_edge_interaction_pairs.Append(ee);

        edge_edge_pairs.Set(nodes,true);
        collisions.Append(col);}
}
namespace PhysBAM{
template<> void COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<VECTOR<float,1> >::Discover_Pruned(){PHYSBAM_NOT_IMPLEMENTED();}
template<> void COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<VECTOR<double,1> >::Discover_Pruned(){PHYSBAM_NOT_IMPLEMENTED();}
}
//#####################################################################
template class COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN<VECTOR<double,3> >;
#endif
