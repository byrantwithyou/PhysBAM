//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_DEFORMABLE_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_POINT_FACE_PULL_IN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_POINT_FACE_PULL_IN<TV>::
COMBINED_COLLISIONS_POINT_FACE_PULL_IN(TRIANGLE_COLLISIONS<TV>& triangle_collisions_input,TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,const T dt_input,const bool update_swept_hierarchies_input)
    :BASE(triangle_collisions_input,triangle_repulsions_input,update_swept_hierarchies_input),dt(dt_input)
{
    point_face_pairs.Remove_All();
    // wipe out point_face_interaction_pairs, refill
    triangle_repulsions.point_face_interaction_pairs.Remove_All();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_POINT_FACE_PULL_IN<TV>::
~COMBINED_COLLISIONS_POINT_FACE_PULL_IN()
{
}
//#####################################################################
// Function Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_POINT_FACE_PULL_IN<TV>::
Discover(const T dt,const T time)
{
    Discover_Pruned();
}
template<class T,class TV>
bool Point_Face_Interaction(POINT_SIMPLEX_1D<T>& face,const TV& X,T& thickness,T& distance)
{
    return false;
}
template<class T_FACE,class T,class TV>
bool Point_Face_Interaction(T_FACE& face,const TV& X,T& thickness,T& distance)
{
    return face.Point_Face_Interaction(X,thickness,distance);
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
template<class TV_INT> TV_INT Sort_Point_Face_Nodes(const TV_INT& nodes)
{
    return nodes.Remove_Index(1).Sorted().Insert(nodes(1),1);
}
//#####################################################################
// Function Discover_Pruned
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_POINT_FACE_PULL_IN<TV>::
Discover_Pruned()
{
    for(int i=flagged_for_removal.m;i>0;i--){
        int index=flagged_for_removal(i);
        point_face_pairs.Set(collisions(index).nodes,false);
        collisions.Remove_Index_Lazy(index);
        triangle_repulsions.point_face_interaction_pairs.Remove_Index_Lazy(index); // lockstep with collisions
    }
    flagged_for_removal.Remove_All();

    // set up swept hierarchies
    PARTICLES<TV>& particles=triangle_collisions.geometry.deformable_body_collection.particles;
    if(update_swept_hierarchies){
        ARRAY<bool> recently_modified(particles.array_collection->Size());ARRAYS_COMPUTATIONS::Fill(recently_modified,true);
        //triangle_collisions.Update_Swept_Hierachies_And_Compute_Pairs(particles.X,triangle_collisions.X_save,recently_modified,1.1*ARRAYS_COMPUTATIONS::Max(triangle_collisions.repulsion_thickness));}
        triangle_collisions.Update_Swept_Hierachies_And_Compute_Pairs(particles.X,triangle_collisions.geometry.X_self_collision_free,recently_modified,(T)(1.1*ARRAYS_COMPUTATIONS::Max(triangle_collisions.repulsion_thickness)));}

    ARRAY_VIEW<const TV> X=particles.X;
    ARRAY_VIEW<TV> V=particles.V;
    ARRAY<TV> V_null;V_null.Append(TV());
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,d>&> V_indirect_null(V_null,VECTOR<int,d>::All_Ones_Vector());
    
    for(int i=1;i<=triangle_collisions.point_face_pairs_internal.m;i++){
        const VECTOR<int,d+1>& nodes=Sort_Point_Face_Nodes(triangle_collisions.point_face_pairs_internal(i));

        if(point_face_pairs.Get_Default(nodes,false)) continue; // we've already got one!

        VECTOR<int,d> face_nodes=nodes.Remove_Index(1);
        T repulsion_thickness=POINT_FACE_REPULSION_PAIR<TV>::Total_Repulsion_Thickness(triangle_collisions.repulsion_thickness,nodes);
        TV normal,weights;
        INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,d>&> V_face(V,face_nodes);
        VECTOR<TV,d> X_face;
        for(int j=1;j<=d;j++) X_face(j)=X(face_nodes(j))-dt*V_face(j);
        typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX face(X_face);
        T collision_time;
        T relative_speed;
        TV X1_original=X(nodes(1))-dt*V(nodes(1));
        T detection_thickness=(T)(repulsion_thickness*1.1);
        T distance=0;
        TV face_normal=face.Normal();

        if(!face.Point_Face_Collision(X1_original,V(nodes(1)),V_face,dt,detection_thickness,collision_time,normal,weights,relative_speed,false)){
            for(int j=1;j<=d;j++) face.X(j)=X(face_nodes(j));
            collision_time=dt;
            if(!face.Point_Face_Interaction(X(nodes(1)),TV(),V_indirect_null,detection_thickness,distance,normal,weights,relative_speed,false,false)) continue;}
        // don't want to exit early...but the normal it returns can be screwed up, since it tries to use a zero distance/distance based on the exact collision

        // need distance as well as collision time, since it may not actually collide

        typename BASE::COLLISION col;
        col.normal=TV::Dot_Product(X1_original-face.X(1),face_normal)<0?-face_normal:face_normal;
        col.nodes=nodes;
        //col.normal=normal;
        col.weights=-weights.Insert(-1,1);
        T hit_fraction=collision_time/dt;
        col.relative_velocity=ARRAYS_COMPUTATIONS::Weighted_Sum(col.weights,V.Subset(col.nodes));
        col.target_relative_velocity=hit_fraction*TV::Dot_Product(col.relative_velocity,col.normal)+(repulsion_thickness-distance)/dt;
        col.total_impulse=TV();

        // TODO: not sure about this check
        // allow repulsions
        //if(TV::Dot_Product(col.relative_velocity,col.normal)>=0) continue;
        //if(nodes(1)==2753)
        //    DEBUG_UTILITIES::Debug_Breakpoint();
        POINT_FACE_REPULSION_PAIR<TV> pf;
        pf.distance=0;
        pf.nodes=nodes;
        pf.weights=weights;
        pf.normal=col.normal;
        triangle_repulsions.point_face_interaction_pairs.Append(pf);

        point_face_pairs.Set(nodes,true);
        collisions.Append(col);}
}
namespace PhysBAM{
template<> void COMBINED_COLLISIONS_POINT_FACE_PULL_IN<VECTOR<float,1> >::Discover_Pruned(){PHYSBAM_NOT_IMPLEMENTED();}
template<> void COMBINED_COLLISIONS_POINT_FACE_PULL_IN<VECTOR<double,1> >::Discover_Pruned(){PHYSBAM_NOT_IMPLEMENTED();}
}
//#####################################################################
template class COMBINED_COLLISIONS_POINT_FACE_PULL_IN<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_POINT_FACE_PULL_IN<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_POINT_FACE_PULL_IN<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_POINT_FACE_PULL_IN<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_POINT_FACE_PULL_IN<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_POINT_FACE_PULL_IN<VECTOR<double,3> >;
#endif
