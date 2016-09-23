//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_COLLISIONS  
//##################################################################### 
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/OPERATION_HASH.h>
#include <Core/Data_Structures/SPARSE_UNION_FIND.h>
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Math_Tools/Robust_Arithmetic.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Geometry/Basic_Geometry/POINT_2D.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS_POINT_FACE_VISITOR.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_COLLISIONS<TV>::
TRIANGLE_COLLISIONS(TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry,const ARRAY<T>& repulsion_thickness)
    :geometry(geometry),repulsion_thickness(repulsion_thickness),final_repulsion_youngs_modulus((T)30),final_repulsion_limiter_fraction((T).1),mpi_solids(0),use_gauss_jacobi(false)
{
    // set parameters 
    Set_Collision_Thickness();Set_Restitution_Coefficient();Set_Gauss_Jacobi();
    // set checking
    Compute_Point_Face_Collisions();Compute_Edge_Edge_Collisions();
    Set_Attempts_For_Rigid_Collisions(false);
    // output
    Output_Collision_Results(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TRIANGLE_COLLISIONS<TV>::
~TRIANGLE_COLLISIONS()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters)
{
    Set_Collision_Thickness(triangle_collision_parameters.collisions_collision_thickness);
    if(triangle_collision_parameters.collisions_output_collision_results) Output_Collision_Results();
    Set_Attempts_For_Nonrigid_Collisions(triangle_collision_parameters.collisions_nonrigid_collision_attempts);
    final_repulsion_youngs_modulus=triangle_collision_parameters.collisions_final_repulsion_youngs_modulus;
    final_repulsion_limiter_fraction=triangle_collision_parameters.collisions_final_repulsion_limiter_fraction;
    use_gauss_jacobi=triangle_collision_parameters.use_gauss_jacobi;
}
//#####################################################################
// Function Update_Swept_Hierachies
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Update_Swept_Hierachies_And_Compute_Pairs(ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> X_self_collision_free,ARRAY_VIEW<bool> recently_modified,const T detection_thickness)
{
    // Update swept hierarchies
    LOG::SCOPE scope("updating swept hierarchies","updating swept hierarchies");
    for(int structure_index=0;structure_index<geometry.structure_geometries.m;structure_index++){
        if(!geometry.structure_geometries(structure_index)) continue;
        STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(structure_index);
        if(d==3 && structure.triangulated_surface){
            BOX_HIERARCHY<TV>& hierarchy=structure.Face_Hierarchy();
            structure.triangulated_surface_modified.Resize(hierarchy.box_hierarchy.m,false,false);
            for(int kk=0;kk<structure.triangulated_surface->mesh.elements.m;kk++){
                const VECTOR<int,3>& nodes=structure.triangulated_surface->mesh.elements(kk);
                structure.triangulated_surface_modified(kk)=VECTOR<bool,3>(recently_modified.Subset(nodes)).Contains(true); // TODO: hacking around compiler bug
                if(structure.triangulated_surface_modified(kk))
                    hierarchy.box_hierarchy(kk)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(X.Subset(nodes)),RANGE<TV>::Bounding_Box(X_self_collision_free.Subset(nodes)));}
            hierarchy.Update_Modified_Nonleaf_Boxes(structure.triangulated_surface_modified);}

        if(structure.segmented_curve){
            SEGMENT_HIERARCHY<TV>& hierarchy=*structure.segmented_curve->hierarchy;
            structure.segmented_curve_modified.Resize(hierarchy.box_hierarchy.m,false,false);
            for(int kk=0;kk<structure.segmented_curve->mesh.elements.m;kk++){
                const VECTOR<int,2>& nodes=structure.segmented_curve->mesh.elements(kk);
                structure.segmented_curve_modified(kk)=VECTOR<bool,2>(recently_modified.Subset(nodes)).Contains(true); // TODO: hacking around compiler bug
                if(structure.segmented_curve_modified(kk))
                    hierarchy.box_hierarchy(kk)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(X.Subset(nodes)),RANGE<TV>::Bounding_Box(X_self_collision_free.Subset(nodes)));}
            hierarchy.Update_Modified_Nonleaf_Boxes(structure.segmented_curve_modified);}

        {PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >& hierarchy=structure.particle_hierarchy;
            structure.point_modified.Resize(hierarchy.box_hierarchy.m,false,false);
            for(int kk=0;kk<hierarchy.leaves;kk++){
                int p=structure.active_indices(kk);
                structure.point_modified(kk)=recently_modified(p);
                if(structure.point_modified(kk)) hierarchy.box_hierarchy(kk)=RANGE<TV>::Bounding_Box(X(p),X_self_collision_free(p));}
            hierarchy.Update_Modified_Nonleaf_Boxes(structure.point_modified);}}
    recently_modified.Fill(false);

    // Compute pairs
    point_face_pairs_internal.Remove_All();edge_edge_pairs_internal.Remove_All();
    point_face_pairs_external.Remove_All();edge_edge_pairs_external.Remove_All();
    for(int pair_i=0;pair_i<geometry.interacting_structure_pairs.m;pair_i++){VECTOR<int,2>& pair=geometry.interacting_structure_pairs(pair_i);
        if(compute_point_face_collisions){
            for(int i=0;i<2;i++){if(i==1 && pair[0]==pair[1]) break;
                STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1=*geometry.structure_geometries(pair[i]);
                STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2=*geometry.structure_geometries(pair[1-i]);
                Get_Moving_Faces_Near_Moving_Points(structure_1,structure_2,point_face_pairs_internal,point_face_pairs_external,detection_thickness);}}
        if(compute_edge_edge_collisions){
            STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1=*geometry.structure_geometries(pair[0]);
            STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2=*geometry.structure_geometries(pair[1]);
            Get_Moving_Edges_Near_Moving_Edges(structure_1,structure_2,edge_edge_pairs_internal,edge_edge_pairs_external,detection_thickness);}}
}
template<class TV,class T,int n,int d>
void Test(ARRAY_VIEW<VECTOR<T,d> > X,ARRAY_VIEW<TV> V,VECTOR<int,n> ind,T dt,T collision_thickness,T small_number,ARRAY_VIEW<TV> X2)
{
}
template<class TV,class T,int n>
void Test(ARRAY_VIEW<VECTOR<T,3> > X,ARRAY_VIEW<TV> V,VECTOR<int,n> ind,T dt,T collision_thickness,T small_number,ARRAY_VIEW<TV> X2)
{
    SEGMENT_3D<T> e(X(ind[0]),X(ind[1]));
    SEGMENT_3D<T> e0(X(ind[2]),X(ind[3]));
    SEGMENT_3D<T> e1(X(ind[3]),X(ind[4]));
    SEGMENT_3D<T> e2(X(ind[2]),X(ind[4]));
    VECTOR<int,3> tri_ind(ind[2],ind[3],ind[4]);
    TRIANGLE_3D<T> tri(X.Subset(tri_ind));

    T collision_time;
    TV normal;
    VECTOR<T,TV::m+1> weights;
    bool ee0=e.Edge_Edge_Collision(e0,V.Subset(ind.Remove_Index(4)),dt,collision_thickness,collision_time,normal,weights,small_number,false);
    bool ee1=e.Edge_Edge_Collision(e0,V.Subset(ind.Remove_Index(2)),dt,collision_thickness,collision_time,normal,weights,small_number,false);
    bool ee2=e.Edge_Edge_Collision(e0,V.Subset(ind.Remove_Index(3)),dt,collision_thickness,collision_time,normal,weights,small_number,false);
    bool pf0=tri.Point_Face_Collision(X(ind[0]),V(ind[0]),V.Subset(tri_ind),dt,collision_thickness,collision_time,normal,weights,false);
    bool pf1=tri.Point_Face_Collision(X(ind[1]),V(ind[1]),V.Subset(tri_ind),dt,collision_thickness,collision_time,normal,weights,false);

    bool int0=INTERSECTION::Intersects(e,tri,small_number);
    SEGMENT_3D<T> ex(X(ind[0])+dt*V(ind[0]),X(ind[1])+dt*V(ind[1]));
    ARRAY<TV> BB(X+dt*V);
    TRIANGLE_3D<T> trix(BB.Subset(tri_ind));
    bool int1=INTERSECTION::Intersects(ex,trix,small_number);
    SEGMENT_3D<T> ey(X2(ind[0]),X2(ind[1]));
    TRIANGLE_3D<T> triy(X2.Subset(tri_ind));
    bool int2=INTERSECTION::Intersects(ey,triy,small_number);

    LOG::cout<<"TEST: "<<ind<<" : "<<VECTOR<bool,5>(pf1,pf0,ee1,ee2,ee0)<<"  int "<<int0<<" "<<int1<<" "<<int2<<std::endl;
}
//#####################################################################
// Function Adjust_Velocity_For_Self_Collisions
//#####################################################################
template<class TV> int TRIANGLE_COLLISIONS<TV>::
Adjust_Velocity_For_Self_Collisions(const T dt,const T time,const bool exit_early,const bool change_positions_only)
{
    LOG::SCOPE scope("collisions","checking collisions");
    geometry.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
    geometry.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
    DEFORMABLE_PARTICLES<TV>& full_particles=geometry.deformable_body_collection.particles;
    ARRAY_VIEW<TV> X(full_particles.X),X_self_collision_free(geometry.X_self_collision_free);ARRAY<bool>& modified_full=geometry.modified_full;
    ARRAY<TV> V_save;
    if(change_positions_only){
        V_save=full_particles.V;
        full_particles.V=(full_particles.X-X_self_collision_free)/dt;}
    int collisions=0,collisions_in_attempt=0,
        point_face_collisions=0,edge_edge_collisions=0;
    SPARSE_UNION_FIND<> union_find(full_particles.Size());

    recently_modified_full.Resize(full_particles.Size(),false,false);
    recently_modified_full.Fill(true);
    ARRAY<TV> X_save;
    // input velocities are average V.  Also want original velocities?  Delta may be sufficient.
    geometry.deformable_body_collection.binding_list.Update_Neighbor_Bindings();

    int attempts=0;bool rigid=false;
    while(!attempts || (!exit_early && collisions_in_attempt)){
        attempts++;if(attempts > nonrigid_collision_attempts) rigid=true;
        if(limit_rigid_collision_attempts && attempts > nonrigid_collision_attempts+rigid_collision_attempts) break; // quit and allow intersections
        LOG::SCOPE scope("collision attempt","collision attempt %d",attempts);

        Update_Swept_Hierachies_And_Compute_Pairs(X,X_self_collision_free,recently_modified_full,collision_thickness);

        T attempt_ratio=(T)attempts/(T)nonrigid_collision_attempts;
        if(mpi_solids){
            //int culled=Prune_Non_Intersecting_Pairs(dt,point_face_pairs,edge_edge_pairs,attempt_ratio);
            //culled=mpi_solids->Reduce_Add_Global(culled);LOG::Stat("Pre-gather collision pairs culled",culled);
            LOG::Time("gathering interaction pairs");
            mpi_solids->Gather_Interaction_Pairs(point_face_pairs_external,edge_edge_pairs_external);} // move all collision pairs to root
  
        int exited_early=1;

        // Make a copy of the particles
        impulse_velocities.Resize(full_particles.Size());for(int i=0;i<full_particles.Size();i++) impulse_velocities(i)=full_particles.V(i);
        pf_target_impulses.Resize(point_face_pairs_internal.Size());pf_target_impulses.Fill(TV());
        ee_target_impulses.Resize(edge_edge_pairs_internal.Size());ee_target_impulses.Fill(TV());
        pf_target_weights.Resize(point_face_pairs_internal.Size());pf_target_weights.Fill(VECTOR<T,d+1>());
        ee_target_weights.Resize(edge_edge_pairs_internal.Size());ee_target_weights.Fill(VECTOR<T,d+1>());
        pf_normals.Resize(point_face_pairs_internal.Size());pf_normals.Fill(TV());
        ee_normals.Resize(edge_edge_pairs_internal.Size());ee_normals.Fill(TV());
        pf_old_speeds.Resize(point_face_pairs_internal.Size());pf_old_speeds.Fill(T());
        ee_old_speeds.Resize(edge_edge_pairs_internal.Size());ee_old_speeds.Fill(T());

        // point face first for stability
        point_face_collisions=0;edge_edge_collisions=0;collisions_in_attempt=0;
        if(mpi_solids && mpi_solids->rank==0){
            point_face_collisions+=Adjust_Velocity_For_Point_Face_Collision(dt,rigid,union_find,point_face_pairs_external,attempt_ratio,false,exit_early);
            PHYSBAM_ASSERT(!exit_early);}
        if(mpi_solids){
            LOG::Time("broadcast");
            mpi_solids->Broadcast_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);}
        point_face_collisions+=Adjust_Velocity_For_Point_Face_Collision(dt,rigid,union_find,point_face_pairs_internal,attempt_ratio,false,exit_early);
        if(exit_early && point_face_collisions) goto EXIT_EARLY_AND_COMMUNICATE;
        if(mpi_solids){
            LOG::Time("gather");
            mpi_solids->Gather_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);}
        // edge edge pairs
        if(mpi_solids && mpi_solids->rank==0){
            edge_edge_collisions+=Adjust_Velocity_For_Edge_Edge_Collision(dt,rigid,union_find,edge_edge_pairs_external,attempt_ratio,false,exit_early);
            PHYSBAM_ASSERT(!exit_early);}
        if(mpi_solids){
            LOG::Time("broadcast");
            mpi_solids->Broadcast_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);}
        edge_edge_collisions+=Adjust_Velocity_For_Edge_Edge_Collision(dt,rigid,union_find,edge_edge_pairs_internal,attempt_ratio,false,exit_early);
        if(exit_early && edge_edge_collisions) goto EXIT_EARLY_AND_COMMUNICATE;
        collisions_in_attempt=edge_edge_collisions+point_face_collisions;
        if(mpi_solids) mpi_solids->Gather_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);
        if(mpi_solids){
            LOG::Time("gather");
            mpi_solids->Gather_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);
            LOG::Time("reduce");
            collisions_in_attempt=mpi_solids->Reduce_Add_Global(collisions_in_attempt);}
        collisions+=collisions_in_attempt;

        if(use_gauss_jacobi){
            Scale_And_Apply_Impulses();

            pf_target_impulses.Fill(TV());ee_target_impulses.Fill(TV());
            pf_target_weights.Fill(VECTOR<T,d+1>());ee_target_weights.Fill(VECTOR<T,d+1>());
            pf_normals.Fill(TV());ee_normals.Fill(TV());
            pf_old_speeds.Fill(T());ee_old_speeds.Fill(T());
            Adjust_Velocity_For_Point_Face_Collision(dt,rigid,union_find,point_face_pairs_internal,attempt_ratio,true,exit_early);
            Adjust_Velocity_For_Edge_Edge_Collision(dt,rigid,union_find,edge_edge_pairs_internal,attempt_ratio,true,exit_early);
            
            Scale_And_Apply_Impulses();}

        // Apply rigid motions
        if(rigid && collisions_in_attempt && (!mpi_solids || mpi_solids->rank==0)) Apply_Rigid_Body_Motions(dt,union_find);
        // Update positions
        if(collisions_in_attempt){
            for(int p=0;p<full_particles.Size();p++) if(modified_full(p)) full_particles.X(p)=X_self_collision_free(p)+dt*full_particles.V(p);
            geometry.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions();
            geometry.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();}

        exited_early=0;
        EXIT_EARLY_AND_COMMUNICATE:;
        if(mpi_solids){
            LOG::Time("broadcasting modified data");
            // communicate all data that has modified set to true
            mpi_solids->Broadcast_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);
            LOG::Stop_Time();}

        if(exited_early) break;

        LOG::Stat("processed collisions",collisions_in_attempt);
    }
    if(exit_early && collisions_in_attempt) collisions*=-1; // flag indicating that the collisions were not resolved

//    Test(X_self_collision_free,full_particles.V,VECTOR<int,5>(550,893,1433,1726,1434),dt,collision_thickness,geometry.small_number,full_particles.X);

    if(change_positions_only){
        full_particles.X=X_self_collision_free+dt*full_particles.V;
        full_particles.V=V_save;}
    return collisions;
}
template<> int TRIANGLE_COLLISIONS<VECTOR<float,1> >::Adjust_Velocity_For_Self_Collisions(const T,const T time,const bool,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_COLLISIONS<VECTOR<double,1> >::Adjust_Velocity_For_Self_Collisions(const T,const T time,const bool,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Scale_And_Apply_Impulses
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Scale_And_Apply_Impulses()
{
    // Go through the computed impulses, and compare them to the actual change in velocity seen. Scale back impulses accordingly
    for(int i=0;i<pf_target_impulses.m;i++){
        if(pf_target_impulses(i) == TV()) continue;
        const VECTOR<int,d+1>& nodes=point_face_pairs_internal(i);
        // Compute actual new relative_speed
        TV relative_velocity=-impulse_velocities.Subset(nodes).Weighted_Sum(pf_target_weights(i));
        T relative_speed=relative_velocity.Dot(pf_normals(i));
        if(relative_speed*pf_old_speeds(i)<0){
            T new_scale=-(relative_speed-pf_old_speeds(i))/pf_old_speeds(i);
            pf_target_impulses(i)/=new_scale;}}
    for(int i=0;i<ee_target_impulses.m;i++){
        if(ee_target_impulses(i) == TV()) continue;
        const VECTOR<int,d+1>& nodes=edge_edge_pairs_internal(i);
        // Compute actual new relative_speed
        TV relative_velocity=-impulse_velocities.Subset(nodes).Weighted_Sum(ee_target_weights(i));
        T relative_speed=relative_velocity.Dot(ee_normals(i));
        if(relative_speed*ee_old_speeds(i)<0){
            T new_scale=-(relative_speed-ee_old_speeds(i))/ee_old_speeds(i);
            ee_target_impulses(i)/=new_scale;}}
    // Apply the newly scaled impulses
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    for(int i=0;i<pf_target_impulses.m;i++){
        if(pf_target_impulses(i) == TV()) continue;
        const VECTOR<int,d+1>& nodes=point_face_pairs_internal(i);
        for(int j=0;j<d+1;j++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(j),-pf_target_weights(i)(j)*pf_target_impulses(i));}
    for(int i=0;i<ee_target_impulses.m;i++){
        if(ee_target_impulses(i) == TV()) continue;
        const VECTOR<int,d+1>& nodes=edge_edge_pairs_internal(i);
        for(int j=0;j<d+1;j++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(j),-ee_target_weights(i)(j)*ee_target_impulses(i));}
}
template<> void TRIANGLE_COLLISIONS<VECTOR<float,1> >::Scale_And_Apply_Impulses(){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_COLLISIONS<VECTOR<double,1> >::Scale_And_Apply_Impulses(){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Get_Moving_Faces_Near_Moving_Points
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Get_Moving_Faces_Near_Moving_Points(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,const T detection_thickness)
{
    if(!structure_2.Face_Mesh_Object()) return;

    int old_total=pairs_internal.m;
    TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV> visitor(pairs_internal,pairs_external,structure_1,structure_2,geometry,detection_thickness,mpi_solids);
    if(mpi_solids){
        BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV> > mpi_visitor(visitor,structure_1.point_processor_masks,structure_2.Face_Processor_Masks());
        structure_1.particle_hierarchy.Intersection_List(structure_2.Face_Hierarchy(),mpi_visitor,detection_thickness);}
    else structure_1.particle_hierarchy.Intersection_List(structure_2.Face_Hierarchy(),visitor,detection_thickness);

    if(geometry.output_number_checked && pairs_internal.m-old_total>0) LOG::Stat("checked point face collisions",pairs_internal.m-old_total);
}
template<> void TRIANGLE_COLLISIONS<VECTOR<float,1> >::Get_Moving_Faces_Near_Moving_Points(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_COLLISIONS<VECTOR<double,1> >::Get_Moving_Faces_Near_Moving_Points(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Get_Moving_Edges_Near_Moving_Edges
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Get_Moving_Edges_Near_Moving_Edges(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,const T detection_thickness)
{
//    LOG::SCOPE scope("computing edge edge collision pairs", "computing edge edge collision pairs");
    if(!structure_1.Has_Edges() || !structure_2.Has_Edges()) return;

    int old_total=pairs_internal.m;
    TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV> visitor(pairs_internal,pairs_external,structure_1,structure_2,geometry,detection_thickness,mpi_solids);
    if(mpi_solids){
        BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV> > mpi_visitor(visitor,structure_1.Edge_Processor_Masks(),structure_2.Edge_Processor_Masks());
        structure_1.Edge_Hierarchy().Intersection_List(structure_2.Edge_Hierarchy(),mpi_visitor,detection_thickness);}
    else structure_1.Edge_Hierarchy().Intersection_List(structure_2.Edge_Hierarchy(),visitor,detection_thickness);

    if(geometry.output_number_checked && pairs_internal.m-old_total>0) LOG::Stat("checked edge edge collisions",pairs_internal.m-old_total);
}
template<> void TRIANGLE_COLLISIONS<VECTOR<float,1> >::Get_Moving_Edges_Near_Moving_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T){}
template<> void TRIANGLE_COLLISIONS<VECTOR<double,1> >::Get_Moving_Edges_Near_Moving_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T){}
template<> void TRIANGLE_COLLISIONS<VECTOR<float,2> >::Get_Moving_Edges_Near_Moving_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,2> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,2> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T){}
template<> void TRIANGLE_COLLISIONS<VECTOR<double,2> >::Get_Moving_Edges_Near_Moving_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,2> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,2> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T){}
//#####################################################################
// Function Adjust_Velocity_For_Point_Face_Collision
//#####################################################################
template<class TV> int TRIANGLE_COLLISIONS<TV>::
Adjust_Velocity_For_Point_Face_Collision(const T dt,const bool rigid,SPARSE_UNION_FIND<>& union_find,const ARRAY<VECTOR<int,d+1> >& pairs,
    const T attempt_ratio,const bool final_repulsion_only,const bool exit_early)
{
    final_point_face_repulsions=final_point_face_collisions=0;
    ARRAY<bool>& modified_full=geometry.modified_full;
    ARRAY<int> parent_list;
    int collisions=0,skipping_already_rigid=0;T collision_time;
    for(int i=0;i<pairs.m;i++){const VECTOR<int,d+1>& nodes=pairs(i);
        GAUSS_JACOBI_DATA pf_data(pf_target_impulses(i),pf_target_weights(i),pf_normals(i),pf_old_speeds(i));
        parent_list.Remove_All();
        for(int j=0;j<nodes.m;j++)
            geometry.deformable_body_collection.binding_list.Parents(parent_list,nodes(j));
        if(rigid){
            int root=union_find.Find(parent_list(0));
            bool ok=true;
            for(int j=1;j<parent_list.m;j++)
                if(union_find.Find(parent_list(j))!=root)
                    ok=false;
            if(ok){skipping_already_rigid++;continue;}}
        bool collided;
        if(final_repulsion_only)
            collided=Point_Face_Final_Repulsion(pf_data,nodes,dt,REPULSION_PAIR<TV>::Total_Repulsion_Thickness(repulsion_thickness,nodes),collision_time,attempt_ratio,
                exit_early||rigid);
        else collided=Point_Face_Collision(pf_data,nodes,dt,REPULSION_PAIR<TV>::Total_Repulsion_Thickness(repulsion_thickness,nodes),collision_time,attempt_ratio,exit_early||rigid);
        if(collided){
            collisions++;
            modified_full.Subset(nodes).Fill(true);
            recently_modified_full.Subset(nodes).Fill(true);
            modified_full.Subset(parent_list).Fill(true);
            recently_modified_full.Subset(parent_list).Fill(true);
            for(int k=0;k<nodes.m;k++){
                modified_full.Subset(geometry.deformable_body_collection.binding_list.neighbor_bindings(nodes(k))).Fill(true);
                recently_modified_full.Subset(geometry.deformable_body_collection.binding_list.neighbor_bindings(nodes(k))).Fill(true);}
            if(exit_early){if(output_collision_results) LOG::cout<<"exiting collision checking early - point face collision"<<std::endl;return collisions;}
            if(rigid) union_find.Union(parent_list);}}
    if(output_collision_results && !final_repulsion_only){
        if(final_point_face_repulsions) LOG::Stat("final point face repulsions",final_point_face_repulsions);
        if(final_point_face_collisions) LOG::Stat("final point face collisions",final_point_face_collisions);
        if(skipping_already_rigid) LOG::Stat("rigid collisions where checking was skipped",skipping_already_rigid);
        if(collisions) LOG::cout<<"repeating position update "<<" - point face collisions = "<<collisions<<std::endl;}
    return collisions;
}
//#####################################################################
// Function Prune_Non_Intersecting_Pairs
//#####################################################################
template<class TV> int TRIANGLE_COLLISIONS<TV>::
Prune_Non_Intersecting_Pairs(const T dt,ARRAY<VECTOR<int,d+1> >& point_face_pairs,ARRAY<VECTOR<int,d+1> >& edge_edge_pairs,const T attempt_ratio)
{
    int culled=0;T collision_time;
    for(int i=point_face_pairs.m-1;i>=0;i--){const VECTOR<int,d+1>& nodes=point_face_pairs(i);
        TV temporary_vector;VECTOR<T,d+1> temporary_weights;T temp_old_speed;
        GAUSS_JACOBI_DATA temp_data(temporary_vector,temporary_weights,temporary_vector,temp_old_speed);
        if(!Point_Face_Collision(temp_data,nodes,dt,REPULSION_PAIR<TV>::Total_Repulsion_Thickness(repulsion_thickness,nodes),collision_time,attempt_ratio,true)){
            culled++;
            point_face_pairs.Remove_Index_Lazy(i);}}
    for(int i=edge_edge_pairs.m-1;i>=0;i--){const VECTOR<int,d+1>& nodes=edge_edge_pairs(i);
        TV temporary_vector;VECTOR<T,d+1> temporary_weights;T temp_old_speed;
        GAUSS_JACOBI_DATA temp_data(temporary_vector,temporary_weights,temporary_vector,temp_old_speed);
        if(!Edge_Edge_Collision(temp_data,nodes,dt,collision_time,attempt_ratio,true)){
            culled++;
            edge_edge_pairs.Remove_Index_Lazy(i);}}
    return culled;
}
template<> int TRIANGLE_COLLISIONS<VECTOR<float,1> >::Prune_Non_Intersecting_Pairs(const T,ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T attempt_ratio)
{PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_COLLISIONS<VECTOR<double,1> >::Prune_Non_Intersecting_Pairs(const T,ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T attempt_ratio)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Point_Face_Collision
//#####################################################################
namespace{
template<class T,class TV> inline SEGMENT_2D<T> Create_Final_Face(const SEGMENT_2D<T>& face,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2>&> V_face,const T dt)
{
    return SEGMENT_2D<T>(face.X.x+dt*V_face(0),face.X.y+dt*V_face(1));
}
template<class T,class TV> inline TRIANGLE_3D<T> Create_Final_Face(const TRIANGLE_3D<T>& face,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,3>&> V_face,const T dt)
{
    return TRIANGLE_3D<T>(face.X.x+dt*V_face(0),face.X.y+dt*V_face(1),face.X.z+dt*V_face(2));
}
}
template<class TV> bool TRIANGLE_COLLISIONS<TV>::
Point_Face_Collision(GAUSS_JACOBI_DATA& pf_data,const VECTOR<int,d+1>& nodes,const T dt,const T repulsion_thickness,T& collision_time,const T attempt_ratio,const bool exit_early)
{
    bool return_type=false;

    ARRAY_VIEW<TV> X(geometry.X_self_collision_free); // TODO: this should be a parameter as we usually want X_self_collision_free but not in Stop_All_Nodes
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_effective_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY_VIEW<TV> V_save(impulse_velocities);

    TV normal;VECTOR<T,TV::m+1> weights;T_FACE face(X.Subset(nodes.Remove_Index(0)));
    if(face.Point_Face_Collision(X(nodes[0]),V(nodes[0]),V.Subset(nodes.Remove_Index(0)),dt,collision_thickness,collision_time,normal,weights,exit_early)){
        if(exit_early) return true;
        VECTOR<T,d+1> one_over_m(one_over_mass.Subset(nodes));
        if(use_gauss_jacobi){
            pf_data.old_speed=-V.Subset(nodes).Weighted_Sum(weights).Dot(normal);
            TV impulse=-(1+restitution_coefficient)*pf_data.old_speed/geometry.deformable_body_collection.binding_list.One_Over_Effective_Mass(nodes,weights)*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse);
            pf_data.target_impulse=impulse;
            pf_data.target_weight=weights;
            pf_data.target_normal=normal;}
        else{
            TV impulse=-(1+restitution_coefficient)*pf_data.old_speed/one_over_m.Weighted_Sum(sqr(weights))*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse);}
        return_type=true;}

    if(!use_gauss_jacobi) return_type|=Point_Face_Final_Repulsion(pf_data,nodes,dt,repulsion_thickness,collision_time,attempt_ratio,exit_early);

    return return_type;
}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,1> >::Point_Face_Collision(GAUSS_JACOBI_DATA&,const VECTOR<int,2>&,const T,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,1> >::Point_Face_Collision(GAUSS_JACOBI_DATA&,const VECTOR<int,2>&,const T,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Point_Face_Pull_In
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Point_Face_Pull_In(const VECTOR<int,d+1>& nodes,ARRAY_VIEW<TV> V,const T dt,const T repulsion_thickness)
{
    ARRAY_VIEW<TV> X(geometry.X_self_collision_free);
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_effective_mass;
    VECTOR<int,d> face_nodes=nodes.Remove_Index(0);
    T_FACE face(X.Subset(face_nodes));
    TV normal;
    VECTOR<T,TV::m+1> weights;
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,d>&> V_face(V,face_nodes);
    T collision_time;
    if(face.Point_Face_Collision(X(nodes[0]),V(nodes[0]),V_face,dt,collision_thickness,collision_time,normal,weights,false)){
        T relative_speed=-V_face.Weighted_Sum(weights).Dot(normal);
        VECTOR<T,d+1> one_over_m(one_over_mass.Subset(nodes));
        // want to stop before the collision time, basically...but let's do it by distance
        // assume distance is positive
        /*
          Want to figure out how much velocity to put back.  Do this by getting relative.
          For pull-in, actually want to take a min of how much each collision wants to pull in (in each component?).  I.e., look at the collision normal for each collision and say we
          can't move more than that much in that direction (but other components are okay).  A series of projections/limits?  Say we pretend we went all the way to start with, then each
          projection cuts off part of the velocity.
          If there's sticking friction it can lose other compoments too

        */
        T hit_fraction=collision_time/dt;
        // leave only hit_fraction of normal direction velocity
        T scalar_impulse=clamp(1-hit_fraction-repulsion_thickness/dt,(T)0,(T)1)*relative_speed;
        //T remove_relative_velocity=(distance-repulsion_thickness)/dt;
        //T scalar_impulse=remove_relative_velocity*one_over_m.Average();
        TV impulse=-scalar_impulse/one_over_m.Weighted_Sum(sqr(weights))*normal;
        for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse);}
}
template<> void TRIANGLE_COLLISIONS<VECTOR<float,1> >::Point_Face_Pull_In(const VECTOR<int,2>&,ARRAY_VIEW<VECTOR<float,1> >,const float,const float){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_COLLISIONS<VECTOR<double,1> >::Point_Face_Pull_In(const VECTOR<int,2>&,ARRAY_VIEW<VECTOR<double,1> >,const double,const double){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Point_Face_Final_Repulsion
//#####################################################################
template<class TV> bool TRIANGLE_COLLISIONS<TV>::
Point_Face_Final_Repulsion(GAUSS_JACOBI_DATA& pf_data,const VECTOR<int,d+1>& nodes,const T dt,const T repulsion_thickness,T& collision_time,const T attempt_ratio,const bool exit_early)
{
    bool return_type=false;

    ARRAY_VIEW<TV> X(geometry.X_self_collision_free); // TODO: this should be a parameter as we usually want X_self_collision_free but not in Stop_All_Nodes
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_effective_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY_VIEW<TV> V_save(impulse_velocities);
    VECTOR<int,d> face_nodes=nodes.Remove_Index(0);INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,d>&> V_face(V,face_nodes);

    // check to see if the final position is too close
    TV normal;VECTOR<T,TV::m+1> weights;T_FACE face(X.Subset(nodes.Remove_Index(0)));
    T distance;T_FACE face2=Create_Final_Face(face,V_face,dt);TV point(X(nodes[0])+dt*V(nodes[0]));
    if(face2.Point_Face_Interaction(point,V(nodes[0]),V_face,collision_thickness,distance,normal,weights,true,exit_early)){
        collision_time=dt;
        if(exit_early) return true;
        VECTOR<T,d+1> one_over_m(one_over_mass.Subset(nodes));
        T relative_speed=-V.Subset(nodes).Weighted_Sum(weights).Dot(normal);
        if(!use_gauss_jacobi && relative_speed<0){
            final_point_face_collisions++;
            TV impulse=-(1+restitution_coefficient)*relative_speed/geometry.deformable_body_collection.binding_list.One_Over_Effective_Mass(nodes,weights)*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse);
            face2=Create_Final_Face(face,V_face,dt);point=X(nodes[0])+dt*V(nodes[0]); // update point and face and see if repulsion is still necessary
            if(!face2.Point_Face_Interaction(point,V(nodes[0]),V_face,collision_thickness,distance,normal,weights,false,exit_early)) return true;}
        final_point_face_repulsions++;
        T final_relative_speed=final_repulsion_limiter_fraction*(repulsion_thickness-distance)/dt;
        if(relative_speed >= final_relative_speed) return true;
        T ym_over_mass_times_length=final_repulsion_youngs_modulus*one_over_m.Average()/repulsion_thickness;
        T scalar_impulse=min(final_relative_speed-relative_speed,dt*ym_over_mass_times_length*(repulsion_thickness-distance));
        if(use_gauss_jacobi){
            pf_data.old_speed=relative_speed;
            TV impulse=-scalar_impulse/one_over_m.Weighted_Sum(sqr(weights))*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse,V_save);
            pf_data.target_impulse=impulse;
            pf_data.target_weight=weights;
            pf_data.target_normal=normal;}
        else{
            TV impulse=-scalar_impulse/one_over_m.Weighted_Sum(sqr(weights))*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse);}
        return_type=true;}

    return return_type;
}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,1> >::Point_Face_Final_Repulsion(GAUSS_JACOBI_DATA&,const VECTOR<int,2>&,const T,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,1> >::Point_Face_Final_Repulsion(GAUSS_JACOBI_DATA&,const VECTOR<int,2>&,const T,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Adjust_Velocity_For_Edge_Edge_Collision
//#####################################################################
template<class TV> int TRIANGLE_COLLISIONS<TV>::
Adjust_Velocity_For_Edge_Edge_Collision(const T dt,const bool rigid,SPARSE_UNION_FIND<>& union_find,const ARRAY<VECTOR<int,d+1> >& pairs,
    const T attempt_ratio,const bool final_repulsion_only,const bool exit_early)
{
    final_edge_edge_repulsions=final_edge_edge_collisions=0;
    ARRAY<bool>& modified_full=geometry.modified_full;
    int collisions=0,skipping_already_rigid=0;T collision_time;

    ARRAY<int> parent_list;
    for(int i=0;i<pairs.m;i++){const VECTOR<int,d+1>& nodes=pairs(i);
        GAUSS_JACOBI_DATA ee_data(ee_target_impulses(i),ee_target_weights(i),ee_normals(i),ee_old_speeds(i));
        parent_list.Remove_All();
        for(int j=0;j<nodes.m;j++)
            geometry.deformable_body_collection.binding_list.Parents(parent_list,nodes(j));
        if(rigid){
            int root=union_find.Find(parent_list(0));
            bool ok=true;
            for(int j=1;j<parent_list.m;j++)
                if(union_find.Find(parent_list(j))!=root)
                    ok=false;
            if(ok){skipping_already_rigid++;continue;}}
        bool collided;
        if(final_repulsion_only)
            collided=Edge_Edge_Final_Repulsion(ee_data,nodes,dt,collision_time,attempt_ratio,(exit_early||rigid));
        else collided=Edge_Edge_Collision(ee_data,nodes,dt,collision_time,attempt_ratio,(exit_early||rigid));
        if(collided){collisions++;
            modified_full.Subset(nodes).Fill(true);
            recently_modified_full.Subset(nodes).Fill(true);
            modified_full.Subset(parent_list).Fill(true);
            recently_modified_full.Subset(parent_list).Fill(true);
            for(int k=0;k<nodes.m;k++){
                modified_full.Subset(geometry.deformable_body_collection.binding_list.neighbor_bindings(nodes(k))).Fill(true);
                recently_modified_full.Subset(geometry.deformable_body_collection.binding_list.neighbor_bindings(nodes(k))).Fill(true);}
            if(exit_early){if(output_collision_results) LOG::cout<<"exiting collision checking early - edge collision"<<std::endl;return collisions;}
            if(rigid) union_find.Union(parent_list);}}
    if(output_collision_results && !final_repulsion_only){
        if(final_edge_edge_repulsions) LOG::Stat("final edge edge repulsions",final_edge_edge_repulsions);
        if(final_edge_edge_collisions) LOG::Stat("final edge edge collisions",final_edge_edge_collisions);
        if(skipping_already_rigid) LOG::Stat("rigid collisions where checking was skipped",skipping_already_rigid);
        if(collisions) LOG::cout<<"repeating position update "<<" - edge edge collisions = "<<collisions<<std::endl;}
    return collisions;
}
//#####################################################################
// Function Edge_Edge_Collision
//#####################################################################
template<class TV> bool TRIANGLE_COLLISIONS<TV>::
Edge_Edge_Collision(GAUSS_JACOBI_DATA& ee_data,const VECTOR<int,d+1>& nodes,const T dt,T& collision_time,const T attempt_ratio,const bool exit_early)
{
    ARRAY_VIEW<TV> X(geometry.X_self_collision_free); // TODO: this should be a parameter as we usually want X_self_collision_free but not in Stop_All_Nodes
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_effective_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY_VIEW<TV> V_save(impulse_velocities);
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,d+1>&> V_edges(V,nodes);

    bool return_type=false;
    T_EDGE edge1=SEGMENT_3D<T>(X(nodes[0]),X(nodes[1])),edge2=SEGMENT_3D<T>(X(nodes[2]),X(nodes[3]));

    TV normal;VECTOR<T,TV::m+1> weights;
    if(edge1.Edge_Edge_Collision(edge2,V_edges,dt,collision_thickness,collision_time,normal,weights,geometry.small_number,exit_early)){if(exit_early) return true;
        T relative_speed=-V_edges.Weighted_Sum(weights).Dot(normal);
        VECTOR<T,d+1> one_over_m(one_over_mass.Subset(nodes));
        if(use_gauss_jacobi){
            ee_data.old_speed=relative_speed;
            TV impulse=-(1+restitution_coefficient)*relative_speed/geometry.deformable_body_collection.binding_list.One_Over_Effective_Mass(nodes,weights)*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse,V_save);
            ee_data.target_impulse=impulse;
            ee_data.target_weight=weights;
            ee_data.target_normal=normal;}
        else{
            TV impulse=-(1+restitution_coefficient)*relative_speed/one_over_m.Weighted_Sum(sqr(weights))*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse);}
        return_type=true;}

    if(!use_gauss_jacobi) return_type|=Edge_Edge_Final_Repulsion(ee_data,nodes,dt,collision_time,attempt_ratio,exit_early);
    return return_type;
}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,1> >::Edge_Edge_Collision(GAUSS_JACOBI_DATA&,const VECTOR<int,2>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,1> >::Edge_Edge_Collision(GAUSS_JACOBI_DATA&,const VECTOR<int,2>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,2> >::Edge_Edge_Collision(GAUSS_JACOBI_DATA&,const VECTOR<int,3>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,2> >::Edge_Edge_Collision(GAUSS_JACOBI_DATA&,const VECTOR<int,3>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Edge_Edge_Final_Repulsion
//#####################################################################
template<class TV> bool TRIANGLE_COLLISIONS<TV>::
Edge_Edge_Final_Repulsion(GAUSS_JACOBI_DATA& ee_data,const VECTOR<int,d+1>& nodes,const T dt,T& collision_time,const T attempt_ratio,const bool exit_early)
{
    ARRAY_VIEW<TV> X(geometry.X_self_collision_free); // TODO: this should be a parameter as we usually want X_self_collision_free but not in Stop_All_Nodes
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_effective_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY_VIEW<TV> V_save(impulse_velocities);
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,d+1>&> V_edges(V,nodes);
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,d+1>&> V_save_edges(V_save,nodes);

    bool return_type=false;
    T_EDGE edge1=SEGMENT_3D<T>(X(nodes[0]),X(nodes[1])),edge2=SEGMENT_3D<T>(X(nodes[2]),X(nodes[3]));
    T total_repulsion_thickness=REPULSION_PAIR<TV>::Total_Repulsion_Thickness(repulsion_thickness,nodes);
    TV normal;VECTOR<T,TV::m+1> weights;

    // check to see if the final position is too close - see if the edge x2-x3 intersects the cylinder around x0-x1
    T distance;
    T_EDGE edge1_final=SEGMENT_3D<T>(edge1.X.x+dt*V_edges(0),edge1.X.y+dt*V_edges(1));
    T_EDGE edge2_final=SEGMENT_3D<T>(edge2.X.x+dt*V_edges(2),edge2.X.y+dt*V_edges(3));
    if(edge1_final.Edge_Edge_Interaction(edge2_final,V_edges,collision_thickness,distance,normal,weights,false,geometry.small_number,exit_early)){
        T relative_speed=-V_edges.Weighted_Sum(weights).Dot(normal);
        collision_time=dt;if(exit_early) return true;
        VECTOR<T,d+1> one_over_m_edges(one_over_mass.Subset(nodes));
        if(!use_gauss_jacobi && relative_speed<0){
            final_edge_edge_collisions++;
            TV impulse=-(1+restitution_coefficient)*relative_speed/geometry.deformable_body_collection.binding_list.One_Over_Effective_Mass(nodes,weights)*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse);
            if(!edge1_final.Edge_Edge_Interaction(edge2_final,V_edges,collision_thickness,distance,normal,weights,false,geometry.small_number)) return true;}
        final_edge_edge_repulsions++;
        T final_relative_speed=final_repulsion_limiter_fraction*(total_repulsion_thickness-distance)/dt;
        if(relative_speed >= final_relative_speed) return true;
        T ym_over_mass_times_length=final_repulsion_youngs_modulus*one_over_m_edges.Average()/total_repulsion_thickness;
        T scalar_impulse=min(final_relative_speed-relative_speed,dt*ym_over_mass_times_length*(total_repulsion_thickness-distance));
        if(use_gauss_jacobi){
            ee_data.old_speed=relative_speed;
            TV impulse=-scalar_impulse/one_over_m_edges.Weighted_Sum(sqr(weights))*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse,V_save);
            ee_data.target_impulse=impulse;
            ee_data.target_weight=weights;
            ee_data.target_normal=normal;}
        else{
            TV impulse=-scalar_impulse/one_over_m_edges.Weighted_Sum(sqr(weights))*normal;
            for(int i=0;i<TV::m+1;i++) geometry.deformable_body_collection.binding_list.Apply_Impulse(nodes(i),-weights(i)*impulse);}
        return_type=true;}
    return return_type;
}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,1> >::Edge_Edge_Final_Repulsion(GAUSS_JACOBI_DATA&,const VECTOR<int,2>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,1> >::Edge_Edge_Final_Repulsion(GAUSS_JACOBI_DATA&,const VECTOR<int,2>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,2> >::Edge_Edge_Final_Repulsion(GAUSS_JACOBI_DATA&,const VECTOR<int,3>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,2> >::Edge_Edge_Final_Repulsion(GAUSS_JACOBI_DATA&,const VECTOR<int,3>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Apply_Rigid_Body_Motions
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Apply_Rigid_Body_Motions(const T dt,const SPARSE_UNION_FIND<>& union_find)
{
    HASHTABLE<int,ARRAY<int> > lists;
    for(typename HASHTABLE<int,int>::ITERATOR it(union_find.parents);it.Valid();it.Next())
        lists.Get_Or_Insert(union_find.Find(it.Data())).Append(it.Key());

    for(typename HASHTABLE<int,ARRAY<int> >::ITERATOR it(lists);it.Valid();it.Next())
        it.Data().Append(it.Key());

    DEFORMABLE_PARTICLES<TV>& full_particles=geometry.deformable_body_collection.particles;
    ARRAY<TV>& X_self_collision_free=geometry.X_self_collision_free;
    if(output_collision_results) LOG::cout<<"TOTAL RIGID GROUPS = "<<lists<<std::endl;
    
    T one_over_dt=1/dt;
    for(typename HASHTABLE<int,ARRAY<int> >::ITERATOR it(lists);it.Valid();it.Next()){
        const ARRAY<int>& list=it.Data();
        ARRAY<int> flat_indices;
        geometry.deformable_body_collection.binding_list.Flatten_Indices(flat_indices,list);
        TV average_velocity,center_of_mass;
        T total_mass=0;
        for(int k=0;k<flat_indices.m;k++){int p=flat_indices(k);
            total_mass+=full_particles.mass(p);
            center_of_mass+=full_particles.mass(p)*X_self_collision_free(p);
            average_velocity+=full_particles.mass(p)*(full_particles.X(p)-X_self_collision_free(p));}
        T one_over_total_mass=1/total_mass;
        center_of_mass*=one_over_total_mass;
        average_velocity*=one_over_dt*one_over_total_mass;
        typename TV::SPIN L; // moment of inertia & angular momentum
        SYMMETRIC_MATRIX<T,TV::SPIN::m> inertia;
        for(int k=0;k<flat_indices.m;k++){int p=flat_indices(k);
            TV radial_vector=X_self_collision_free(p)-center_of_mass;
            L+=TV::Cross_Product(radial_vector,full_particles.mass(p)*(full_particles.X(p)-X_self_collision_free(p)));
            inertia+=full_particles.mass(p)*MATRIX<T,TV::SPIN::m,TV::m>::Cross_Product_Matrix(radial_vector).Transposed().Cross_Product_Matrix_Times_With_Symmetric_Result(radial_vector);}
        L*=one_over_dt;
        typename TV::SPIN omega=inertia.Inverse_Times(L);
        ROTATION<TV> R=ROTATION<TV>::From_Rotation_Vector(dt*omega);
        for(int k=0;k<flat_indices.m;k++){int p=flat_indices(k);
            TV new_position=center_of_mass+dt*average_velocity+R.Rotate(X_self_collision_free(p)-center_of_mass);
            full_particles.V(p)=one_over_dt*(new_position-X_self_collision_free(p));
            assert(geometry.modified_full(p));recently_modified_full(p)=true;}
        for(int k=0;k<list.m;k++){
            int p=list(k);
            full_particles.X(p)=geometry.deformable_body_collection.binding_list.Embedded_Position(p);
            full_particles.V(p)=geometry.deformable_body_collection.binding_list.Embedded_Velocity(p);}}
}
//####################################################################
template class TRIANGLE_COLLISIONS<VECTOR<float,1> >;
template class TRIANGLE_COLLISIONS<VECTOR<float,2> >;
template class TRIANGLE_COLLISIONS<VECTOR<float,3> >;
template class TRIANGLE_COLLISIONS<VECTOR<double,1> >;
template class TRIANGLE_COLLISIONS<VECTOR<double,2> >;
template class TRIANGLE_COLLISIONS<VECTOR<double,3> >;
}
