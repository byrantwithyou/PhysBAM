//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_POINT_FACE_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS<TV>::
TRIANGLE_REPULSIONS(TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry)
    :geometry(geometry),youngs_modulus((T)30),spring_limiter_fraction((T).1),perform_attractions(true),attractions_threshold((T)-2),
    hierarchy_repulsion_thickness_multiplier(1),repulsion_thickness_detection_multiplier((T)1.1),mpi_solids(0),use_gauss_jacobi(false)
{
    // set parameters
    Set_Repulsion_Thickness();Set_Friction_Coefficient();
    // set checking
    Compute_Point_Face_Friction();Compute_Edge_Edge_Friction();
    Compute_Point_Face_Inelastic_Collision_Repulsion();Compute_Edge_Edge_Inelastic_Collision_Repulsion();
    Compute_Point_Face_Repulsion();Compute_Edge_Edge_Repulsion();
    // output
    Output_Repulsion_Results(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS<TV>::
~TRIANGLE_REPULSIONS()
{}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Clean_Memory()
{
    repulsion_thickness.Clean_Memory();
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters)
{
    Clean_Memory();
    Set_Repulsion_Thickness(triangle_collision_parameters.collisions_repulsion_thickness);
    perform_attractions=triangle_collision_parameters.perform_repulsion_pair_attractions;
    LOG::cout<<"Solids Parameters set perform_attractions "<<perform_attractions<<std::endl;
    attractions_threshold=triangle_collision_parameters.repulsion_pair_attractions_threshold;
    LOG::cout<<"Solids Parameters set attractions_threshold "<<attractions_threshold<<std::endl;
    if(triangle_collision_parameters.clamp_repulsion_thickness) Clamp_Repulsion_Thickness_With_Meshes(triangle_collision_parameters.collisions_repulsion_clamp_fraction);
    Output_Repulsion_Results(triangle_collision_parameters.collisions_output_repulsion_results);
    Set_Friction_Coefficient(triangle_collision_parameters.self_collision_friction_coefficient);
    if(triangle_collision_parameters.collisions_disable_repulsions_based_on_proximity_factor){
        Turn_Off_Repulsions_Based_On_Current_Proximity(triangle_collision_parameters.collisions_disable_repulsions_based_on_proximity_factor);}
    youngs_modulus=triangle_collision_parameters.repulsions_youngs_modulus;
    spring_limiter_fraction=triangle_collision_parameters.repulsions_limiter_fraction;
    Set_Gauss_Jacobi(triangle_collision_parameters.use_gauss_jacobi);
}
//#####################################################################
// Function Clamp_Repulsion_Thickness_With_Meshes
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Clamp_Repulsion_Thickness_With_Meshes(ARRAY_VIEW<const TV> X,const T scale)
{
    for(int k=0;k<geometry.structure_geometries.m;k++){
        STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(k);
        if(structure.segmented_curve) for(int s=0;s<structure.segmented_curve->mesh.elements.m;s++){
            int i,j;structure.segmented_curve->mesh.elements(s).Get(i,j);T d=scale*(X(i)-X(j)).Magnitude();
            repulsion_thickness(i)=min(repulsion_thickness(i),d);repulsion_thickness(j)=min(repulsion_thickness(j),d);}}
}
template<> void TRIANGLE_REPULSIONS<VECTOR<float,1> >::Clamp_Repulsion_Thickness_With_Meshes(ARRAY_VIEW<const VECTOR<T,1> >,const T){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_REPULSIONS<VECTOR<double,1> >::Clamp_Repulsion_Thickness_With_Meshes(ARRAY_VIEW<const VECTOR<T,1> >,const T){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Turn_Off_Repulsions_Based_On_Current_Proximity
//#####################################################################
static inline VECTOR<int,3> Flip_Edges(const VECTOR<int,3>& nodes)
{
    PHYSBAM_FATAL_ERROR();
}
static inline VECTOR<int,4> Flip_Edges(const VECTOR<int,4>& nodes)
{
    return VECTOR<int,4>(nodes[2],nodes[3],nodes[0],nodes[1]);
}
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Turn_Off_Repulsions_Based_On_Current_Proximity(const T extra_factor_on_distance)
{
    if(geometry.intersecting_point_face_pairs.Size() || geometry.intersecting_edge_edge_pairs.Size()) PHYSBAM_FATAL_ERROR();
    const ARRAY<TV>& X_self_collision_free=geometry.X_self_collision_free;
    bool output_number_checked_save=geometry.output_number_checked;geometry.output_number_checked=false;
    repulsion_thickness_detection_multiplier*=extra_factor_on_distance; // enlarge to catch repulsions at a further distance away
    Update_Faces_And_Hierarchies_With_Collision_Free_Positions(0);
    // find interacting pairs
    point_face_interaction_pairs.Remove_All();edge_edge_interaction_pairs.Remove_All();
    // TODO: consider checking cross-structure interactions as well
    for(int a=0;a<geometry.interacting_structure_pairs.m;a++){VECTOR<int,2>& pair=geometry.interacting_structure_pairs(a);
        if(pair[0]!=pair[1]) continue;
        STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(pair[0]);
        Get_Faces_Near_Points(structure,structure,X_self_collision_free,false);
        Get_Edges_Near_Edges(structure,structure,X_self_collision_free,false);}
    // construct omit hashtables
    omit_point_face_repulsion_pairs.Remove_All();
    for(int k=0;k<point_face_interaction_pairs.m;k++) omit_point_face_repulsion_pairs.Insert(point_face_interaction_pairs(k).nodes);
    omit_edge_edge_repulsion_pairs.Remove_All();
    for(int k=0;k<edge_edge_interaction_pairs.m;k++){const VECTOR<int,d+1>& nodes=edge_edge_interaction_pairs(k).nodes;
        omit_edge_edge_repulsion_pairs.Insert(nodes); // add edge edge pairs in both orders since the hierarchy can change after this is called
        omit_edge_edge_repulsion_pairs.Insert(Flip_Edges(nodes));}
    // cleanup
    repulsion_thickness_detection_multiplier/=extra_factor_on_distance;
    geometry.output_number_checked=output_number_checked_save;
    LOG::Stat("omitted point face repulsion pairs",point_face_interaction_pairs.m);
    LOG::Stat("omitted edge edge repulsion pairs",edge_edge_interaction_pairs.m);
    point_face_interaction_pairs.Remove_All();edge_edge_interaction_pairs.Remove_All();
}
template<> void TRIANGLE_REPULSIONS<VECTOR<float,1> >::Turn_Off_Repulsions_Based_On_Current_Proximity(const T){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_REPULSIONS<VECTOR<double,1> >::Turn_Off_Repulsions_Based_On_Current_Proximity(const T){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Update_Faces_And_Hierarchies_With_Collision_Free_Positions
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Update_Faces_And_Hierarchies_With_Collision_Free_Positions(const ARRAY_VIEW<TV>* X_other)
{
    LOG::SCOPE scope("Update Triangles and Hierarchies","Update Triangles and Hierarchies");
    const ARRAY_VIEW<TV>& X=X_other?*X_other:geometry.X_self_collision_free;
    T multiplier=repulsion_thickness_detection_multiplier*hierarchy_repulsion_thickness_multiplier;
    for(int k=0;k<geometry.structure_geometries.m;k++)
        geometry.structure_geometries(k)->Update_Faces_And_Hierarchies_With_Collision_Free_Positions(repulsion_thickness,multiplier,X);
}
//#####################################################################
// Function Compute_Interaction_Pairs
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Compute_Interaction_Pairs(ARRAY_VIEW<const TV> X_other)
{
    LOG::SCOPE scope("computing repulsion pairs", "computing repulsion pairs");
    point_face_interaction_pairs.Remove_All();edge_edge_interaction_pairs.Remove_All();
    for(int pair_i=0;pair_i<geometry.interacting_structure_pairs.m;pair_i++){VECTOR<int,2>& pair=geometry.interacting_structure_pairs(pair_i);
        for(int i=0;i<2;i++){if(i==1 && pair[0]==pair[1]) break;
            if(compute_point_face_friction || compute_point_face_inelastic_collision_repulsion || compute_point_face_repulsion){
                Get_Faces_Near_Points(*geometry.structure_geometries(pair[i]),*geometry.structure_geometries(pair[1-i]),X_other,true);}}
        if(compute_edge_edge_friction || compute_edge_edge_inelastic_collision_repulsion || compute_edge_edge_repulsion){
            Get_Edges_Near_Edges(*geometry.structure_geometries(pair[0]),*geometry.structure_geometries(pair[1]),X_other,true);}}

    if(mpi_solids){
        mpi_solids->Gather_Interaction_Pairs(point_face_interaction_pairs,edge_edge_interaction_pairs);

        mpi_solids->Distribute_Repulsion_Pairs(point_face_interaction_pairs,point_face_send_particles,point_face_receive_particles,point_face_boundary_pairs,point_face_internal_pairs);
        
        mpi_solids->Distribute_Repulsion_Pairs(edge_edge_interaction_pairs,edge_edge_send_particles,edge_edge_receive_particles,edge_edge_boundary_pairs,edge_edge_internal_pairs);}
}
//#####################################################################
// Function Adjust_Velocity_For_Self_Repulsion
//#####################################################################
template<class T,class TV,class T_ARRAY> void
Edge_Edge_Interaction_Data_Helper(ARRAY_VIEW<const VECTOR<T,2> > X,REPULSION_PAIR<TV>& pair,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,3>&> V_edges,const T& small_number)
{
    PHYSBAM_FATAL_ERROR();
}
template<class T,class TV,class T_ARRAY> void
Edge_Edge_Interaction_Data_Helper(ARRAY_VIEW<const VECTOR<T,3> > X,REPULSION_PAIR<TV>& pair,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V_edges,const T& small_number)
{
    SEGMENT_3D<T> segment1(X.Subset(VECTOR<int,2>(pair.nodes[0],pair.nodes[1]))),segment2(X.Subset(VECTOR<int,2>(pair.nodes[2],pair.nodes[3])));
    segment1.Edge_Edge_Interaction_Data(segment2,V_edges,pair.distance,pair.normal,pair.weights,small_number);
}
template<class TV> int TRIANGLE_REPULSIONS<TV>::
Adjust_Velocity_For_Self_Repulsion(const T dt,bool use_saved_pairs)
{
    PHYSBAM_ASSERT(!mpi_solids); // Only per time step repulsions are supported in MPI

    LOG::SCOPE scope("repulsions","checking repulsions");
    ARRAY_VIEW<const TV> X_self_collision_free(geometry.X_self_collision_free);
    ARRAY<bool>& modified_full=geometry.modified_full;
    ARRAY_VIEW<const TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY<REPULSION_PAIR<TV> > point_face_pairs(point_face_interaction_pairs);
    ARRAY<REPULSION_PAIR<TV> > edge_edge_pairs(edge_edge_interaction_pairs);

    for(int pair_index=0;pair_index<point_face_pairs.m;pair_index++){
        REPULSION_PAIR<TV>& pair=point_face_pairs(pair_index);
        T_FACE face(X_self_collision_free.Subset(pair.nodes.Remove_Index(0)));
        face.Point_Face_Interaction_Data(X_self_collision_free(pair.nodes[0]),pair.distance,pair.normal,pair.weights,perform_attractions);
        modified_full.Subset(pair.nodes).Fill(true);}

    // TODO: do we need update binding here?
    // TODO: MPI check fragments in super fragment for whether this fragment is in the processor?
    geometry.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();

    for(int pair_index=0;pair_index<edge_edge_pairs.m;pair_index++){
        REPULSION_PAIR<TV>& pair=edge_edge_pairs(pair_index);
        // Note: don't call this, all it does is mess up the normal and has already been called when the pairs are created
        //Edge_Edge_Interaction_Data_Helper(X_self_collision_free,pair,V.Subset(pair.nodes),geometry.small_number);
        INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,TV::m+1>&> modified_subset=modified_full.Subset(pair.nodes);
        modified_subset.Fill(true);}

    int repulsions=Apply_Repulsions_To_Velocities(dt,point_face_pairs,edge_edge_pairs,true,use_saved_pairs);
    LOG::Stat("adjusting velocity for repulsions",repulsions);
    if(repulsions) for(int p=0;p<geometry.deformable_body_collection.particles.Size();p++) if(modified_full(p)) geometry.deformable_body_collection.particles.X(p)=X_self_collision_free(p)+dt*V(p);
// TODO: this is broken
//    for(int i=0;i<geometry.deformable_body_collection.rigid_body_particles.Size();i++){
//        geometry.deformable_body_collection.rigid_body_particles.Frame(i)=geometry.rigid_body_particle_state_collision_free(i).frame;
//        geometry.deformable_body_collection.rigid_body_particles.Euler_Step_Position(VECTOR<int,1>(i),dt);}
    return repulsions;
}
template<> int TRIANGLE_REPULSIONS<VECTOR<float,1> >::Adjust_Velocity_For_Self_Repulsion(const T,bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<double,1> >::Adjust_Velocity_For_Self_Repulsion(const T,bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Pair_Is_Separating
//#####################################################################
template<class TV> bool Pair_Is_Separating(REPULSION_PAIR<TV>& pair,ARRAY_VIEW<const TV> V)
{
    TV relative_velocity=-V.Subset(pair.nodes).Weighted_Sum(pair.weights);
    return TV::Dot_Product(relative_velocity,pair.normal)>=0;
}
//#####################################################################
// Function Update_Repulsion_Pairs_Using_History
//#####################################################################
template<class T,class TV> bool Edge_Edge_Interaction_Helper(ARRAY_VIEW<const VECTOR<T,2> > X,REPULSION_PAIR<TV>& pair,const ARRAY<T>& repulsion_thickness,
    const T repulsion_thickness_detection_multiplier)
{
    pair.normal=X(pair.nodes[0])-X(pair.nodes[1]);
    pair.distance=pair.normal.Magnitude();
    T total_repulsion_thickness=repulsion_thickness_detection_multiplier*pair.Total_Repulsion_Thickness(repulsion_thickness);
    return pair.distance<=total_repulsion_thickness;
}
template<class T,class TV> bool Edge_Edge_Interaction_Helper(ARRAY_VIEW<const VECTOR<T,3> > X,REPULSION_PAIR<TV>& pair,const ARRAY<T>& repulsion_thickness,
    const T repulsion_thickness_detection_multiplier)
{
    SEGMENT_3D<T> segment(X.Subset(VECTOR<int,2>(pair.nodes[0],pair.nodes[1])));
    T total_repulsion_thickness=repulsion_thickness_detection_multiplier*pair.Total_Repulsion_Thickness(repulsion_thickness);
    return segment.Edge_Edge_Interaction(SEGMENT_3D<T>(X(pair.nodes[2]),X(pair.nodes[3])),total_repulsion_thickness,pair.distance,pair.normal,pair.weights,false);
}
template<class TV> template<class T_ARRAY1,class T_ARRAY2> void TRIANGLE_REPULSIONS<TV>::
Update_Repulsion_Pairs_Using_History(T_ARRAY1& point_face_pairs,T_ARRAY2& edge_edge_pairs,bool prune_separating)
{
    ARRAY_VIEW<const TV> X(geometry.deformable_body_collection.particles.X),V(geometry.deformable_body_collection.particles.V);
    for(int pair_index=point_face_pairs.Size()-1;pair_index>=0;pair_index--){
        REPULSION_PAIR<TV>& pair=point_face_pairs(pair_index);VECTOR<int,d> face_nodes=pair.nodes.Remove_Index(0);
        T_FACE face(X.Subset(face_nodes));
        if(!face.Point_Face_Interaction(X(pair.nodes(0)),repulsion_thickness_detection_multiplier*pair.Total_Repulsion_Thickness(repulsion_thickness),false,pair.distance) ||
            (prune_separating && Pair_Is_Separating(pair,V))){
            point_face_pairs.Remove_Index_Lazy(pair_index);}
        else face.Point_Face_Interaction_Data(X(pair.nodes[0]),pair.distance,pair.normal,pair.weights,perform_attractions);}

    // TODO: do we need update binding here? (what about the other fragments?)
    geometry.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();

    for(int pair_index=edge_edge_pairs.Size()-1;pair_index>=0;pair_index--){
        REPULSION_PAIR<TV>& pair=edge_edge_pairs(pair_index);
        if(!Edge_Edge_Interaction_Helper(X,pair,repulsion_thickness,repulsion_thickness_detection_multiplier) ||
            (prune_separating && Pair_Is_Separating(pair,V))){
            edge_edge_pairs.Remove_Index_Lazy(pair_index);}
        else{
            Edge_Edge_Interaction_Data_Helper(X,pair,V.Subset(pair.nodes),geometry.small_number);
            if(perform_attractions && TV::Dot_Product(pair.normal,pair.collision_free_normal)<attractions_threshold){pair.distance*=-1;pair.normal*=-1;}}}
}
//#####################################################################
// Function Adjust_Velocity_For_Self_Repulsion_Using_History
//#####################################################################
template<class TV> int TRIANGLE_REPULSIONS<TV>::
Adjust_Velocity_For_Self_Repulsion_Using_History(const T dt,const bool use_repulsions,bool use_saved_pairs)
{
   // TODO: MPI
    LOG::SCOPE scope("history repulsions","checking history repulsions");
    if(use_repulsions) geometry.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true); // otherwise assume we're in the velocity update, and positions are correct
    geometry.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Velocities(true);
    int repulsions=0;
    if(mpi_solids){
        ARRAY_VIEW<TV> X(geometry.deformable_body_collection.particles.X),V(geometry.deformable_body_collection.particles.V),X_self_collision_free(geometry.X_self_collision_free);
        mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);
        mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);
        Update_Repulsion_Pairs_Using_History(point_face_boundary_pairs,edge_edge_boundary_pairs,false);
        Update_Repulsion_Pairs_Using_History(point_face_internal_pairs,edge_edge_internal_pairs,false);
        repulsions+=Apply_Repulsions_To_Velocities(dt,point_face_boundary_pairs,edge_edge_boundary_pairs,point_face_internal_pairs,edge_edge_internal_pairs,use_repulsions);
    }
    else{
        Update_Repulsion_Pairs_Using_History(point_face_interaction_pairs,edge_edge_interaction_pairs,false);
        repulsions+=Apply_Repulsions_To_Velocities(dt,point_face_interaction_pairs,edge_edge_interaction_pairs,use_repulsions,use_saved_pairs);}

    LOG::Stat("adjusting velocity for repulsions using history",repulsions);
    return repulsions;
}
template<> int TRIANGLE_REPULSIONS<VECTOR<float,1> >::Adjust_Velocity_For_Self_Repulsion_Using_History(const T,const bool,bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<double,1> >::Adjust_Velocity_For_Self_Repulsion_Using_History(const T,const bool,bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Adjust_Velocity_For_Self_Repulsion
//#####################################################################
template<class TV> template<class T_ARRAY1,class T_ARRAY2> int TRIANGLE_REPULSIONS<TV>::
Apply_Repulsions_To_Velocities(const T dt,T_ARRAY1& point_face_interaction_pairs,T_ARRAY2& edge_edge_interaction_pairs,const bool use_repulsions,bool use_saved_pairs)
{
    int repulsions=0;

    if(compute_point_face_friction){
        Adjust_Velocity_For_Repulsion(dt,point_face_interaction_pairs,use_repulsions,true,use_repulsions,pf_target_impulses,
            pf_old_speeds,pf_normals,point_face_inelastic_collision_repulsion_attempts); // if not using repulsions, only do friction for inelastic component
        repulsions+=point_face_interaction_pairs.Size();
        if(output_repulsion_results) LOG::Stat("total point face friction",point_face_interaction_pairs.Size());}
    if(compute_edge_edge_friction){
        Adjust_Velocity_For_Repulsion(dt,edge_edge_interaction_pairs,use_repulsions,true,use_repulsions,ee_target_impulses,
            ee_old_speeds,ee_normals,edge_edge_inelastic_collision_repulsion_attempts); // if not using repulsions, only do friction for inelastic component
        repulsions+=edge_edge_interaction_pairs.Size();
        if(output_repulsion_results) LOG::Stat("total edge edge friction",edge_edge_interaction_pairs.Size());}
    if(compute_point_face_inelastic_collision_repulsion){
        Adjust_Velocity_For_Repulsion(dt,point_face_interaction_pairs,false,false,use_repulsions,pf_target_impulses,
            pf_old_speeds,pf_normals,point_face_inelastic_collision_repulsion_attempts);
        repulsions+=point_face_inelastic_collision_repulsion_attempts*point_face_interaction_pairs.Size();
        if(output_repulsion_results) LOG::Stat("total point face inelastic collision repulsion",point_face_inelastic_collision_repulsion_attempts*point_face_interaction_pairs.Size());}
    if(compute_edge_edge_inelastic_collision_repulsion){
        Adjust_Velocity_For_Repulsion(dt,edge_edge_interaction_pairs,false,false,use_repulsions,ee_target_impulses,
            ee_old_speeds,ee_normals,edge_edge_inelastic_collision_repulsion_attempts);
        repulsions+=edge_edge_inelastic_collision_repulsion_attempts*edge_edge_interaction_pairs.Size();
        if(output_repulsion_results) LOG::Stat("total edge edge inelastic collision repulsion",edge_edge_inelastic_collision_repulsion_attempts*edge_edge_interaction_pairs.Size());}
    if(use_repulsions){
        if(compute_point_face_repulsion){
            Adjust_Velocity_For_Repulsion(dt,point_face_interaction_pairs,true,false,use_repulsions,pf_target_impulses,
                pf_old_speeds,pf_normals,point_face_inelastic_collision_repulsion_attempts);
            repulsions+=point_face_interaction_pairs.Size();
            if(output_repulsion_results) LOG::Stat("total point face repulsion",point_face_interaction_pairs.Size());}
        if(compute_edge_edge_repulsion){
            Adjust_Velocity_For_Repulsion(dt,edge_edge_interaction_pairs,true,false,use_repulsions,ee_target_impulses,
                ee_old_speeds,ee_normals,edge_edge_inelastic_collision_repulsion_attempts);
            repulsions+=edge_edge_interaction_pairs.Size();
            if(output_repulsion_results) LOG::Stat("total edge edge repulsion",edge_edge_interaction_pairs.Size());}}
    return repulsions;
}
template<class TV> template<class T_ARRAY1,class T_ARRAY2> int TRIANGLE_REPULSIONS<TV>::
Apply_Repulsions_To_Velocities(const T dt,T_ARRAY1& point_face_boundary_pairs,T_ARRAY2& edge_edge_boundary_pairs,
    T_ARRAY1& point_face_internal_pairs,T_ARRAY2& edge_edge_internal_pairs,const bool use_repulsions)
{
    ARRAY_VIEW<TV> X(geometry.deformable_body_collection.particles.X),V(geometry.deformable_body_collection.particles.V),X_self_collision_free(geometry.X_self_collision_free);
    int repulsions=0;bool used=false;
    if(compute_point_face_friction){
        used=true;
        Adjust_Velocity_For_Repulsion(dt,point_face_boundary_pairs,use_repulsions,true,use_repulsions,pf_target_impulses,
            pf_old_speeds,pf_normals,point_face_inelastic_collision_repulsion_attempts); // if not using repulsions, only do friction for inelastic component
        mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);
        Adjust_Velocity_For_Repulsion(dt,point_face_internal_pairs,use_repulsions,true,use_repulsions,pf_target_impulses,
            pf_old_speeds,pf_normals,point_face_inelastic_collision_repulsion_attempts); // if not using repulsions, only do friction for inelastic component
        int applied_repulsions=point_face_boundary_pairs.Size()+point_face_internal_pairs.Size();
        repulsions+=applied_repulsions;
        if(output_repulsion_results) LOG::Stat("total point face friction",applied_repulsions);}
    if(compute_edge_edge_friction){
        if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);used=true;
        Adjust_Velocity_For_Repulsion(dt,edge_edge_boundary_pairs,use_repulsions,true,use_repulsions,ee_target_impulses,
            ee_old_speeds,ee_normals,edge_edge_inelastic_collision_repulsion_attempts); // if not using repulsions, only do friction for inelastic component
        mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);
        Adjust_Velocity_For_Repulsion(dt,edge_edge_internal_pairs,use_repulsions,true,use_repulsions,ee_target_impulses,
            ee_old_speeds,ee_normals,edge_edge_inelastic_collision_repulsion_attempts); // if not using repulsions, only do friction for inelastic component
        int applied_repulsions=edge_edge_boundary_pairs.Size()+edge_edge_internal_pairs.Size();
        repulsions+=applied_repulsions;
        if(output_repulsion_results) LOG::Stat("total edge edge friction",applied_repulsions);}
    if(compute_point_face_inelastic_collision_repulsion){
        if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);used=true;
        Adjust_Velocity_For_Repulsion(dt,point_face_boundary_pairs,false,false,use_repulsions,pf_target_impulses,
            pf_old_speeds,pf_normals,point_face_inelastic_collision_repulsion_attempts);
        mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);
        Adjust_Velocity_For_Repulsion(dt,point_face_internal_pairs,false,false,use_repulsions,pf_target_impulses,
            pf_old_speeds,pf_normals,point_face_inelastic_collision_repulsion_attempts);
        int applied_repulsions=point_face_inelastic_collision_repulsion_attempts*(point_face_boundary_pairs.Size()+point_face_internal_pairs.Size());
        repulsions+=applied_repulsions;
        if(output_repulsion_results) LOG::Stat("total point face inelastic collision repulsion",applied_repulsions);}
    if(compute_edge_edge_inelastic_collision_repulsion){
        if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);used=true;
        Adjust_Velocity_For_Repulsion(dt,edge_edge_boundary_pairs,false,false,use_repulsions,ee_target_impulses,
            ee_old_speeds,ee_normals,edge_edge_inelastic_collision_repulsion_attempts);
        mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);
        Adjust_Velocity_For_Repulsion(dt,edge_edge_internal_pairs,false,false,use_repulsions,ee_target_impulses,
            ee_old_speeds,ee_normals,edge_edge_inelastic_collision_repulsion_attempts);
        int applied_repulsions=edge_edge_inelastic_collision_repulsion_attempts*(edge_edge_boundary_pairs.Size()+edge_edge_internal_pairs.Size());
        repulsions+=applied_repulsions;
        if(output_repulsion_results) LOG::Stat("total edge edge inelastic collision repulsion",applied_repulsions);}
    if(use_repulsions){
        if(compute_point_face_repulsion){
            if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);used=true;
            Adjust_Velocity_For_Repulsion(dt,point_face_boundary_pairs,true,false,use_repulsions,pf_target_impulses,
                pf_old_speeds,pf_normals,point_face_inelastic_collision_repulsion_attempts);
            mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);
            Adjust_Velocity_For_Repulsion(dt,point_face_internal_pairs,true,false,use_repulsions,pf_target_impulses,
                pf_old_speeds,pf_normals,point_face_inelastic_collision_repulsion_attempts);
            int applied_repulsions=point_face_boundary_pairs.Size()+point_face_internal_pairs.Size();
            repulsions+=applied_repulsions;
            if(output_repulsion_results) LOG::Stat("total point face repulsion",applied_repulsions);}
        if(compute_edge_edge_repulsion){
            if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);used=true;
            Adjust_Velocity_For_Repulsion(dt,edge_edge_boundary_pairs,true,false,use_repulsions,ee_target_impulses,
                ee_old_speeds,ee_normals,edge_edge_inelastic_collision_repulsion_attempts);
            mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);
            Adjust_Velocity_For_Repulsion(dt,edge_edge_internal_pairs,true,false,use_repulsions,ee_target_impulses,
                ee_old_speeds,ee_normals,edge_edge_inelastic_collision_repulsion_attempts);
            int applied_repulsions=edge_edge_boundary_pairs.Size()+edge_edge_internal_pairs.Size();
            repulsions+=applied_repulsions;
            if(output_repulsion_results) LOG::Stat("total edge edge repulsion",applied_repulsions);}}
    return repulsions;
}
//#####################################################################
// Function Get_Faces_Near_Points
//#####################################################################
template<class TV> int TRIANGLE_REPULSIONS<TV>::
Get_Faces_Near_Points(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY_VIEW<const TV> X_other,const bool use_processor_cull)
{
    if(!structure_2.Face_Mesh_Object()) return 0;

    int start_count=point_face_interaction_pairs.m;int pruned=0;
    TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<TV> visitor(point_face_interaction_pairs,structure_1,structure_2,X_other,*this,pruned);
    if(use_processor_cull && mpi_solids){
        BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<TV> > mpi_visitor(visitor,structure_1.point_processor_masks,structure_2.Face_Processor_Masks());
        structure_1.particle_hierarchy.Intersection_List(structure_2.Face_Hierarchy(),mpi_visitor,ZERO());}
    else structure_1.particle_hierarchy.Intersection_List(structure_2.Face_Hierarchy(),visitor,ZERO());

    int checked=point_face_interaction_pairs.m-start_count;
    if(checked){
//        LOG::cout<<"pruned "<<pruned<<" out of "<<checked+pruned<<std::endl;
        if(geometry.output_number_checked) LOG::Stat("checked point face repulsions",checked);}
    return checked;
}
template<> int TRIANGLE_REPULSIONS<VECTOR<float,1> >::Get_Faces_Near_Points(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY_VIEW<const VECTOR<T,1> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<double,1> >::Get_Faces_Near_Points(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY_VIEW<const VECTOR<T,1> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Output_Interaction_Pairs
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Output_Interaction_Pairs(const STREAM_TYPE stream_type,const std::string& filename) const
{
    FILE_UTILITIES::Write_To_File(stream_type,filename,point_face_interaction_pairs,edge_edge_interaction_pairs);
}
template<> void TRIANGLE_REPULSIONS<VECTOR<float,1> >::Output_Interaction_Pairs(const STREAM_TYPE,const std::string&) const {PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_REPULSIONS<VECTOR<double,1> >::Output_Interaction_Pairs(const STREAM_TYPE,const std::string&) const {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
// Function Get_Edges_Near_Edges
//#####################################################################
template<class TV> int TRIANGLE_REPULSIONS<TV>::
Get_Edges_Near_Edges(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY_VIEW<const TV> X_other,const bool use_processor_cull)
{
    if(!structure_1.Has_Edges() || !structure_2.Has_Edges()) return 0;

    int start_count=edge_edge_interaction_pairs.m;int pruned=0;
    TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<TV> visitor(edge_edge_interaction_pairs,structure_1,structure_2,X_other,*this,pruned);
    if(use_processor_cull && mpi_solids){
        BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<TV> > mpi_visitor(visitor,structure_1.Edge_Processor_Masks(),structure_2.Edge_Processor_Masks());
        structure_1.Edge_Hierarchy().Intersection_List(structure_2.Edge_Hierarchy(),mpi_visitor,ZERO());}
    else structure_1.Edge_Hierarchy().Intersection_List(structure_2.Edge_Hierarchy(),visitor,ZERO());

    int checked=edge_edge_interaction_pairs.m-start_count;
    if(checked){
        //LOG::cout<<"pruned "<<pruned<<" out of "<<checked+pruned<<std::endl;
        if(geometry.output_number_checked) LOG::Stat("checked edge edge repulsions",checked);}
    return checked;
}
template<> int TRIANGLE_REPULSIONS<VECTOR<float,1> >::Get_Edges_Near_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY_VIEW<const VECTOR<T,1> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<double,1> >::Get_Edges_Near_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY_VIEW<const VECTOR<T,1> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<float,2> >::Get_Edges_Near_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,2> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,2> >&,
    ARRAY_VIEW<const VECTOR<T,2> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<double,2> >::Get_Edges_Near_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,2> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,2> >&,
    ARRAY_VIEW<const VECTOR<T,2> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Repulsion_Impulse
//#####################################################################
template<class TV> template<class T_PAIR> inline typename TV::SCALAR TRIANGLE_REPULSIONS<TV>::
Repulsion_Impulse(TV& direction,const T dt,const T_PAIR& pair,const TV& relative_velocity,const bool elastic_repulsion,const bool friction)
{
    const T spring_limiter_fraction_over_dt=spring_limiter_fraction/dt;
    const T relative_speed=TV::Dot_Product(relative_velocity,pair.normal);

    // compute scalar impulse
    T scalar_impulse;
    if(elastic_repulsion){
        T total_thickness=pair.Total_Repulsion_Thickness(repulsion_thickness);
        T repulsion_thickness_minus_distance=total_thickness-pair.distance;
        T final_relative_speed=spring_limiter_fraction_over_dt*repulsion_thickness_minus_distance;
        if(relative_speed>=final_relative_speed) return 0;
        T ym_over_mass_times_length=youngs_modulus*geometry.deformable_body_collection.particles.one_over_effective_mass.Subset(pair.nodes).Average()/total_thickness;
        T spring_impulse=dt*ym_over_mass_times_length*max((T)0,repulsion_thickness_minus_distance);
        if(friction) scalar_impulse=min(final_relative_speed-relative_speed,spring_impulse-min((T)0,relative_speed));
        else scalar_impulse=min(final_relative_speed-relative_speed,spring_impulse);}
    else if(relative_speed>=0) return 0;
    else{
        T total_thickness=pair.Total_Repulsion_Thickness(repulsion_thickness);
        T final_relative_speed=total_thickness>pair.distance?0:spring_limiter_fraction_over_dt*(total_thickness-pair.distance);
        if(relative_speed>=final_relative_speed) return 0;
        scalar_impulse=final_relative_speed-relative_speed;}

    // compute friction if necessary and set direction
    if(friction){
        TV tangent=relative_velocity.Projected_Orthogonal_To_Unit_Direction(pair.normal);
        T relative_tangent_velocity_magnitude=tangent.Normalize(),friction_based_velocity_change=friction_coefficient*scalar_impulse;
        scalar_impulse=-1;
        if(friction_based_velocity_change<relative_tangent_velocity_magnitude)
            scalar_impulse=-friction_based_velocity_change/relative_tangent_velocity_magnitude;
        direction=tangent;}
    else direction=pair.normal;

    return scalar_impulse;
}
//#####################################################################
// Function Adjust_Velocity_For_Repulsion
//#####################################################################
template<class TV> template<class T_ARRAY> void TRIANGLE_REPULSIONS<TV>::
Adjust_Velocity_For_Repulsion(const T dt,const T_ARRAY& pairs,const bool elastic_repulsion,const bool friction,const bool use_repulsions,
    ARRAY<TV>& target_impulses,ARRAY<T>& old_speeds,ARRAY<TV>& normals,int inelastic_collision_repulsion_attempts)
{
    DEFORMABLE_PARTICLES<TV>& particles=geometry.deformable_body_collection.particles;
    int attempts=0,total_attempts=1;
    if(!elastic_repulsion && !friction) total_attempts=inelastic_collision_repulsion_attempts;
    int inverted_pairs=0,applied_impulses=0;
    ARRAY_VIEW<TV> V(particles.V);
    ARRAY_VIEW<const T> one_over_effective_mass(particles.one_over_effective_mass);
    while(++attempts<=total_attempts){
        impulse_velocities.Resize(particles.Size());impulse_velocities=V;
        target_impulses.Resize(pairs.Size());target_impulses.Fill(TV());
        normals.Resize(pairs.Size());normals.Fill(TV());
        old_speeds.Resize(pairs.Size());old_speeds.Fill(T());

        for(int pair_index=0;pair_index<pairs.Size();pair_index++){
            const REPULSION_PAIR<TV>& pair=pairs(pair_index);
            if(pair.distance<0) inverted_pairs++;

            TV relative_velocity=-V.Subset(pair.nodes).Weighted_Sum(pair.weights);
            TV direction;T scalar_impulse=Repulsion_Impulse(direction,dt,pair,relative_velocity,elastic_repulsion,friction);
            if(scalar_impulse){
                applied_impulses++;
                T one_over_mass=geometry.deformable_body_collection.binding_list.One_Over_Effective_Mass(pair.nodes,pair.weights);
                TV impulse=Pseudo_Divide(scalar_impulse*direction,one_over_mass);
                if(use_gauss_jacobi && !friction && !elastic_repulsion){
                    for(int i=0;i<TV::m+1;i++) impulse_velocities(pair.nodes(i))-=pair.weights(i)*one_over_effective_mass(pair.nodes(i))*impulse;
                    old_speeds(pair_index)=TV::Dot_Product(relative_velocity,direction);
                    target_impulses(pair_index)=impulse;
                    normals(pair_index)=direction;} // tangential for friction.
                else for(int i=0;i<TV::m+1;i++) V(pair.nodes(i))-=pair.weights(i)*one_over_effective_mass(pair.nodes(i))*impulse;}}

        if(use_gauss_jacobi) Scale_And_Apply_Impulses(pairs,target_impulses,old_speeds,normals);}
    if(inverted_pairs) LOG::Stat("inverted repulsion pairs",inverted_pairs);
    if(applied_impulses) LOG::Stat("applied repulsion impulses",applied_impulses);
}
//#####################################################################
// Function Project_All_Moving_Constraints
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Project_All_Moving_Constraints(const ARRAY<PRECOMPUTE_PROJECT<VECTOR<T,3> > >& point_face_precomputed,
    const ARRAY<PRECOMPUTE_PROJECT<VECTOR<T,3> > >& edge_edge_precomputed,ARRAY_VIEW<VECTOR<T,3> >& field)
{
    for(int i=0;i<point_face_precomputed.m;i++) point_face_precomputed(i).Project(field.Subset(point_face_precomputed(i).nodes));
    for(int i=0;i<edge_edge_precomputed.m;i++) edge_edge_precomputed(i).Project(field.Subset(edge_edge_precomputed(i).nodes));
    for(int i=edge_edge_precomputed.m-2;i>=0;i--) edge_edge_precomputed(i).Project(field.Subset(edge_edge_precomputed(i).nodes));
    for(int i=point_face_precomputed.m-1;i>=0;i--) point_face_precomputed(i).Project(field.Subset(point_face_precomputed(i).nodes));
}
//#####################################################################
// Function Set_Collision_Pairs
//#####################################################################
template<class TV> template<class TV2> void TRIANGLE_REPULSIONS<TV>::
Set_Collision_Pairs(ARRAY<PRECOMPUTE_PROJECT<VECTOR<T,3> > >& point_face_precomputed,
    ARRAY<PRECOMPUTE_PROJECT<VECTOR<T,3> > >& edge_edge_precomputed,ARRAY<REPULSION_PAIR<VECTOR<T,3> > >& point_face_pairs,
    ARRAY<REPULSION_PAIR<TV2> >& edge_edge_pairs,const T repulsion_thickness_multiplier)
{
    STATIC_ASSERT((IS_SAME<VECTOR<T,3>,TV2>::value));
    point_face_pairs.Remove_All();
    edge_edge_pairs.Remove_All();
    point_face_pairs.Append_Elements(point_face_interaction_pairs);
    edge_edge_pairs.Append_Elements(edge_edge_interaction_pairs);
    repulsion_thickness*=repulsion_thickness_multiplier;
    Update_Repulsion_Pairs_Using_History(point_face_pairs,edge_edge_pairs,true); // TODO: Something is wrong here...
    repulsion_thickness/=repulsion_thickness_multiplier;
    point_face_precomputed.Resize(point_face_pairs.m);
    edge_edge_precomputed.Resize(edge_edge_pairs.m);
    for(int i=0;i<point_face_pairs.m;i++){const REPULSION_PAIR<VECTOR<T,3> >& pr=point_face_pairs(i);
        point_face_precomputed(i).Precompute(geometry.deformable_body_collection.particles.one_over_mass.Subset(pr.nodes),pr.weights,pr.normal);}
    for(int i=0;i<edge_edge_pairs.m;i++){const REPULSION_PAIR<VECTOR<T,3> >& pr=edge_edge_pairs(i);
        edge_edge_precomputed(i).Precompute(geometry.deformable_body_collection.particles.one_over_mass.Subset(pr.nodes),pr.weights,pr.normal);}
    internal_point_face_precomputed=point_face_precomputed;
    internal_edge_edge_precomputed=edge_edge_precomputed;
}
//#####################################################################
// Function Scale_And_Apply_Point_Face_Impulses
//#####################################################################
template<class TV> template<class T_ARRAY> void TRIANGLE_REPULSIONS<TV>::
Scale_And_Apply_Impulses(const T_ARRAY& pairs,ARRAY<TV>& target_impulses,ARRAY<T>& old_speeds,ARRAY<TV>& normals)
{
    // Go through the computed impulses, and compare them to the actual change in velocity seen.  Scale back impulses accordingly
    for(int i=0;i<target_impulses.m;i++){
        if(target_impulses(i)==TV()) continue;
        const REPULSION_PAIR<TV>& pair=pairs(i);
        // Compute actual new relative_speed
        TV relative_velocity=-impulse_velocities.Subset(pair.nodes).Weighted_Sum(pair.weights);
        T relative_speed=TV::Dot_Product(relative_velocity,normals(i));
        if(relative_speed*old_speeds(i)<0){
            T new_scale=-(relative_speed-old_speeds(i))/old_speeds(i);
            target_impulses(i)/=new_scale;}} // TODO: should we be doing this, or storing a maximum scale factor at each node?
    // Apply the newly scaled impulses
    ARRAY_VIEW<T>& one_over_effective_mass=geometry.deformable_body_collection.particles.one_over_effective_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    for(int i=0;i<target_impulses.m;i++){
        if(target_impulses(i)==TV()) continue;
        const REPULSION_PAIR<TV>& pair=pairs(i);
        for(int j=0;j<TV::m+1;j++) V(pair.nodes(j))-=pair.weights(j)*one_over_effective_mass(pair.nodes(j))*target_impulses(i);}
}
//#####################################################################
// Precompute
//#####################################################################
template<class TV> void PRECOMPUTE_PROJECT<TV>::
Precompute(const INDIRECT_ARRAY<ARRAY_VIEW<T>,VECTOR<int,TV::m+1>&> one_over_mass,const VECTOR<T,TV::m+1>& weights_input,const TV& normal_input)
{
    // TODO: handle bindings
    weights=weights_input;normal=normal_input;
    T tau=one_over_mass.Weighted_Sum(sqr(weights));
    for(int i=0;i<TV::m+1;i++)
        v_scaled_normals(i)=tau*weights(i)*one_over_mass(i)*normal;
    nodes=one_over_mass.indices;
}
//####################################################################te
template class TRIANGLE_REPULSIONS<VECTOR<float,1> >;
template class TRIANGLE_REPULSIONS<VECTOR<float,2> >;
template class TRIANGLE_REPULSIONS<VECTOR<float,3> >;
template void TRIANGLE_REPULSIONS<VECTOR<float,3> >::Set_Collision_Pairs<VECTOR<float,3> >(ARRAY<PRECOMPUTE_PROJECT<VECTOR<float,3> >,int>&,
    ARRAY<PRECOMPUTE_PROJECT<VECTOR<float,3> >,int>&,ARRAY<REPULSION_PAIR<VECTOR<float,3> >,int>&,
    ARRAY<REPULSION_PAIR<VECTOR<float,3> >,int>&,float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGLE_REPULSIONS<VECTOR<double,1> >;
template class TRIANGLE_REPULSIONS<VECTOR<double,2> >;
template class TRIANGLE_REPULSIONS<VECTOR<double,3> >;
template void TRIANGLE_REPULSIONS<VECTOR<double,3> >::Set_Collision_Pairs<VECTOR<double,3> >(ARRAY<PRECOMPUTE_PROJECT<VECTOR<double,3> >,int>&,
    ARRAY<PRECOMPUTE_PROJECT<VECTOR<double,3> >,int>&,ARRAY<REPULSION_PAIR<VECTOR<double,3> >,int>&,
    ARRAY<REPULSION_PAIR<VECTOR<double,3> >,int>&,double);
#endif
}
