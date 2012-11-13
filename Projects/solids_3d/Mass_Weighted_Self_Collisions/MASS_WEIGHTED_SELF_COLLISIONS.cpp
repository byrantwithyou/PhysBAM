//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include "MASS_WEIGHTED_SELF_COLLISIONS.h"
using namespace PhysBAM;
//#####################################################################
// Function MASS_WEIGHTED_SELF_COLLISIONS
//#####################################################################
template<class T_input> MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
MASS_WEIGHTED_SELF_COLLISIONS(const STREAM_TYPE stream_type)
    :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection)
{
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Register_Options()
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Parse_Options()
{
    BASE::Parse_Options();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Initialize_Bodies()
{
    //helper references
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    ARRAY<float> masses;

    output_directory=STRING_UTILITIES::string_sprintf("Mass_Weighted_Self_Collisions/Test_%d",test_number);
    
    last_frame=1000;
    frame_rate=24;
    solids_parameters.enforce_repulsions_in_cg=false;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T)1e-4;
    solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
    solids_parameters.triangle_collision_parameters.total_collision_loops=4;
    comparator=new COLLISION_PAIR_COMPARATOR(&solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry);
    
    switch(test_number){
        case 1:{
            particles.Add_Elements(4);
            particles.X(0)=TV(0,0,0);particles.mass(0)=(T).00333333;particles.V(0)=TV(0,-(T).1,0);
            particles.X(1)=TV((T)0.01,0,0);particles.mass(1)=(T).00333333;particles.V(1)=TV(0,-(T).1,0);
            particles.X(2)=TV(0,0,(T)0.01);particles.mass(2)=(T).00333333;particles.V(2)=TV(0,-(T).1,0);
            //particles.X(3)=TV((T).0025,-(T).01,(T).0025);// THIS IS SET BELOW
            particles.mass(3)=(T).01;particles.V(3)=TV(0,(T).1,0); 
            TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);
            surface.mesh.elements.Append(VECTOR<int,3>(0,1,2));
            FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create(particles);
            particles.X(3)=surface.Get_Element(0).Center()+TV(0,(T)-.01,0);
            free_particles.nodes.Append(3);
            deformable_body_collection.deformable_geometry.Add_Structure(&surface);
            deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);
        }break;
        case 2:{
            particles.Add_Elements(4);
            particles.X(0)=TV(0,0,(T).005);particles.mass(0)=(T).00333333;particles.V(0)=TV(0,(T).1,0);
            particles.X(1)=TV((T)0.01,0,(T).005);particles.mass(1)=(T).00333333;particles.V(1)=TV(0,(T).1,0);
            particles.X(2)=TV((T).005,(T).01,0);particles.mass(2)=(T).00333333;particles.V(2)=TV(0,-(T).1,0);
            particles.X(3)=TV((T).005,(T).01,(T).01);particles.mass(3)=(T).00333333;particles.V(3)=TV(0,-(T).1,0); 
            SEGMENTED_CURVE<TV>& curve=*SEGMENTED_CURVE<TV>::Create(particles);
            curve.mesh.elements.Append(VECTOR<int,2>(0,1));
            curve.mesh.elements.Append(VECTOR<int,2>(2,3));
            deformable_body_collection.deformable_geometry.Add_Structure(&curve);
        }break;
        case 3:{
            particles.Add_Elements(4);
            particles.X(0)=TV((T).1,(T).1,0);particles.mass(0)=(T).00333333;particles.V(0)=TV(0,0,0);
            particles.X(1)=TV(0,(T).1,0);particles.mass(1)=(T).00333333;particles.V(1)=TV(0,0,0);
            particles.X(2)=TV((T).1,(T).15,(T).1);particles.mass(2)=(T).00333333;particles.V(2)=TV(0,0,-(T).1);
            particles.X(3)=TV(0,(T).09,(T).1);particles.mass(3)=(T).00333333;particles.V(3)=TV(0,0,-(T).1); 
            SEGMENTED_CURVE<TV>& curve=*SEGMENTED_CURVE<TV>::Create(particles);
            curve.mesh.elements.Append(VECTOR<int,2>(0,1));
            curve.mesh.elements.Append(VECTOR<int,2>(2,3));
            deformable_body_collection.deformable_geometry.Add_Structure(&curve);
        }break;
        case 4:{
            int max=5;
            particles.Add_Elements(4*max);
            for (int i=0;i<max;i++) {
                particles.X(4*i)=TV(0,(T).02*i,0);particles.mass(4*i)=(T).00333333;particles.V(4*i)=TV(0,-(T).1,0);
                particles.X(4*i+1)=TV((T)0.01,(T).02*i,0);particles.mass(4*i+1)=(T).00333333;particles.V(4*i+1)=TV(0,-(T).1,0);
                particles.X(4*i+2)=TV(0,(T).02*i,(T)0.01);particles.mass(4*i+2)=(T).00333333;particles.V(4*i+2)=TV(0,-(T).1,0);
                particles.mass(4*i+3)=(T).01;particles.V(4*i+3)=TV(0,-(T).1,0); 
            }
            particles.V(0)=TV(0,(T).1,0);particles.V(1)=TV(0,(T).1,0);particles.V(2)=TV(0,(T).1,0);
            TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);
            for (int i=0;i<max;i++) surface.mesh.elements.Append(VECTOR<int,3>(4*i,4*i+1,4*i+2));
            FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create(particles);
            for (int i=0;i<max;i++) {
                particles.X(4*i+3)=surface.Get_Element(i+1).Center()+TV(0,(T).01,0);
                free_particles.nodes.Append(4*i+3);
            }
            deformable_body_collection.deformable_geometry.Add_Structure(&surface);
            deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);
        }break;
    }


    //################################################################
    // Mass Update Phase
    //################################################################
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();



    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.collision_body_list.Add_Bodies(*solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list);

    // number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();
    
    // Forces
    switch(test_number){
        case 1:{
            TRIANGULATED_SURFACE<T>& surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>(0);
            //solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,false));
            solid_body_collection.Add_Force(Create_Edge_Springs(surface,(T)1000,(T)2)); // were *2 and *10
        }break;
        case 2:{
            SEGMENTED_CURVE<TV>& curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>&>(0);
            solid_body_collection.Add_Force(Create_Edge_Springs(curve,(T)1000,(T)2)); // were *2 and *10
        }break;
        case 3:{
            SEGMENTED_CURVE<TV>& curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>&>(0);
            solid_body_collection.Add_Force(Create_Edge_Springs(curve,(T)1000,(T)2)); // were *2 and *10
        }break;
        case 4:{
            TRIANGULATED_SURFACE<T>& surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>(0);
            //solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,false));
            solid_body_collection.Add_Force(Create_Edge_Springs(surface,(T)1000,(T)2)); // were *2 and *10
        }break;
    }
}
//#####################################################################
// Function Point_Face_Mass
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,VECTOR<T,4>& one_over_mass)
{
    T factor = (T)0.5-sqr(attempt_ratio)/(T)2.;    
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry=solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry;
    const ARRAY<TV>& X=geometry.X_self_collision_free;
    VECTOR<T,3> X_embedded=weights[0]*X(nodes[1])+weights[1]*X(nodes[2])+weights[2]*X(nodes[3]);
    if(X_embedded.y<X(nodes[0]).y) {one_over_mass[1]*=factor;one_over_mass[2]*=factor;one_over_mass[3]*=factor;}
    else one_over_mass[0]*=factor;
}
//#####################################################################
// Function Point_Face_Mass
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,ARRAY_VIEW<T>& one_over_mass)
{
    T factor = (T)0.5-sqr(attempt_ratio)/(T)2.;    
    saved_mass=one_over_mass.Subset(nodes);
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry=solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry;
    const ARRAY<TV>& X=geometry.X_self_collision_free;
    VECTOR<T,3> X_embedded=weights[0]*X(nodes[1])+weights[1]*X(nodes[2])+weights[2]*X(nodes[3]);
    if(X_embedded.y<X(nodes[0]).y) {one_over_mass(nodes[1])*=factor;one_over_mass(nodes[2])*=factor;one_over_mass(nodes[3])*=factor;}
    else one_over_mass(nodes[0])*=factor;
}
//#####################################################################
// Function Point_Face_Mass_Revert
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Point_Face_Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass)
{one_over_mass(nodes[0])=saved_mass[0];one_over_mass(nodes[1])=saved_mass[1];one_over_mass(nodes[2])=saved_mass[2];one_over_mass(nodes[3])=saved_mass[3];}
//#####################################################################
// Function Edge_Edge_Mass
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,VECTOR<T,4>& one_over_mass)
{
    T factor = (T)0.5-sqr(attempt_ratio)/(T)2.;    
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry=solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry;
    const ARRAY<TV>& X=geometry.X_self_collision_free;
    VECTOR<T,3> X_embedded1=(1-weights[0])*X(nodes[0])+weights[0]*X(nodes[1]),X_embedded2=(1-weights[1])*X(nodes[2])+weights[1]*X(nodes[3]);
    if(X_embedded1.y<X_embedded2.y){one_over_mass[0]*=factor;one_over_mass[1]*=factor;}
    else{one_over_mass[2]*=factor;one_over_mass[3]*=factor;}
}
//#####################################################################
// Function Edge_Edge_Mass
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,ARRAY_VIEW<T>& one_over_mass)
{
    T factor = (T)0.5-sqr(attempt_ratio)/(T)2.;    
    saved_mass=one_over_mass.Subset(nodes);
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry=solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry;
    const ARRAY<TV>& X=geometry.X_self_collision_free;
    VECTOR<T,3> X_embedded1=(1-weights[0])*X(nodes[0])+weights[0]*X(nodes[1]),X_embedded2=(1-weights[1])*X(nodes[2])+weights[1]*X(nodes[3]);
    if(X_embedded1.y<X_embedded2.y){one_over_mass(nodes[0])*=factor;one_over_mass(nodes[1])*=factor;}
    else{one_over_mass(nodes[2])*=factor;one_over_mass(nodes[3])*=factor;}
}
//#####################################################################
// Function Edge_Edge_Mass_Revert
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Edge_Edge_Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass)
{one_over_mass(nodes[0])=saved_mass[0];one_over_mass(nodes[1])=saved_mass[1];one_over_mass(nodes[2])=saved_mass[2];one_over_mass(nodes[3])=saved_mass[3];}
//#####################################################################
// Function Reorder_Pairs
//#####################################################################
template<class T_input> void MASS_WEIGHTED_SELF_COLLISIONS<T_input>::
Reorder_Pairs(ARRAY<VECTOR<int,4> >& edge_edge_pairs,ARRAY<VECTOR<int,4> >& point_face_pairs) {
    edge_edge_pairs.Sort(*comparator);
    point_face_pairs.Sort(*comparator);
}
//#####################################################################
template class MASS_WEIGHTED_SELF_COLLISIONS<float>;
template class MASS_WEIGHTED_SELF_COLLISIONS<double>;
