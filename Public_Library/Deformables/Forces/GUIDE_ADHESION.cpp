//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GUIDE_ADHESION
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/QUEUE.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Deformables/Forces/GUIDE_ADHESION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GUIDE_ADHESION<TV>::
GUIDE_ADHESION(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& mesh_input,SEGMENT_MESH& guide_mesh_input,
    ARRAY<HAIR_ID>& particle_to_spring_id_input,ARRAY<int,HAIR_ID>& roots_input)
    :DEFORMABLES_FORCES<TV>(particles),mesh(mesh_input),guide_mesh(guide_mesh_input),curve(mesh,particles),guide_curve(guide_mesh,particles),
    thickness(10),max_connections(2000),particle_to_spring_id(particle_to_spring_id_input),roots(roots_input),one_over_effective_mass(particles.one_over_effective_mass),
    springs(new T_SPRING_HASH()),segments_with_springs(guide_mesh.elements.m)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GUIDE_ADHESION<TV>::
~GUIDE_ADHESION()
{
    delete springs;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    this->cfl_initialized=false;
}
//#####################################################################
// Function Set_Parameters
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Set_Parameters(const T youngs_modulus_input,const T overdamping_fraction_input,const T thickness_input,const int max_connections_input)
{
    this->youngs_modulus=youngs_modulus_input;
    this->overdamping_fraction=overdamping_fraction_input;
    this->thickness=thickness_input;
    this->max_connections=max_connections_input;
}
//#####################################################################
// Function Update_Springs
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Update_Springs(const bool search_hierarchy)
{
    // prune broken springs
    ARRAY<VECTOR<int,2> > deletion_list;
    T_SPRING_HASH* new_springs=search_hierarchy?new T_SPRING_HASH():0; // many deletes can hurt open addressed hash table, use "generational" garbage collection
    for(typename HASHTABLE<VECTOR<int,2>,SPRING_STATE>::ITERATOR i(*springs);i.Valid();i.Next()){
        const VECTOR<int,2> &segment1_nodes=guide_curve.mesh.elements(i.Key()[0]),&segment2_nodes=curve.mesh.elements(i.Key()[1]);
        SPRING_STATE& state=i.Data();
        SEGMENT_3D<T> segment1(guide_curve.particles.X.Subset(segment1_nodes)),segment2(curve.particles.X.Subset(segment2_nodes));
        TV X0=(1-state.weights[0])*segment1.X.x+state.weights[0]*segment1.X.y;
        TV X1=(1-state.weights[1])*segment2.X.x+state.weights[1]*segment2.X.y;
        state.normal=X0-X1; 
        state.distance=state.normal.Normalize();
        if(state.distance>thickness){
            deletion_list.Append(i.Key());
            segments_with_springs(i.Key()[0])--;}}
    //if(!new_springs) for(int i=0;i<deletion_list.m;i++)springs->Delete(deletion_list(i));

    if(search_hierarchy){ 
        // Initialize/Update Hierarchy
        if(!guide_curve.hierarchy) guide_curve.Initialize_Hierarchy(false);
        guide_curve.hierarchy->Update_Leaf_Boxes(guide_curve.particles.X);
        guide_curve.hierarchy->Update_Nonleaf_Boxes();

        // find new pairs
        for(int i=0;i<curve.mesh.elements.m;i++){
            const VECTOR<int,2> &segment_nodes=curve.mesh.elements(i);
            SEGMENT_3D<T> segment=curve.particles.X.Subset(segment_nodes);
            RANGE<TV> box(segment.X.x,segment.X.y);
            ARRAY<int> intersections;
            guide_curve.hierarchy->Intersection_List(box,intersections,2*thickness);
            SPRING_STATE state;
            int min_index=-2;
            VECTOR<T,2> weights;
            TV normal;
            T distance=0;
            for(int j=0;j<intersections.m;j++){
                SEGMENT_3D<T> guide_segment=guide_curve.particles.X.Subset(guide_curve.mesh.elements(intersections(j)));
                normal=segment.Shortest_Vector_Between_Segments(guide_segment,weights);
                if((distance>normal.Magnitude()||min_index==-1)&&segments_with_springs(intersections(j))<2*max_connections){
                    min_index=intersections(j);
                    distance=normal.Magnitude();
                    state.weights=weights;
                }
            }
            if(min_index==-1) continue;
            const VECTOR<int,2> &segment1_nodes=guide_curve.mesh.elements(min_index),&segment2_nodes=curve.mesh.elements(i);
            SEGMENT_3D<T> segment1(guide_curve.particles.X.Subset(segment1_nodes)),segment2(curve.particles.X.Subset(segment2_nodes));
            TV X0=(1-state.weights[0])*segment1.X.x+state.weights[0]*segment1.X.y;
            TV X1=(1-state.weights[1])*segment2.X.x+state.weights[1]*segment2.X.y;
            assert(particle_to_spring_id(segment1_nodes[0])==particle_to_spring_id(segment1_nodes[1]));
            assert(particle_to_spring_id(segment2_nodes[0])==particle_to_spring_id(segment2_nodes[1]));
            TV guide_root=guide_curve.particles.X(roots(particle_to_spring_id(segment1_nodes[0])));
            TV curve_root=curve.particles.X(roots(particle_to_spring_id(segment2_nodes[0])));
            state.normal=X0-X1;
            state.distance=state.normal.Normalize();
            state.restlength=(guide_root-curve_root).Magnitude();
            state.nodes=guide_curve.mesh.elements(min_index).Append_Elements(segment_nodes);
            T harmonic_mass=Pseudo_Inverse((T).25*(one_over_effective_mass(state.nodes[0])+one_over_effective_mass(state.nodes[1])
                +one_over_effective_mass(state.nodes[2])+one_over_effective_mass(state.nodes[3])));
            state.damping=overdamping_fraction*2*sqrt(youngs_modulus*state.restlength*harmonic_mass);
            new_springs->Insert(VECTOR<int,2>(min_index,i),state);
            segments_with_springs(min_index)++;
        }
        delete springs;
        springs=new_springs;}
}
//#####################################################################
// Function Update_Springs
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    Update_Springs(false);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(const auto& i:*springs){
        const SPRING_STATE& state=i.data;
        TV force=youngs_modulus*(state.distance/state.restlength-(T)1)*state.normal;
        F(state.nodes[0])-=((T)1-state.weights[0])*force;F(state.nodes[1])-=state.weights[0]*force;
        F(state.nodes[2])+=((T)1-state.weights[1])*force;F(state.nodes[3])+=state.weights[1]*force;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(const auto& i:*springs){
        const SPRING_STATE& state=i.data;
        TV force=state.damping/state.restlength*TV::Dot_Product(((T)1-state.weights[0])*V(state.nodes[0])+state.weights[0]*V(state.nodes[1])
                                       -((T)1-state.weights[1])*V(state.nodes[2])-state.weights[1]*V(state.nodes[3]),state.normal)*state.normal;
        F(state.nodes[0])-=((T)1-state.weights[0])*force;F(state.nodes[1])-=state.weights[0]*force;
        F(state.nodes[2])+=((T)1-state.weights[1])*force;F(state.nodes[3])+=state.weights[1]*force;}
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR GUIDE_ADHESION<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{   
    // TODO: Implement
}
//#####################################################################
// Function Write_To_File
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Write_State(STREAM_TYPE stream_type,const std::string& filename)
{   
    Write_To_File(stream_type,filename,*springs);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void GUIDE_ADHESION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
namespace PhysBAM{
template class GUIDE_ADHESION<VECTOR<float,3> >;
template class GUIDE_ADHESION<VECTOR<double,3> >;
}
