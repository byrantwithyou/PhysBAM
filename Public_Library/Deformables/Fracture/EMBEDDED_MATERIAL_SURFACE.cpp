//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_MATERIAL_SURFACE
//##################################################################### 
#include <Core/Log/LOG.h>
#include <Deformables/Bindings/LINEAR_BINDING.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Fracture/EMBEDDED_MATERIAL_SURFACE.h>
#include <Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <Deformables/Fracture/TRIANGLES_OF_MATERIAL.h>
#include <Deformables/Particles/FREE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> EMBEDDED_MATERIAL_SURFACE<TV,d>::
EMBEDDED_MATERIAL_SURFACE(T_EMBEDDED_OBJECT& embedded_object)
    :EMBEDDING<TV>(embedded_object.particles),embedded_object(embedded_object),embedded_particles(embedded_object.embedded_particles),
    parent_particles(embedded_object.parent_particles),interpolation_fraction(embedded_object.interpolation_fraction),number_of_soft_bound_particles(0),
    need_destroy_embedded_object(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> EMBEDDED_MATERIAL_SURFACE<TV,d>::
~EMBEDDED_MATERIAL_SURFACE()
{
    if(need_destroy_embedded_object) delete &embedded_object;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> typename EMBEDDING_POLICY<TV,d>::EMBEDDING* EMBEDDED_MATERIAL_SURFACE<TV,d>::
Create()
{
    typename EMBEDDING_POLICY<TV,d>::EMBEDDING* embedding=new typename EMBEDDING_POLICY<TV,d>::EMBEDDING(*T_EMBEDDED_OBJECT::Create());
    embedding->need_destroy_embedded_object=true;return embedding;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> typename EMBEDDING_POLICY<TV,d>::EMBEDDING* EMBEDDED_MATERIAL_SURFACE<TV,d>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    typename EMBEDDING_POLICY<TV,d>::EMBEDDING* embedding=new typename EMBEDDING_POLICY<TV,d>::EMBEDDING(*T_EMBEDDED_OBJECT::Create(particles));
    embedding->need_destroy_embedded_object=true;return embedding;
}
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV,int d> STRUCTURE<TV>* EMBEDDED_MATERIAL_SURFACE<TV,d>::
Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const
{
    typename EMBEDDING_POLICY<TV,d>::EMBEDDING* embedding=Create(new_particles);
    int offset=new_particles.Size();
    new_particles.Append(particles);
    if(particle_indices) for(int p=0;p<particles.Size();p++) particle_indices->Append(p+offset);
    embedding->embedded_object.embedded_particles.active_indices=embedded_object.embedded_particles.active_indices;
    embedding->embedded_object.embedded_particles.active_indices+=offset;
    embedding->embedded_object.embedded_particles.Update_Subset_Index_From_Element_Index();
    embedding->embedded_object.parent_particles=embedded_object.parent_particles+offset;
    embedding->embedded_object.interpolation_fraction=embedded_object.interpolation_fraction; 
    embedding->embedded_object.interpolation_fraction_threshold=embedded_object.interpolation_fraction_threshold; 
    embedding->embedded_object.average_interpolation_fractions=embedded_object.average_interpolation_fractions;
    embedding->embedded_object.simplicial_object.mesh.Initialize_Mesh_With_Particle_Offset(embedded_object.simplicial_object.mesh,offset);
    embedding->embedded_object.embedded_mesh.Initialize_Mesh_With_Particle_Offset(embedded_object.embedded_mesh,offset);
    embedding->embedded_object.node_in_simplex_is_material=embedded_object.node_in_simplex_is_material;
    embedding->embedded_object.orientation_index=embedded_object.orientation_index;
    embedding->material_surface_mesh.Initialize_Mesh_With_Particle_Offset(material_surface_mesh,offset);
    embedding->previously_perturbed=previously_perturbed;
    embedding->Update_Number_Nodes();
    return embedding;
}
//#####################################################################
// Function Create_Material_Surface
//#####################################################################
template<class TV,int d> void EMBEDDED_MATERIAL_SURFACE<TV,d>::
Create_Material_Surface(const bool verbose)
{
    material_surface_mesh.Clean_Memory();material_surface.Clean_Memory();
    Construct_Material_Surface_Mesh();
    material_surface_mesh.number_nodes=embedded_object.simplicial_object.mesh.number_nodes;
    if(verbose) LOG::cout<<"Triangles Of Material Surface Mesh: "<<material_surface_mesh.elements.m<<std::endl;
    previously_perturbed.Resize(embedded_object.embedded_particles.active_indices.m);
}
//#####################################################################
// Function Update_Binding_List_And_Particles_From_Embedding
//#####################################################################
template<class TV,int d> void EMBEDDED_MATERIAL_SURFACE<TV,d>::
Update_Binding_List_From_Embedding(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection,bool add_soft_bindings)
{
    BINDING_LIST<TV>& binding_list=deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=deformable_body_collection.soft_bindings;
    binding_list.Clean_Memory();soft_bindings.Clean_Memory();
    FREE_PARTICLES<TV>* free_particles=deformable_body_collection.template Find_Structure<FREE_PARTICLES<TV>*>();
    if(!free_particles){
        free_particles=new FREE_PARTICLES<TV>;
        deformable_body_collection.Add_Structure(free_particles);
        deformable_body_collection.collisions.collision_structures.Append(free_particles);}
    else free_particles->nodes.Remove_All();
    for(int p=0;p<embedded_particles.active_indices.m;p++) binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(particles),embedded_particles.active_indices(p),
        parent_particles(p),VECTOR<T,2>((T)1-interpolation_fraction(p),interpolation_fraction(p))));
    if(add_soft_bindings){
        ARRAY<bool> used(deformable_body_collection.particles.Size());
        used.Subset(material_surface_mesh.elements.Flattened()).Fill(true);
        int j=-1;
        for(int p=0;p<embedded_particles.active_indices.m;p++){
            if(used(embedded_particles.active_indices(p))){
                if(++j>=soft_particles.m) soft_particles.Append(deformable_body_collection.particles.Add_Element());
                soft_bindings.Add_Binding(VECTOR<int,2>(soft_particles(j),embedded_particles.active_indices(p)),true);
                free_particles->nodes.Append(soft_particles(j));}}
        number_of_soft_bound_particles=embedded_particles.active_indices.m;
        soft_bindings.Clamp_Particles_To_Embedded_Positions();soft_bindings.Clamp_Particles_To_Embedded_Velocities();
        soft_bindings.particles.Compute_Auxiliary_Attributes(soft_bindings);
        soft_bindings.Set_Mass_From_Effective_Mass();}
}
//#####################################################################
namespace PhysBAM{
template class EMBEDDED_MATERIAL_SURFACE<VECTOR<float,2>,2>;
template class EMBEDDED_MATERIAL_SURFACE<VECTOR<float,3>,2>;
template class EMBEDDED_MATERIAL_SURFACE<VECTOR<float,3>,3>;
template class EMBEDDED_MATERIAL_SURFACE<VECTOR<double,2>,2>;
template class EMBEDDED_MATERIAL_SURFACE<VECTOR<double,3>,2>;
template class EMBEDDED_MATERIAL_SURFACE<VECTOR<double,3>,3>;
}
