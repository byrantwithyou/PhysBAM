//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINDING_LIST
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/PARTICLE_BINDING.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BINDING_LIST<TV>::
BINDING_LIST(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection)
    :particles(deformable_body_collection.particles),deformable_body_collection(deformable_body_collection)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BINDING_LIST<TV>::
~BINDING_LIST()
{
    bindings.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clean_Memory()
{
    bindings.Delete_Pointers_And_Clean_Memory();
    binding_index_from_particle_index.Clean_Memory();
}
//#####################################################################
// Function Add_Binding
//#####################################################################
template<class TV> int BINDING_LIST<TV>::
Add_Binding(BINDING<TV>* binding)
{
    int id=bindings.Append(binding);
    if(binding_index_from_particle_index.m<=binding->particle_index)
        binding_index_from_particle_index.Resize(binding->particle_index+1,use_init,-1);
    binding_index_from_particle_index(binding->particle_index)=id;
    return id;
}
//#####################################################################
// Function Update_Binding_Index_From_Particle_Index
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Update_Binding_Index_From_Particle_Index()
{
    binding_index_from_particle_index.Clean_Memory();
    int max_particle_index=-1;for(int b=0;b<bindings.m;b++) max_particle_index=max(max_particle_index,bindings(b)->particle_index);
    binding_index_from_particle_index.Resize(max_particle_index+1,use_init,-1);
    for(int b=0;b<bindings.m;b++){assert(binding_index_from_particle_index(bindings(b)->particle_index)<0);binding_index_from_particle_index(bindings(b)->particle_index)=b;}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int b=0;b<bindings.m;b++) bindings(b)->Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Compute_Dependency_Closure_Based_On_Embedding
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Compute_Dependency_Closure_Based_On_Embedding(SEGMENT_MESH& dependency_mesh) const
{
    ARRAY<VECTOR<int,2> > dependencies;dependency_mesh.elements.Exchange(dependencies);dependency_mesh.Refresh_Auxiliary_Structures();
    for(int i=0;i<dependencies.m;i++){VECTOR<int,2>& dependency=dependencies(i);
        ARRAY<int> parents1,parents2;
        Parents(parents1,dependency.x);
        Parents(parents2,dependency.y);
        for(int i=0;i<parents1.m;i++)
            for(int j=0;j<parents2.m;j++)
                dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(parents1(i),parents2(j)));}
}
//#####################################################################
// Function Compute_Particle_Closure_Based_On_Embedding
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Compute_Particle_Closure_Based_On_Embedding(ARRAY<int>& particle_set) const
{
    ARRAY<bool> particle_is_present(particles.Size());
    particle_is_present.Subset(particle_set).Fill(true);
    for(int b=0;b<bindings.m;b++){BINDING<TV>& binding=*bindings(b);
        ARRAY<int> parents;
        binding.Parents(parents);
        if(!particle_is_present.Subset(parents).Contains(false) && !particle_is_present(binding.particle_index))
            particle_set.Append(binding.particle_index);}
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Positions
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Positions() const
{
    for(int i=0;i<bindings.m;i++) bindings(i)->Clamp_To_Embedded_Position();
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Velocities() const
{
    for(int i=0;i<bindings.m;i++) bindings(i)->Clamp_To_Embedded_Velocity();
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Positions
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Positions(ARRAY_VIEW<TV> X) const
{
    for(int i=0;i<bindings.m;i++) X(bindings(i)->particle_index)=bindings(i)->Embedded_Position(X);
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Velocities(ARRAY_VIEW<TV> V) const
{
    for(int i=0;i<bindings.m;i++) V(bindings(i)->particle_index)=bindings(i)->Embedded_Velocity(V);
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Velocities(ARRAY_VIEW<TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const
{
    for(int i=0;i<bindings.m;i++) V(bindings(i)->particle_index)=bindings(i)->Embedded_Velocity(V,twist);
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full) const
{
    for(int i=0;i<bindings.m;i++) bindings(i)->Distribute_Force_To_Parents(F_full,F_full(bindings(i)->particle_index));
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full) const
{
    for(int i=0;i<bindings.m;i++) bindings(i)->Distribute_Force_To_Parents(F_full,wrench_full,F_full(bindings(i)->particle_index));
}
//#####################################################################
// Function Distribute_Mass_To_Parents
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Distribute_Mass_To_Parents() const
{
    for(int b=0;b<bindings.m;b++) bindings(b)->Distribute_Mass_To_Parents(particles.mass);
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Read(TYPED_ISTREAM input)
{
    Clean_Memory();
    int m;Read_Binary(input,m);bindings.Preallocate(m);
    for(int k=0;k<m;k++){BINDING<TV>* binding=BINDING<TV>::Create(input,particles);
        Add_Binding(binding);}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Write(TYPED_OSTREAM output) const
{
    Write_Binary(output,bindings.m);
    for(int k=0;k<bindings.m;k++) bindings(k)->Write(output);
}
//#####################################################################
// Function Update_Neighbor_Bindings
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Update_Neighbor_Bindings()
{
    ARRAY<ARRAY<int> > children(particles.number);
    ARRAY<int> parents;
    for(int k=0;k<bindings.m;k++){
        parents.Remove_All();
        bindings(k)->Parents(parents);
        for(int i=0;i<parents.m;i++)
            children(parents(i)).Append(bindings(k)->particle_index);}

    neighbor_bindings.Resize(particles.number);
    for(int k=0;k<bindings.m;k++){
        parents.Remove_All();
        bindings(k)->Parents(parents);
        for(int i=0;i<parents.m;i++)
            neighbor_bindings(bindings(k)->particle_index).Append_Elements(children(parents(i)));
        Prune_Duplicates(neighbor_bindings(bindings(k)->particle_index));}
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Apply_Impulse(const int particle_index,const TV& impulse,bool update_neighbors)
{
    if(BINDING<TV>* binding=Binding(particle_index)){
        binding->Apply_Impulse(impulse);
        if(update_neighbors)
            for(int i=0;i<neighbor_bindings(binding->particle_index).m;i++)
                particles.V(neighbor_bindings(binding->particle_index)(i))=V(neighbor_bindings(binding->particle_index)(i));}
    else particles.V(particle_index)+=particles.one_over_mass(particle_index)*impulse;
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Apply_Impulse(const int particle_index,const TV& impulse,ARRAY_VIEW<TV> V_input,bool update_neighbors) const
{
    if(BINDING<TV>* binding=Binding(particle_index)){
        binding->Apply_Impulse(impulse,V_input);
        if(update_neighbors)
            for(int i=0;i<neighbor_bindings(binding->particle_index).m;i++)
                V_input(neighbor_bindings(binding->particle_index)(i))=V(neighbor_bindings(binding->particle_index)(i),V_input);}
    else V_input(particle_index)+=particles.one_over_mass(particle_index)*impulse;
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Apply_Impulse(const int particle_index,const TV& impulse,ARRAY_VIEW<TV> V_input,ARRAY_VIEW<TWIST<TV> > rigid_V,bool update_neighbors) const
{
    if(BINDING<TV>* binding=Binding(particle_index)){
        binding->Apply_Impulse(impulse,V_input,rigid_V);
        if(update_neighbors)
            for(int i=0;i<neighbor_bindings(binding->particle_index).m;i++)
                V_input(neighbor_bindings(binding->particle_index)(i))=V(neighbor_bindings(binding->particle_index)(i),V_input,rigid_V);}
    else V_input(particle_index)+=particles.one_over_mass(particle_index)*impulse;
}
//#####################################################################
// Function Remove_Bindings
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Remove_Bindings(const ARRAY<int>& bound_particles)
{
    ARRAY<bool> remove_binding(bindings.m);
    remove_binding.Subset(binding_index_from_particle_index.Subset(bound_particles)).Fill(true);

    // recreate the binding list without the input binding
    ARRAY<BINDING<TV>*> bindings_old=bindings;
    bindings.Remove_All();
    binding_index_from_particle_index.Remove_All();

    // add bindings that aren't being removed
    for(int b=0;b<bindings_old.m;b++)
    {
        if(remove_binding(b)) delete bindings_old(b);
        else Add_Binding(bindings_old(b));
    }
}
//#####################################################################
namespace PhysBAM{
template class BINDING_LIST<VECTOR<float,1> >;
template class BINDING_LIST<VECTOR<float,2> >;
template class BINDING_LIST<VECTOR<float,3> >;
template class BINDING_LIST<VECTOR<double,1> >;
template class BINDING_LIST<VECTOR<double,2> >;
template class BINDING_LIST<VECTOR<double,3> >;
}
