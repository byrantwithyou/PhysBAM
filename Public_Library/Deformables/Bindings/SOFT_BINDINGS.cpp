//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOFT_BINDINGS
//#####################################################################
#include <Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <Tools/Log/LOG.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOFT_BINDINGS<TV>::
SOFT_BINDINGS(BINDING_LIST<TV>& binding_list_input)
    :binding_list(binding_list_input),particles(binding_list_input.particles),binding_mesh(0),last_read(-1),is_stale(true),frame_list(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOFT_BINDINGS<TV>::
~SOFT_BINDINGS()
{
    delete binding_mesh;
}
//#####################################################################
// Function Initialize_Binding_Mesh
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Initialize_Binding_Mesh(const bool exclude_particles_using_impulses_for_collisions)
{
    delete binding_mesh;
    if(exclude_particles_using_impulses_for_collisions){
        binding_mesh=new SEGMENT_MESH;binding_mesh->Set_Number_Nodes(particles.Size());
        for(int b=0;b<bindings.m;b++) if(!use_impulses_for_collisions(b)) binding_mesh->elements.Append(bindings(b));}
    else binding_mesh=new SEGMENT_MESH(particles.Size(),bindings);
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int i=0;i<bindings.m;i++) dependency_mesh.Add_Element_If_Not_Already_There(bindings(i));
}
//#####################################################################
// Function Set_Mass_From_Effective_Mass
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Set_Mass_From_Effective_Mass()
{
    for(int b=0;b<bindings.m;b++){
        particles.mass(bindings(b).x)=particles.effective_mass(bindings(b).x);
        particles.one_over_mass(bindings(b).x)=particles.one_over_effective_mass(bindings(b).x);}
}
//#####################################################################
// Function Remove_Soft_Bound_Particles
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Remove_Soft_Bound_Particles(ARRAY<int>& particles) const
{
    for(int i=particles.m-1;i>=0;i--) if(Particle_Is_Bound(particles(i))) particles.Remove_Index_Lazy(i);
}
//#####################################################################
// Function Need_Bindings_Mapped
//#####################################################################
template<class TV> bool SOFT_BINDINGS<TV>::
Need_Bindings_Mapped() const
{
    return bindings_using_forces_for_collisions.m!=0;
}
//#####################################################################
// Function Map_Forces_From_Parents
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Map_Forces_From_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<const TWIST<TV> > wrench_full) const
{
    for(int b=0;b<bindings.m;b++){ // TODO: MPI
        BINDING<TV>* hard_binding=binding_list.Binding(bindings(b).y);
        if(hard_binding) F_full(bindings(b).x)+=particles.mass(bindings(b).x)*hard_binding->Embedded_Acceleration(F_full,wrench_full);
        else F_full(bindings(b).x)+=particles.mass(bindings(b).x)*particles.one_over_mass(bindings(b).y)*F_full(bindings(b).y);}
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Positions
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Clamp_Particles_To_Embedded_Positions(const bool bindings_using_impulses_for_collisions_only) const
{
    if(bindings_using_impulses_for_collisions_only){
        for(int i=0;i<bindings_using_impulses_for_collisions.m;i++){
            const VECTOR<int,2>& binding=bindings(bindings_using_impulses_for_collisions(i));particles.X(binding.x)=particles.X(binding.y);}}
    else for(int i=0;i<bindings.m;i++){particles.X(bindings(i).x)=particles.X(bindings(i).y);}
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Clamp_Particles_To_Embedded_Velocities(const bool bindings_using_impulses_for_collisions_only) const
{
    if(bindings_using_impulses_for_collisions_only){
        for(int i=0;i<bindings_using_impulses_for_collisions.m;i++){
            const VECTOR<int,2>& binding=bindings(bindings_using_impulses_for_collisions(i));particles.V(binding.x)=particles.V(binding.y);}}
    else for(int i=0;i<bindings.m;i++){particles.V(bindings(i).x)=particles.V(bindings(i).y);}
}
//#####################################################################
// Function Update_Binding_Index_From_Particle_Index
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Update_Binding_Index_From_Particle_Index()
{
    binding_index_from_particle_index.Clean_Memory();
    int max_particle_index=-1;for(int b=0;b<bindings.m;b++) max_particle_index=max(max_particle_index,bindings(b).x);
    binding_index_from_particle_index.Resize(max_particle_index+1);
    for(int b=0;b<bindings.m;b++){assert(!binding_index_from_particle_index(bindings(b).x));binding_index_from_particle_index(bindings(b).x)=b;}
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Read(TYPED_ISTREAM& input)
{
    SEGMENT_MESH* binding_mesh_save=binding_mesh;binding_mesh=0; // save binding mesh so we don't kill it in Clean_Memory
    Clean_Memory();int backward_compatible;Read_Binary(input,backward_compatible,bindings);
    if(bindings.m) Read_Binary(input,use_impulses_for_collisions); // TODO: remove this backwards compatibility hack after Siggraph
    Update_Binding_Index_From_Particle_Index();
    if(binding_mesh_save){binding_mesh=binding_mesh_save;binding_mesh->Initialize_Mesh(particles.Size(),bindings);}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void SOFT_BINDINGS<TV>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,2,bindings,use_impulses_for_collisions);
}
//#####################################################################
namespace PhysBAM{
template class SOFT_BINDINGS<VECTOR<float,1> >;
template class SOFT_BINDINGS<VECTOR<float,2> >;
template class SOFT_BINDINGS<VECTOR<float,3> >;
template class SOFT_BINDINGS<VECTOR<double,1> >;
template class SOFT_BINDINGS<VECTOR<double,2> >;
template class SOFT_BINDINGS<VECTOR<double,3> >;
}
