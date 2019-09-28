//#####################################################################
// Copyright 2007-2008, Michael Lentine, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/IDENTITY_ARRAY.h>
#include <Geometry/Topology/TETRAHEDRON_MESH.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/POINTWISE_FORCE.h>
using namespace PhysBAM;
template<class TV> POINTWISE_FORCE<TV>::
POINTWISE_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_particles_input,
    ARRAY<int>* influenced_rigid_body_particles_input)
    :SOLIDS_FORCES<TV>(particles_input,rigid_body_collection_input),influenced_particles(influenced_particles_input),
    influenced_rigid_body_particles(influenced_rigid_body_particles_input),
    need_destroy_influenced_particles(false),need_destroy_influenced_rigid_body_particles(false),influence_all_particles(false),influence_all_rigid_body_particles(false),mpi_solids(0)
{
}
template<class TV> POINTWISE_FORCE<TV>::
POINTWISE_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_particles_input,
    const bool influence_all_rigid_body_particles_input)
    :SOLIDS_FORCES<TV>(particles_input,rigid_body_collection_input),influenced_particles(0),influenced_rigid_body_particles(0),need_destroy_influenced_particles(true),
    need_destroy_influenced_rigid_body_particles(true),influence_all_particles(influence_all_particles_input),influence_all_rigid_body_particles(influence_all_rigid_body_particles_input),
    mpi_solids(0)
{
}
template<class TV> template<class T_MESH> POINTWISE_FORCE<TV>::
POINTWISE_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const T_MESH& mesh,ARRAY<int>* influenced_rigid_body_particles_input)
    :SOLIDS_FORCES<TV>(particles_input,rigid_body_collection_input),influenced_particles(new ARRAY<int>),
    influenced_rigid_body_particles(influenced_rigid_body_particles_input),
    need_destroy_influenced_particles(true),need_destroy_influenced_rigid_body_particles(false),influence_all_particles(false),influence_all_rigid_body_particles(false)
{
    Get_Unique(*influenced_particles,mesh.elements.Flattened());
}
template<class TV> POINTWISE_FORCE<TV>::
~POINTWISE_FORCE()
{
    if(need_destroy_influenced_particles) delete influenced_particles;
    if(need_destroy_influenced_rigid_body_particles) delete influenced_rigid_body_particles;
}
//#####################################################################
// Function Get_Rigid_Body_Particle_List
//#####################################################################
template<class TV> template<class T_ARRAY> ARRAY<int> POINTWISE_FORCE<TV>::
Get_Rigid_Body_Particle_List(const T_ARRAY& array)
{
    LOG::cout<<"attempting to add "<<array.Size()<<" rigid body particles"<<std::endl;
    ARRAY<int> list;
    for(int i=0;i<array.Size();i++){
        int p=array(i);
        if(rigid_body_collection.Is_Active(p) && !rigid_body_collection.Rigid_Body(p).Has_Infinite_Inertia()) list.Append(p);}
    return list;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void POINTWISE_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    if(influence_all_particles)
        Update_Force_Particles(force_particles,
            IDENTITY_ARRAY<>(particle_is_simulated.m),
            particle_is_simulated,false);
    else if(influenced_particles)
        Update_Force_Particles(force_particles,*influenced_particles,
            particle_is_simulated,false);

    auto g=[&](const auto& a)
        {
            for(int i:a)
                if(rigid_particle_is_simulated(i) &&
                    rigid_body_collection.Is_Active(i) &&
                    !rigid_body_collection.Rigid_Body(i).Has_Infinite_Inertia())
                    force_rigid_body_particles.Append(i);
        };

    if(influence_all_rigid_body_particles)
        g(IDENTITY_ARRAY<>(rigid_particle_is_simulated.m));
    else if(influenced_rigid_body_particles)
        g(*influenced_rigid_body_particles);
}
//#####################################################################
namespace PhysBAM{
template class POINTWISE_FORCE<VECTOR<float,1> >;
template class POINTWISE_FORCE<VECTOR<float,2> >;
template class POINTWISE_FORCE<VECTOR<float,3> >;
template POINTWISE_FORCE<VECTOR<float,3> >::POINTWISE_FORCE(DEFORMABLE_PARTICLES<VECTOR<float,3> >&,RIGID_BODY_COLLECTION<VECTOR<float,3> >&,TETRAHEDRON_MESH const&,ARRAY<int,int>*);
template POINTWISE_FORCE<VECTOR<float,3> >::POINTWISE_FORCE(DEFORMABLE_PARTICLES<VECTOR<float,3> >&,RIGID_BODY_COLLECTION<VECTOR<float,3> >&,TRIANGLE_MESH const&,ARRAY<int,int>*);
template class POINTWISE_FORCE<VECTOR<double,1> >;
template class POINTWISE_FORCE<VECTOR<double,2> >;
template class POINTWISE_FORCE<VECTOR<double,3> >;
template POINTWISE_FORCE<VECTOR<double,3> >::POINTWISE_FORCE(DEFORMABLE_PARTICLES<VECTOR<double,3> >&,RIGID_BODY_COLLECTION<VECTOR<double,3> >&,TETRAHEDRON_MESH const&,ARRAY<int,int>*);
template POINTWISE_FORCE<VECTOR<double,3> >::POINTWISE_FORCE(DEFORMABLE_PARTICLES<VECTOR<double,3> >&,RIGID_BODY_COLLECTION<VECTOR<double,3> >&,TRIANGLE_MESH const&,ARRAY<int,int>*);
}
