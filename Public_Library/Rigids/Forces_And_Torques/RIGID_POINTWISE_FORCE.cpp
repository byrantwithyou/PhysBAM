//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/IDENTITY_ARRAY.h>
#include <Rigids/Forces_And_Torques/RIGID_POINTWISE_FORCE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Get_Rigid_Body_Particle_List
//#####################################################################
template<class TV> template<class T_ARRAY> ARRAY<int> RIGID_POINTWISE_FORCE<TV>::
Get_Rigid_Body_Particle_List(const T_ARRAY& array)
{
    LOG::cout<<"attempting to add "<<array.Size()<<" rigid body particles"<<std::endl;
    ARRAY<int> list;
    for(int i=0;i<array.Size();i++){
        int p=array(i);if(rigid_body_collection.Is_Active(p) && !rigid_body_collection.Rigid_Body(p).Has_Infinite_Inertia()) list.Append(p);}
    return list;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void RIGID_POINTWISE_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated)
{
    auto g=[&](const auto& a)
        {
            for(int i:a)
                if(particle_is_simulated(i) &&
                    rigid_body_collection.Is_Active(i) &&
                    !rigid_body_collection.Rigid_Body(i).Has_Infinite_Inertia())
                    force_rigid_body_particles.Append(i);
        };

    if(influenced_rigid_body_particles)
        g(*influenced_rigid_body_particles);
    else g(IDENTITY_ARRAY<>(particle_is_simulated.m));
}
//#####################################################################
namespace PhysBAM{
template class RIGID_POINTWISE_FORCE<VECTOR<float,1> >;
template class RIGID_POINTWISE_FORCE<VECTOR<float,2> >;
template class RIGID_POINTWISE_FORCE<VECTOR<float,3> >;
template class RIGID_POINTWISE_FORCE<VECTOR<double,1> >;
template class RIGID_POINTWISE_FORCE<VECTOR<double,2> >;
template class RIGID_POINTWISE_FORCE<VECTOR<double,3> >;
}
