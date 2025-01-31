//#####################################################################
// Copyright 2007-2009, Michael Lentine, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINTWISE_FORCE
//#####################################################################
#ifndef __POINTWISE_FORCE__
#define __POINTWISE_FORCE__

#include <Core/Vectors/VECTOR.h>
#include <Tools/Particles/PARTICLES_FORWARD.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <Solids/Forces_And_Torques/SOLIDS_FORCES.h>
namespace PhysBAM{

template<class TV>
class POINTWISE_FORCE:public SOLIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef SOLIDS_FORCES<TV> BASE;
    using BASE::particles;using BASE::rigid_body_collection;
    typedef typename BASE::RIGID_FREQUENCY_DATA RIGID_FREQUENCY_DATA;
    typedef typename BASE::DEFORMABLE_FREQUENCY_DATA DEFORMABLE_FREQUENCY_DATA;
protected:
    ARRAY<int> *influenced_particles,*influenced_rigid_body_particles;
    bool need_destroy_influenced_particles,need_destroy_influenced_rigid_body_particles;
    bool influence_all_particles,influence_all_rigid_body_particles;
    MPI_SOLIDS<TV>* mpi_solids;
public:
    ARRAY<int> force_particles;
    ARRAY<int> force_rigid_body_particles;

    POINTWISE_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_particles_input,
        ARRAY<int>* influenced_rigid_body_particles_input);
    POINTWISE_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_particles_input,
        const bool influence_all_rigid_body_particles_input);
    template<class T_MESH>
    POINTWISE_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const T_MESH& mesh,ARRAY<int>* influenced_rigid_body_particles_input);
    virtual ~POINTWISE_FORCE();

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override
    {}

    T CFL_Strain_Rate() const override
    {return FLT_MAX;}

    void Initialize_CFL(ARRAY_VIEW<DEFORMABLE_FREQUENCY_DATA> frequency,ARRAY_VIEW<RIGID_FREQUENCY_DATA> rigid_frequency) override
    {}

protected:
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;

    template<class T_ARRAY>
    ARRAY<int> Get_Rigid_Body_Particle_List(const T_ARRAY& array);
//#####################################################################
};
}
#endif
