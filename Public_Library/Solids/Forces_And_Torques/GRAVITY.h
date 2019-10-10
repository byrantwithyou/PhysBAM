//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRAVITY
//#####################################################################
#ifndef __GRAVITY__
#define __GRAVITY__

#include <Solids/Forces_And_Torques/POINTWISE_FORCE.h>
namespace PhysBAM{

template<class TV>
class GRAVITY:public POINTWISE_FORCE<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef POINTWISE_FORCE<TV> BASE;
    using BASE::particles;using BASE::rigid_body_collection;using BASE::influenced_particles;
    using BASE::mpi_solids;using BASE::force_particles;using BASE::force_rigid_body_particles;

    TV gravity;
public:

    GRAVITY(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_particles_input,
        ARRAY<int>* influenced_rigid_body_particles_input,const TV& gravity_input=-(T)9.8*TV::Axis_Vector(1-(TV::m==1)))
        :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,influenced_particles_input,influenced_rigid_body_particles_input),gravity(gravity_input)
    {}

    GRAVITY(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_particles_input,
        const bool influence_all_rigid_body_particles_input,const TV& gravity_input=-(T)9.8*TV::Axis_Vector(1-(TV::m==1)))
        :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,influence_all_particles_input,influence_all_rigid_body_particles_input),
        gravity(gravity_input)
    {}

    template<class T_MESH>
    GRAVITY(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const T_MESH& mesh,ARRAY<int>* influenced_rigid_body_particles_input,
        const TV& gravity_input=-(T)9.8*TV::Axis_Vector(1-(TV::m==1)))
        :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,mesh,influenced_rigid_body_particles_input),gravity(gravity_input)
    {
        Get_Unique(*influenced_particles,mesh.elements.Flattened());
    }

    virtual ~GRAVITY()
    {}

    void Add_Velocity_Dependent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time) const override
    {}

    void Enforce_Definiteness(const bool enforce_definiteness_input) override
    {}

    void Add_Implicit_Velocity_Independent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time,bool transpose=false) const override
    {}

    int Velocity_Dependent_Forces_Size() const override
    {return 0;}

    void Add_Velocity_Dependent_Forces_First_Half(const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> aggregate,const T time) const override
    {}

    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,GENERALIZED_VELOCITY<TV>& F,const T time) const override
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const override;
    T Potential_Energy(const T time) const override;
//#####################################################################
};
}
#endif
