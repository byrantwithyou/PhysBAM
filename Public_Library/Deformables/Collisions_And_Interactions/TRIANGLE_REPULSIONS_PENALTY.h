//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS_PENALTY
//#####################################################################
#ifndef __TRIANGLE_REPULSIONS_PENALTY__
#define __TRIANGLE_REPULSIONS_PENALTY__

#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Forces/LAGGED_FORCE.h>
namespace PhysBAM{

template<class TV>
class TRIANGLE_REPULSIONS_PENALTY:public LAGGED_FORCE<TV>
{
private:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    using LAGGED_FORCE<TV>::particles;using typename LAGGED_FORCE<TV>::FREQUENCY_DATA;
    ARRAY<REPULSION_PAIR<TV> >& interaction_pairs;

    T stiffness;

    ARRAY<int> bad_pairs;
    ARRAY<T> volume;
    T pe;
    ARRAY<VECTOR<TV,TV::m+1> > grad_pe;
    ARRAY<VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1> > H_pe;

    TRIANGLE_REPULSIONS_PENALTY(DEFORMABLE_PARTICLES<TV>& particles_input,ARRAY<REPULSION_PAIR<TV> >& interaction_pairs,T stiffness=1e4)
        :LAGGED_FORCE<TV>(particles_input),interaction_pairs(interaction_pairs),stiffness(stiffness)
    {}

    virtual ~TRIANGLE_REPULSIONS_PENALTY()
    {}

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Lagged_Update_Position_Based_State(const T time) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    T Potential_Energy(const T time) const override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    T CFL_Strain_Rate() const override;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) override;
    void Penalty(T original_volume,const VECTOR<int,4>& nodes, const ARRAY_VIEW<TV, int>&X, T& e, VECTOR<TV,4>& de, VECTOR<VECTOR<MATRIX<T,TV::m>,4>,4>& he);
//#####################################################################
};
}
#endif
