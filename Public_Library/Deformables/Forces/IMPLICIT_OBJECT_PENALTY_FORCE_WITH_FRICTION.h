//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION
//#####################################################################
#ifndef __IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION__
#define __IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <functional>
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;

    ARRAY<IMPLICIT_OBJECT<TV>*> ios;
    T stiffness_coefficient;
    T friction;
    
    struct COLLISION_PAIR
    {
        int p; // colliding particle
        int o; // colliding implicit object
        TV X; // original attachment point
        TV Y; // relaxed attachment point
        MATRIX<T,TV::m> dYdZ; // Dependence of Y on X(p)
        bool active;
    };
    ARRAY<COLLISION_PAIR> collision_pairs;
    HASHTABLE<PAIR<int,int> > hash;
    std::function<void()> get_candidates=0; // Call Add_Pair on collision candidates.

    IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input,
        T stiffness_coefficient,T friction);
    virtual ~IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION();

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const override;
    void Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency) override;
    T CFL_Strain_Rate() const override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    T Potential_Energy(const T time) const override;
    void Relax_Attachment(int cp);
    void Update_Attachments_And_Prune_Pairs();
    void Add_Pair(int p,int b);
//#####################################################################
};

template<class TV>
struct RELAX_ATTACHMENT_HELPER
{
    TV Y;
    MATRIX<typename TV::SCALAR,TV::m> dYdZ,dYdX,dYdW;
};
template<class TV,class T> RELAX_ATTACHMENT_HELPER<TV>
Relax_Attachment_Helper(const TV& Z,const TV& X,const TV& W,T mu);
}
#endif
