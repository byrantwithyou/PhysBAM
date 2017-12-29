//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION
//#####################################################################
#ifndef __RIGID_DEFORMABLE_PENALTY_WITH_FRICTION__
#define __RIGID_DEFORMABLE_PENALTY_WITH_FRICTION__

#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Solids/Forces_And_Torques/SOLIDS_FORCES.h>
#include <functional>
namespace PhysBAM{

template<class TV>
class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION:public SOLIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef SOLIDS_FORCES<TV> BASE;
    using BASE::particles;using BASE::rigid_body_collection;

    T stiffness_coefficient;
    T friction;
    
    struct COLLISION_PAIR
    {
        int p; // colliding particle
        int b; // colliding rigid body
        TV X; // original attachment point (object space)
        TV Y; // relaxed attachment point (world space)
        bool active;
    };
    ARRAY<COLLISION_PAIR> collision_pairs;
    HASHTABLE<PAIR<int,int> > hash;
    std::function<void()> get_candidates=0; // Call Add_Pair on collision candidates.

    RIGID_DEFORMABLE_PENALTY_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,T stiffness_coefficient,
        T friction);
    virtual ~RIGID_DEFORMABLE_PENALTY_WITH_FRICTION();

//#####################################################################
    void Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input=true) override;
    void Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input=true,const T max_strain_per_time_step_input=.1) override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    void Update_Position_Based_State(const T time) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    T CFL_Strain_Rate() const override;
    void Initialize_CFL(ARRAY_VIEW<typename BASE::DEFORMABLE_FREQUENCY_DATA> frequency,ARRAY_VIEW<typename BASE::RIGID_FREQUENCY_DATA> rigid_frequency) override;
    T Potential_Energy(const T time) const override;
    void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const override;
    void Relax_Attachment(int cp);
    void Update_Attachments_And_Prune_Pairs();
    void Add_Pair(int p,int b);
//#####################################################################
};
}
#endif
