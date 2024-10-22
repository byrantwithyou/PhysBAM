//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION
//#####################################################################
#ifndef __RIGID_DEFORMABLE_PENALTY_WITH_FRICTION__
#define __RIGID_DEFORMABLE_PENALTY_WITH_FRICTION__

#include <Solids/Forces_And_Torques/SOLIDS_FORCES.h>
#include <functional>
namespace PhysBAM{

template<class TV> class MOVE_RIGID_BODY_DIFF;
template<class TV>
class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION:public SOLIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef SOLIDS_FORCES<TV> BASE;
    typedef int HAS_TYPED_READ_WRITE;
    using BASE::particles;using BASE::rigid_body_collection;

    T stiffness_coefficient=0;
    T friction=0;
    int num_dynamic=0,num_stick=0;
    bool use_bisection=false;
    
    struct COLLISION_PAIR
    {
        int p; // colliding particle
        int b; // colliding rigid body
        TV X; // original attachment point (object space)
        TV Y; // relaxed attachment point (world space)
        int e,e0;
        TV w,w0; // barycentric coords of X, for CCD
        MATRIX<T,TV::m> dYdZ; // Dependence of Y on X(p)
        MATRIX<T,TV::m> dYdL; // Dependence of Y on twist.linear
        MATRIX<T,TV::m,TV::SPIN::m> dYdA; // Dependence of Y on twist.angular
        bool active;

        template<class RW> void Write(std::ostream& output) const
        {Write_Binary<RW>(output,p,b,X);}

        template<class RW> void Read(std::istream& input)
        {Read_Binary<RW>(input,p,b,X);}
    };
    ARRAY<COLLISION_PAIR> collision_pairs;
    HASHTABLE<PAIR<int,int> > hash;
    std::function<void()> get_candidates=0; // Call Add_Pair on collision candidates.
    const ARRAY<MOVE_RIGID_BODY_DIFF<TV> >& move_rb_diff;

    RIGID_DEFORMABLE_PENALTY_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        const ARRAY<MOVE_RIGID_BODY_DIFF<TV> >& move_rb_diff);
    virtual ~RIGID_DEFORMABLE_PENALTY_WITH_FRICTION();

//#####################################################################
    void Update_Position_Based_State(const T time) override;
    void Add_Implicit_Velocity_Independent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time,bool transpose=false) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    void Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const override;
    void Add_Velocity_Dependent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,GENERALIZED_VELOCITY<TV>& F,const T time) const override;
    T CFL_Strain_Rate() const override;
    void Initialize_CFL(ARRAY_VIEW<typename BASE::DEFORMABLE_FREQUENCY_DATA> frequency,ARRAY_VIEW<typename BASE::RIGID_FREQUENCY_DATA> rigid_frequency) override;
    void Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input=true) override;
    void Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input=true,const T max_strain_per_time_step_input=.1) override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    T Potential_Energy(const T time) const override;
    void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const override;
    void Relax_Attachment(int cp);
    void Relax_Attachment_Implicit(int cp);
    void Relax_Attachment_Mesh(int cp);
    void Update_Attachments_And_Prune_Pairs();
    void Add_Pair(int p,int b);
    void Add_Pair(int p,int b,int e,const TV& X0,const FRAME<TV>& f,T thickness); // CCD
    void Read(TYPED_ISTREAM input);
    void Write(TYPED_OSTREAM output) const;
//#####################################################################
};
}
#endif
