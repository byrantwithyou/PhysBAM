//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_PENALTY_WITH_FRICTION
//#####################################################################
#ifndef __RIGID_PENALTY_WITH_FRICTION__
#define __RIGID_PENALTY_WITH_FRICTION__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Matrices/MATRIX.h>
#include <Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
#include <functional>
namespace PhysBAM{

template<class TV> class MOVE_RIGID_BODY_DIFF;
template<class TV> class IMPLICIT_OBJECT;
template<class TV>
class RIGID_PENALTY_WITH_FRICTION:public RIGIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef RIGIDS_FORCES<TV> BASE;
    typedef int HAS_TYPED_READ_WRITE;
    using BASE::rigid_body_collection;

    T stiffness_coefficient=0;
    T friction=0;
    bool use_bisection=false;
    
    struct COLLISION_PAIR
    {
        int bs,v; // colliding vertex (simplicial mesh)
        int bi; // colliding rigid body (implicit object)
        TV X; // original attachment point (object space)
        TV Y; // relaxed attachment point (world space)
        TV Z; // simplicial point
        MATRIX<T,TV::m> dYdLs,dYdLi,dZdLs; // Dependence on twist.linear
        MATRIX<T,TV::m,TV::SPIN::m> dYdAs,dYdAi,dZdAs; // Dependence on twist.angular
        bool active;

        template<class RW> void Write(std::ostream& output) const
        {Write_Binary<RW>(output,bs,v,bi,X);}

        template<class RW> void Read(std::istream& input)
        {Read_Binary<RW>(input,bs,v,bi,X);}
    };
    ARRAY<COLLISION_PAIR> collision_pairs;
    HASHTABLE<TRIPLE<int,int,int> > hash;
    std::function<void()> get_candidates=0; // Call Add_Pair on collision candidates.
    const ARRAY<MOVE_RIGID_BODY_DIFF<TV> >& move_rb_diff;

    RIGID_PENALTY_WITH_FRICTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        const ARRAY<MOVE_RIGID_BODY_DIFF<TV> >& move_rb_diff);
    virtual ~RIGID_PENALTY_WITH_FRICTION();

//#####################################################################
    void Update_Position_Based_State(const T time) override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time,bool transpose=false) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const override;
    T CFL_Strain_Rate() const override;
    void Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> rigid_frequency) override;
    void Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input=true) override;
    void Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input=true,const T max_strain_per_time_step_input=.1) override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated) override;
    T Potential_Energy(const T time) const override;
    void Relax_Attachment(int cp);
    void Update_Attachments_And_Prune_Pairs();
    void Add_Pair(int b1,int v,int b2);
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
//#####################################################################
};
template<class TV>
struct MOVING_LEVEL_SET_HELPER
{
    typedef typename TV::SCALAR T;
    MATRIX<T,TV::m> dUdX,dUdL,dndU;
    MATRIX<T,TV::m,TV::SPIN::m> dUdA,dndA;
    TV U,N,n;
    T phi;

    bool Init(const MOVE_RIGID_BODY_DIFF<TV>& mr,const IMPLICIT_OBJECT<TV>* io,const TV& X,bool exit_if_sep);
};
template<class TV,class T> void
Project_Attachment_To_Surface(TV& W,MOVING_LEVEL_SET_HELPER<TV>& s,
    const TV& X,MATRIX<T,TV::m>& dWdX,
    MATRIX<T,TV::m>& dWdL,MATRIX<T,TV::m,TV::SPIN::m>& dWdA);
}
#endif
