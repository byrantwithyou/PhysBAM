//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SOLVE_CONTACT
//##################################################################### 
#ifndef __SOLVE_CONTACT__
#define __SOLVE_CONTACT__

#include <Tools/Data_Structures/HASHTABLE.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISIONS.h>

namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGID_BODY_COLLISION_PARAMETERS;
template<class TV> class RIGID_BODY_COLLISIONS;
template<class TV> class RIGIDS_COLLISION_CALLBACKS;
template<class TV> class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class RIGID_BODY_CONTACT_GRAPH;
template<class TV> class MULTIBODY_LEVELSET_IMPLICIT_OBJECT;

namespace SOLVE_CONTACT
{
template<class TV> struct ANALYTICS
{
    typedef typename TV::SCALAR T;
    typedef bool (*UPDATE_ANALYTIC_CONTACT_PAIR_T)(RIGID_BODY_COLLISIONS<TV>&,RIGIDS_COLLISION_CALLBACKS<TV>&,const int,const int,IMPLICIT_OBJECT<TV>*,IMPLICIT_OBJECT<TV>*,const bool,const int,const T,const T,const T,const bool);
};

template<class TV> class CONTACT
{
public:
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::dimension,s=T_SPIN::dimension};

    VECTOR<int,2> id;

    TV location;
    TV normal;
    VECTOR<TV,d-1> tangent;
    T distance;
    T coefficient_of_friction;

    VECTOR<TWIST<TV>,2> normal_constraint;
    VECTOR<TWIST<TV>,2> inverse_mass_times_normal_constraint;
    T normal_relative_velocity;

    VECTOR<VECTOR<TWIST<TV>,d-1>,2> tangent_constraint;
    VECTOR<VECTOR<TWIST<TV>,d-1>,2> inverse_mass_times_tangent_constraint;
    
    T normal_diagonal;
    VECTOR<T,d-1> tangent_diagonal;
    
    CONTACT();
    CONTACT(RIGID_BODY<TV>& body_1,RIGID_BODY<TV>& body_2,TV& _location,TV& _normal,T _distance,T dt);
    void Compute_Normal_Constraints();
    void Compute_Tangential_Constraints();
};

template<class TV>
void Push_Out(RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact);
template<class TV> void Solve(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,
    RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters,const bool correct_contact_energy,const bool use_saved_pairs,const typename TV::SCALAR dt,
    const typename TV::SCALAR time,ARRAY<TWIST<TV> >& mpi_rigid_velocity_save,ARRAY<typename TV::SPIN>& mpi_rigid_angular_momentum_save);
template<class TV> bool Update_Analytic_Multibody_Contact(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,
    const int id_1,const int id_2,MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>& multibody,IMPLICIT_OBJECT<TV>& levelset,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost);
template<class TV> bool Update_Analytic_Multibody_Contact(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,
    RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost);
template<class TV> bool Update_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,HASHTABLE<VECTOR<std::string,2>,
    typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,const int id_1,const int id_2,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost);
template<class TV> void Update_Contact_Pair_Helper(RIGID_BODY_COLLISIONS<TV>& rigid_body_collision,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id_1,const int id_2,
    const typename TV::SCALAR dt,const typename TV::SCALAR time,const typename TV::SCALAR epsilon_scale,const TV& collision_location,const TV& collision_normal,
    const TV& collision_relative_velocity,const bool correct_contact_energy,const bool rolling_friction,const bool mpi_one_ghost);
template<class TV> void Euler_Step_Position(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id,
    const typename TV::SCALAR dt,const typename TV::SCALAR time);
template<class TV> void Register_Analytic_Contacts(HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry);
template<class TV> bool Solve_Projected_Gauss_Seidel(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact,typename TV::SCALAR desired_separation_distance,typename TV::SCALAR contact_proximity,
    typename TV::SCALAR dt,typename TV::SCALAR tolerance);
template<class TV> void Get_Contact_Points(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters,ARRAY<VECTOR<int,2> >& pairs,
    ARRAY<CONTACT<TV> >& contacts,typename TV::SCALAR contact_proximity,typename TV::SCALAR dt,const bool stagger_points,const bool use_old_states);

}
}
#endif
