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

template<class TV> class CONTACT;
template<class TV> void Compute_Tangent_Helper(CONTACT<TV>& contact);

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
    
    CONTACT()
        :id(0,0)
    {}

    CONTACT(RIGID_BODY<TV>& body_1,RIGID_BODY<TV>& body_2,TV& _location,TV& _normal,T _distance,T dt)
        :location(_location),normal(_normal),distance(_distance),coefficient_of_friction(RIGID_BODY<TV>::Coefficient_Of_Friction(body_1,body_2))
    {
        id(0)=body_1.particle_index;
        id(1)=body_2.particle_index;

        TV r_1=location-body_1.Frame().t;
        TV r_2=location-body_2.Frame().t;

        normal_constraint(0).linear=-normal;
        normal_constraint(0).angular=TV::Cross_Product(normal,r_1);

        normal_constraint(1).linear=normal;
        normal_constraint(1).angular=-TV::Cross_Product(normal,r_2);

        inverse_mass_times_normal_constraint(0)=body_1.Inertia_Inverse_Times(normal_constraint(0));
        inverse_mass_times_normal_constraint(1)=body_2.Inertia_Inverse_Times(normal_constraint(1));

        normal_diagonal=0;
        if(!body_1.Has_Infinite_Inertia())
            normal_diagonal+=TWIST<TV>::Dot_Product(normal_constraint(0),inverse_mass_times_normal_constraint(0));
        if(!body_2.Has_Infinite_Inertia())
            normal_diagonal+=TWIST<TV>::Dot_Product(normal_constraint(1),inverse_mass_times_normal_constraint(1));
        normal_relative_velocity=-distance/dt;

        Compute_Tangent_Helper(*this);
        
        for(int i=1;i<d;i++)
        {
            tangent_constraint(0)(i).linear=-tangent(i);
            tangent_constraint(0)(i).angular=TV::Cross_Product(tangent(i),r_1);

            tangent_constraint(1)(i).linear=tangent(i);
            tangent_constraint(1)(i).angular=-TV::Cross_Product(tangent(i),r_2);

            inverse_mass_times_tangent_constraint(0)(i)=body_1.Inertia_Inverse_Times(tangent_constraint(0)(i));
            inverse_mass_times_tangent_constraint(1)(i)=body_2.Inertia_Inverse_Times(tangent_constraint(1)(i));

            tangent_diagonal(i)=0;
            if(!body_1.Has_Infinite_Inertia())
                tangent_diagonal(i)+=TWIST<TV>::Dot_Product(tangent_constraint(0)(i),inverse_mass_times_tangent_constraint(0)(i));
            if(!body_2.Has_Infinite_Inertia())
                tangent_diagonal(i)+=TWIST<TV>::Dot_Product(tangent_constraint(1)(i),inverse_mass_times_tangent_constraint(1)(i));
        }
    }

    void Compute_Normal_Constraints();
    void Compute_Tangential_Constraints();
};

template<class T> void Compute_Tangent_Helper(CONTACT<VECTOR<T,1> >& contact) {}
template<class T> void Compute_Tangent_Helper(CONTACT<VECTOR<T,2> >& contact)
{
    contact.tangent(0)=contact.normal.Unit_Orthogonal_Vector();
}
template<class T> void Compute_Tangent_Helper(CONTACT<VECTOR<T,3> >& contact)
{
    contact.tangent(0)=contact.normal.Unit_Orthogonal_Vector();
    contact.tangent(1)=VECTOR<T,3>::Cross_Product(contact.normal,contact.tangent(0)).Normalized();
}

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
    RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const bool correct_contact_energy,const int max_iterations,
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
