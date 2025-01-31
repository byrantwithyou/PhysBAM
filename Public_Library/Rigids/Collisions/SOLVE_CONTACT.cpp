//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/SCOPE.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Collisions/BOX_BOX_CONTACT_PAIR.h>
#include <Rigids/Collisions/BOX_PLANE_CONTACT_PAIR.h>
#include <Rigids/Collisions/CONTACT_PRECONDITIONER.h>
#include <Rigids/Collisions/LEVELSET_CONTACT_PAIR.h>
#include <Rigids/Collisions/PARTICLES_IN_IMPLICIT_OBJECT.h>
#include <Rigids/Collisions/PARTICLES_IN_PROXIMITY.h>
#include <Rigids/Collisions/PROJECTED_GAUSS_SEIDEL.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <Rigids/Collisions/RIGID_BODY_CONTACT_GRAPH.h>
#include <Rigids/Collisions/RIGID_BODY_SKIP_COLLISION_CHECK.h>
#include <Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <Rigids/Collisions/SOLVE_CONTACT.h>
#include <Rigids/Collisions/SPHERE_PLANE_CONTACT_PAIR.h>
#include <Rigids/Collisions/SPHERE_SPHERE_CONTACT_PAIR.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>

namespace PhysBAM
{
namespace SOLVE_CONTACT
{
template<class T> static void Compute_Tangent_Helper(CONTACT<VECTOR<T,1> >& contact) {}
template<class T> static void Compute_Tangent_Helper(CONTACT<VECTOR<T,2> >& contact)
{
    contact.tangent(0)=contact.normal.Unit_Orthogonal_Vector();
}
template<class T> static void Compute_Tangent_Helper(CONTACT<VECTOR<T,3> >& contact)
{
    contact.tangent(0)=contact.normal.Unit_Orthogonal_Vector();
    contact.tangent(1)=VECTOR<T,3>::Cross_Product(contact.normal,contact.tangent(0)).Normalized();
}

//#####################################################################
// Constructor
//#####################################################################
template<class TV> CONTACT<TV>::
CONTACT()
    :id(0,0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> CONTACT<TV>::
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
        
    for(int i=0;i<d-1;i++)
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
template<class TV>
bool Solve_Projected_Gauss_Seidel(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact,typename TV::SCALAR desired_separation_distance,typename TV::SCALAR contact_proximity,typename TV::SCALAR dt,typename TV::SCALAR tolerance,int iteration_maximum);
//#####################################################################
// Function Solve
//#####################################################################
template<class TV>
void Solve(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,
    RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters,const bool correct_contact_energy,const bool use_saved_pairs,const typename TV::SCALAR dt,const typename TV::SCALAR time,
    ARRAY<TWIST<TV> >& mpi_rigid_velocity_save,ARRAY<typename TV::SPIN>& mpi_rigid_angular_momentum_save)
{
    typedef typename TV::SCALAR T;

    LOG::SCOPE scope("rigid body contact kernel");

    ARRAY<ARRAY<VECTOR<int,2> > >& contact_pairs_for_level=use_saved_pairs?rigid_body_collisions.saved_contact_pairs_for_level:rigid_body_collisions.precomputed_contact_pairs_for_level;
    if(!parameters.use_projected_gauss_seidel){
        LOG::SCOPE scope("guendelman-bridson-fedkiw contact");

        HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T> analytic_contact_registry;
        if(parameters.use_analytic_collisions){Register_Analytic_Contacts<TV>(analytic_contact_registry);}

        ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body=&rigid_body_collisions.rigid_body_collection.articulated_rigid_body;
        rigid_body_collisions.skip_collision_check.Reset();bool need_another_iteration=true;int iteration=0;T epsilon_scale=1;
        while(need_another_iteration && iteration++<parameters.contact_iterations){
            need_another_iteration=false;
            if(parameters.use_epsilon_scaling) epsilon_scale=(T)iteration/parameters.contact_iterations;
            for(int level=0;level<rigid_body_collisions.contact_graph.Number_Of_Levels();level++){
                ARRAY<VECTOR<int,2> >& pairs=contact_pairs_for_level(level);
                bool need_another_level_iteration=true;int level_iteration=0;
                while(need_another_level_iteration && level_iteration++<rigid_body_collisions.contact_level_iterations){need_another_level_iteration=false;
                    if(parameters.use_epsilon_scaling_for_level) epsilon_scale=(T)iteration*level_iteration/(parameters.contact_iterations*rigid_body_collisions.contact_level_iterations);
                    for(int i=0;i<pairs.m;i++){int id_1=pairs(i)(0),id_2=pairs(i)(1);
                        if(rigid_body_collisions.skip_collision_check.Skip_Pair(id_1,id_2)) continue;
                        if(rigid_body_collisions.prune_contact_using_velocity){
                            TRIPLE<int,int,int>& pair_scale=rigid_body_collisions.pairs_scale.Get(pairs(i).Sorted());
                            if((rigid_body_collisions.contact_level_iterations-level_iteration)%pair_scale.y || (parameters.contact_iterations-iteration)%pair_scale.x) continue;
                            if(Update_Contact_Pair(rigid_body_collisions,collision_callbacks,analytic_contact_registry,id_1,id_2,correct_contact_energy,
                                    (int)(rigid_body_collisions.contact_pair_iterations/pair_scale.z),epsilon_scale,dt,time,false)){
                                if(!use_saved_pairs) rigid_body_collisions.saved_contact_pairs_for_level(level).Append_Unique(pairs(i)); // TODO: Potentially inefficient
                                need_another_level_iteration=true;need_another_iteration=true;}}
                        else{
                            if(Update_Contact_Pair(rigid_body_collisions,collision_callbacks,analytic_contact_registry,id_1,id_2,correct_contact_energy,
                                    rigid_body_collisions.contact_pair_iterations,epsilon_scale,dt,time,false)){
                                if(!use_saved_pairs) rigid_body_collisions.saved_contact_pairs_for_level(level).Append_Unique(pairs(i));
                                need_another_level_iteration=true;need_another_iteration=true;}}}
                    if(articulated_rigid_body)
                        for(int i=0;i<articulated_rigid_body->contact_level_iterations;i++) for(int j=0;j<articulated_rigid_body->process_list(level).m;j++){
                            rigid_body_collisions.Apply_Prestabilization_To_Joint(dt,time,*articulated_rigid_body,articulated_rigid_body->process_list(level)(j),epsilon_scale);
                            need_another_level_iteration=need_another_iteration=true;}}}}
        if(rigid_body_collisions.prune_stacks_from_contact) rigid_body_collisions.Apply_Stacking_Contact();}
    else{
        LOG::SCOPE scope("projected-gauss-seidel contact");

        if(!use_saved_pairs) rigid_body_collisions.saved_contact_pairs_for_level=contact_pairs_for_level;

        ARRAY<VECTOR<int,2> > pairs;
        for(int level=0;level<rigid_body_collisions.contact_graph.Number_Of_Levels();level++)
            pairs.Append_Unique_Elements(contact_pairs_for_level(level));
        
        int iteration_maximum=parameters.contact_iterations*rigid_body_collisions.contact_level_iterations*rigid_body_collisions.contact_pair_iterations;

        if(Solve_Projected_Gauss_Seidel(rigid_body_collection,collision_callbacks,pairs,rigid_body_collisions.pairs_processed_by_contact,rigid_body_collisions.desired_separation_distance,parameters.contact_proximity,dt,parameters.projected_gauss_seidel_tolerance,iteration_maximum))
        {
            /*LOG::cout << "----- end states -----" << std::endl;
            for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++)
            {
                LOG::cout << "twist " << i << " " << rigid_body_collection.Rigid_Body(i).Twist() << std::endl;
            }*/

            LOG::SCOPE scope_reevolve("reevolving bodies");
            for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++)
            {
                if(rigid_body_collection.Is_Active(i))
                {
                    RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
                    if(!body.Has_Infinite_Inertia())
                    {
                        body.Update_Angular_Momentum();
                        collision_callbacks.Restore_Position(i);
                        collision_callbacks.Euler_Step_Position(i,dt,time);
                        //body.Update_Angular_Momentum();
                    }
                }
            }
        }
    }
}
//#####################################################################
// Function Update_Analytic_Multibody_Contact
//#####################################################################
template<class TV> bool 
Update_Analytic_Multibody_Contact(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,
    const int id_1,const int id_2,MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>& multibody,IMPLICIT_OBJECT<TV>& levelset,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    bool return_val=false;
    IMPLICIT_OBJECT<TV>* levelset_base=&levelset;
    if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* levelset_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(levelset_base)) 
        levelset_base=levelset_transformed->object_space_implicit_object;
    for(int i=0;i<multibody.levelsets->m;i++){
        VECTOR<std::string,2> key(typeid(*dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >&>(*(*multibody.levelsets)(i)).object_space_implicit_object).name(),typeid(*levelset_base).name());
        if(typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T* contact_function=analytic_contact_registry.Get_Pointer(key.Sorted()))
            return_val|=(*contact_function)(rigid_body_collisions,collision_callbacks,id_1,id_2,(*multibody.levelsets)(i),&levelset,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
        else PHYSBAM_FATAL_ERROR();}
    return return_val;
}
//#####################################################################
// Function Update_Analytic_Multibody_Contact
//#####################################################################
template<class TV> bool 
Update_Analytic_Multibody_Contact(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,
    HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,
    RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    bool return_val=false;
    MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* multibody0=dynamic_cast<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>*>(body0.implicit_object->object_space_implicit_object);
    MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* multibody1=dynamic_cast<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>*>(body1.implicit_object->object_space_implicit_object);
    if(multibody0 && multibody1) for(int i=0;i<multibody1->levelsets->m;i++) 
        return_val|=Update_Analytic_Multibody_Contact(rigid_body_collisions,collision_callbacks,analytic_contact_registry,body0.particle_index,body1.particle_index,*multibody0,*(*multibody1->levelsets)(i),correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
    else if(multibody0) 
        return_val=Update_Analytic_Multibody_Contact(rigid_body_collisions,collision_callbacks,analytic_contact_registry,body0.particle_index,body1.particle_index,*multibody0,*body1.implicit_object->object_space_implicit_object,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
    else 
        return_val=Update_Analytic_Multibody_Contact(rigid_body_collisions,collision_callbacks,analytic_contact_registry,body0.particle_index,body1.particle_index,*multibody1,*body0.implicit_object->object_space_implicit_object,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
    return return_val;
}
//#####################################################################
// Function Update_Contact_Pair
//#####################################################################
template<class TV>
bool Update_Contact_Pair(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,HASHTABLE<VECTOR<std::string,2>,
    typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,const int id_1,const int id_2,const bool correct_contact_energy,const int max_iterations,
    const typename TV::SCALAR epsilon_scale,const typename TV::SCALAR dt,const typename TV::SCALAR time,const bool mpi_one_ghost)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;
    RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters=rigid_body_collisions.parameters;

    RIGID_BODY<TV>& body0=rigid_body_collection.Rigid_Body(id_1),&body1=rigid_body_collection.Rigid_Body(id_2);
    VECTOR<std::string,2> key(typeid(*body0.implicit_object->object_space_implicit_object).name(),typeid(*body1.implicit_object->object_space_implicit_object).name());
    if(key.x==typeid(MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>).name()||key.y==typeid(MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>).name()){
        if(!body0.simplicial_object && !body1.simplicial_object) 
            return Update_Analytic_Multibody_Contact(rigid_body_collisions,collision_callbacks,analytic_contact_registry,body0,body1,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);}
    if(typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T* contact_function=analytic_contact_registry.Get_Pointer(key.Sorted()))
        return (*contact_function)(rigid_body_collisions,collision_callbacks,id_1,id_2,body0.implicit_object->object_space_implicit_object,body1.implicit_object->object_space_implicit_object,correct_contact_energy,max_iterations,epsilon_scale,dt,time,mpi_one_ghost);
    return CONTACT_PAIRS::Update_Levelset_Contact_Pair(rigid_body_collisions,collision_callbacks,id_1,id_2,correct_contact_energy,max_iterations,epsilon_scale,dt,time,
        parameters.rigid_collisions_use_triangle_hierarchy,parameters.rigid_collisions_use_edge_intersection,parameters.rigid_collisions_use_triangle_hierarchy_center_phi_test,mpi_one_ghost);
}
//#####################################################################
// Function Update_Contact_Pair_Helper
//#####################################################################
template<class TV>
void Update_Contact_Pair_Helper(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id_1,const int id_2,
    const typename TV::SCALAR dt,const typename TV::SCALAR time,const typename TV::SCALAR epsilon_scale,const TV& collision_location,const TV& collision_normal,
    const TV& collision_relative_velocity,const bool correct_contact_energy,const bool rolling_friction,const bool mpi_one_ghost)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=rigid_body_collisions.rigid_body_collection;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings=rigid_body_collisions.rigid_body_cluster_bindings;

    RIGID_BODY<TV>& body0=rigid_body_collection.Rigid_Body(id_1);
    RIGID_BODY<TV>& body1=rigid_body_collection.Rigid_Body(id_2);
    int parent_id_1=rigid_body_cluster_bindings.Get_Parent(body0).particle_index;
    int parent_id_2=rigid_body_cluster_bindings.Get_Parent(body1).particle_index;
    RIGID_BODY<TV>& parent_body_1=rigid_body_collection.Rigid_Body(parent_id_1);
    RIGID_BODY<TV>& parent_body_2=rigid_body_collection.Rigid_Body(parent_id_2);
    rigid_body_collisions.rigid_body_particles_intersections.Set(Tuple(body0.particle_index,body1.particle_index,collision_location));
    RIGID_BODY<TV>::Apply_Collision_Impulse(parent_body_1,parent_body_2,body0.Frame().r,body1.Frame().r,collision_location,collision_normal,collision_relative_velocity,
        -1+epsilon_scale,RIGID_BODY<TV>::Coefficient_Of_Friction(parent_body_1,parent_body_2),false,rolling_friction,correct_contact_energy,mpi_one_ghost);
    collision_callbacks.Save_Position(parent_id_1);collision_callbacks.Save_Position(parent_id_2); // fix saved values & re-evolve bodies
    Euler_Step_Position(rigid_body_collisions,collision_callbacks,parent_id_1,dt,time);Euler_Step_Position(rigid_body_collisions,collision_callbacks,parent_id_2,dt,time);
    if(parent_id_2!=id_2 || parent_id_1!=id_1){
        rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Positions();
        rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();}
}
//#####################################################################
// Function Euler_Step_Position
//#####################################################################
template<class TV>
void Euler_Step_Position(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,const int id,const typename TV::SCALAR dt,const typename TV::SCALAR time)
{
    collision_callbacks.Euler_Step_Position(id,dt,time);
    rigid_body_collisions.rigid_body_collection.Rigid_Body(id).Update_Bounding_Box();rigid_body_collisions.skip_collision_check.Set_Last_Moved(id);
}
//#####################################################################
// Function Register_Analytic_Collisions
//#####################################################################
template<class TV>
void Register_Analytic_Contacts(HASHTABLE<VECTOR<std::string,2>,typename ANALYTICS<TV>::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry)
{
    const char* sphere=typeid(ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >).name();
    const char* plane=typeid(ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >).name();
    const char* box=typeid(ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >).name();
    analytic_contact_registry.Set(VECTOR<std::string,2>(sphere,sphere),CONTACT_PAIRS::Update_Sphere_Sphere_Contact_Pair<TV>);
    analytic_contact_registry.Set(VECTOR<std::string,2>(box,box),CONTACT_PAIRS::Update_Box_Box_Contact_Pair<TV>);
    analytic_contact_registry.Set(VECTOR<std::string,2>(sphere,plane).Sorted(),CONTACT_PAIRS::Update_Sphere_Plane_Contact_Pair<TV>);
    analytic_contact_registry.Set(VECTOR<std::string,2>(box,plane).Sorted(),CONTACT_PAIRS::Update_Box_Plane_Contact_Pair<TV>);
}
//#####################################################################
// Function Solve_Projected_Gauss_Seidel
//#####################################################################
template<class TV>
bool Solve_Projected_Gauss_Seidel(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact,typename TV::SCALAR desired_separation_distance,typename TV::SCALAR contact_proximity,typename TV::SCALAR dt,typename TV::SCALAR tolerance,int iteration_maximum)
{
    LOG::SCOPE scope("Solve_Projected_Gauss_Seidel");

    ARRAY<CONTACT<TV> > contacts;
    Get_Contact_Points(rigid_body_collection,collision_callbacks,pairs,contacts,contact_proximity,dt,true,true);
    
    pairs_processed_by_contact.Remove_All();
    for(int i=0;i<contacts.m;i++)
        pairs_processed_by_contact.Set(contacts(i).id.Sorted());
    
    if(contacts.m==0)
        return false;

    LOG::cout << "contacts " << contacts.m << std::endl;

    /*CONTACT_PRECONDITIONER::PRECONDITIONER<TV> preconditioner;
    preconditioner.Build_Fine_Problem(rigid_body_collection,contacts);
    preconditioner.Solve(tolerance,iteration_maximum);
    //preconditioner.Solve(tolerance,iteration_maximum);
    //preconditioner.Solve(tolerance,iteration_maximum);
    preconditioner.Update_Fine_Velocities(rigid_body_collection);

    exit(0);

    return true;*/

    return PROJECTED_GAUSS_SEIDEL::Solve(rigid_body_collection,contacts,tolerance,iteration_maximum);
}

//#####################################################################
// Function Get_Contact_Points
//#####################################################################
template<class TV>
void Get_Contact_Points(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,ARRAY<VECTOR<int,2> >& pairs,ARRAY<CONTACT<TV> >& contacts,typename TV::SCALAR contact_proximity,typename TV::SCALAR dt,const bool stagger_points,const bool use_old_states)
{
    LOG::SCOPE scope("Get_Contact_Points");

    if(use_old_states)
    {
        LOG::SCOPE scope_tn("restoring tn states");
        //restore tn states
        for(int i=0;i<rigid_body_collection.simulated_rigid_body_particles.m;i++)
            collision_callbacks.Swap_State(rigid_body_collection.simulated_rigid_body_particles(i));
        scope_tn.Pop();
    }

    //HASHTABLE<VECTOR<int,2> > found_pairs;

    LOG::SCOPE scope_contacts("creating contacts");
    for(int i=0;i<pairs.m;i++)
    {
        int id_1=pairs(i)(0),id_2=pairs(i)(1);
        ARRAY<TV> locations,normals;
        ARRAY<typename TV::SCALAR> distances;
        RIGID_BODY<TV>& body_1=rigid_body_collection.Rigid_Body(id_1);
        RIGID_BODY<TV>& body_2=rigid_body_collection.Rigid_Body(id_2);
        //LOG::SCOPE scope_apip("all particles in proximity");
        PARTICLES_IN_PROXIMITY::All_Particles_In_Proximity(body_1,body_2,locations,normals,distances,contact_proximity,stagger_points);
        //scope_apip.Pop();

        for(int j=0;j<locations.m;j++)
            contacts.Append(CONTACT<TV>(body_1,body_2,locations(j),normals(j),distances(j),dt));

        //if(locations.m)
        //{
        //    found_pairs.Set(pairs(i).Sorted());
        //}
    }
    scope_contacts.Pop();

    //ARRAY<VECTOR<int,2> > found_pairs_keys;
    //found_pairs.Get_Keys(found_pairs_keys);
    //LOG::cout << "found_pairs " << found_pairs_keys << std::endl;
    
    if(use_old_states)
    {
        LOG::SCOPE scope_tn1("restoring tn+1 states");
        //restore tn+1 states
        for(int i=0;i<rigid_body_collection.simulated_rigid_body_particles.m;i++)
            collision_callbacks.Swap_State(rigid_body_collection.simulated_rigid_body_particles(i));
        scope_tn1.Pop();
    }
}
//#####################################################################
// Function Push_Out
//#####################################################################
template<class TV>
void Push_Out(RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact)
{
    LOG::SCOPE scope("Projected Gauss Seidel Push_Out");

    typedef typename TV::SCALAR T;

    ARRAY<VECTOR<int,2> > pairs;
    pairs_processed_by_contact.Get_Keys(pairs);

    ARRAY<CONTACT<TV> > contacts;
    Get_Contact_Points(rigid_body_collection,collision_callbacks,pairs,contacts,0,1,false,false);
    for(int i=1;i<contacts.m;i++)
        contacts(i).coefficient_of_friction=0;

    //store linear velocity/update rotational momentum
    //set velocity to 0
    ARRAY<TV> linear_velocities(rigid_body_collection.simulated_rigid_body_particles.m);
    for(int i=0;i<rigid_body_collection.simulated_rigid_body_particles.m;i++)
    {
        RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(rigid_body_collection.simulated_rigid_body_particles(i));
        linear_velocities(i)=body.Twist().linear;
        body.Update_Angular_Momentum();
        body.Twist()=TWIST<TV>();
    }

    T tolerance=(T)1e-3;
    int iteration_maximum=0;

    PROJECTED_GAUSS_SEIDEL::Solve(rigid_body_collection,contacts,tolerance,iteration_maximum);

    for(int i=0;i<rigid_body_collection.simulated_rigid_body_particles.m;i++)
        if(!rigid_body_collection.rigid_body_particles.kinematic(rigid_body_collection.simulated_rigid_body_particles(i)))
            collision_callbacks.Euler_Step_Position(i,1,0); //should pass the actual time.

    //restore linear velocity/update rotational velocity
    for(int i=0;i<rigid_body_collection.simulated_rigid_body_particles.m;i++)
    {
        RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(rigid_body_collection.simulated_rigid_body_particles(i));
        body.Twist().linear=linear_velocities(i);
        body.Update_Angular_Velocity();
    }
}
//#####################################################################

#define INSTANTIATION_HELPER(T,d) \
    template void Solve<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        RIGID_BODY_COLLECTION<VECTOR<T,d> >& rigid_body_collection,RIGID_BODY_COLLISION_PARAMETERS<VECTOR<T,d> >& parameters,const bool correct_contact_energy,const bool use_saved_pairs, \
        const T dt,const T time,ARRAY<TWIST<VECTOR<T,d> > >& mpi_rigid_velocity_save,ARRAY<VECTOR<T,d>::SPIN>& mpi_rigid_angular_momentum_save); \
    template bool Update_Analytic_Multibody_Contact<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        HASHTABLE<VECTOR<std::string,2>,ANALYTICS<VECTOR<T,d> >::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry, \
        const int id_1,const int id_2,MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<T,d> >& multibody,IMPLICIT_OBJECT<VECTOR<T,d> >& levelset,const bool correct_contact_energy, \
        const int max_iterations,const T epsilon_scale,const T dt,const T time,const bool mpi_one_ghost); \
    template bool Update_Analytic_Multibody_Contact<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        HASHTABLE<VECTOR<std::string,2>,ANALYTICS<VECTOR<T,d> >::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry, \
        RIGID_BODY<VECTOR<T,d> >& body0,RIGID_BODY<VECTOR<T,d> >& body1,const bool correct_contact_energy,const int max_iterations, \
        const T epsilon_scale,const T dt,const T time,const bool mpi_one_ghost); \
    template bool Update_Contact_Pair<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        HASHTABLE<VECTOR<std::string,2>,ANALYTICS<VECTOR<T,d> >::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry,const int id_1,const int id_2,const bool correct_contact_energy, \
        const int max_iterations,const T epsilon_scale,const T dt,const T time,const bool mpi_one_ghost); \
    template void Update_Contact_Pair_Helper<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collection,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        const int id_1,const int id_2,const T dt,const T time,const T epsilon_scale,const VECTOR<T,d> & collision_location,const VECTOR<T,d> & collision_normal, \
        const VECTOR<T,d> & collision_relative_velocity,const bool correct_contact_energy,const bool rolling_friction,const bool mpi_one_ghost); \
    template void Euler_Step_Position<VECTOR<T,d> >(RIGID_BODY_COLLISIONS<VECTOR<T,d> >& rigid_body_collisions,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        const int id,const T dt,const T time); \
    template void Register_Analytic_Contacts<VECTOR<T,d> >(HASHTABLE<VECTOR<std::string,2>,ANALYTICS<VECTOR<T,d> >::UPDATE_ANALYTIC_CONTACT_PAIR_T>& analytic_contact_registry); \
    template void Push_Out(RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks,RIGID_BODY_COLLECTION<VECTOR<T,d> >& rigid_body_collection,RIGID_BODY_COLLISION_PARAMETERS<VECTOR<T,d> >& parameters,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact); \
    template bool Solve_Projected_Gauss_Seidel<VECTOR<T,d> >(RIGID_BODY_COLLECTION<VECTOR<T,d> >& rigid_body_collection,RIGIDS_COLLISION_CALLBACKS<VECTOR<T,d> >& collision_callbacks, \
        ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& pairs_processed_by_contact,T desired_separation_distance,T contact_proximity,T dt,T tolerance,int iteration_maximum);

INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(float,3)
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
INSTANTIATION_HELPER(double,3)
    }
}
