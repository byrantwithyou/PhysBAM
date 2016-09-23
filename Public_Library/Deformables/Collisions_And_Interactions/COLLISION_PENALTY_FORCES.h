//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Joseph Teran, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_PENALTY_FORCES
//#####################################################################
#ifndef __COLLISION_PENALTY_FORCES__
#define __COLLISION_PENALTY_FORCES__

#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV> class COLLISION_BODY_COLLECTION;

template<class TV_input>
class COLLISION_PENALTY_FORCES:public DEFORMABLES_FORCES<TV_input>
{
    typedef TV_input TV;typedef typename TV::SCALAR T;
    typedef DEFORMABLES_FORCES<TV_input> BASE;
public:
    using BASE::particles;

    COLLISION_BODY_COLLECTION<TV>* collision_body_list;
    ARRAY<bool,COLLISION_GEOMETRY_ID> skip_collision_body;
    ARRAY<int> check_collision; // TODO: If Append is being called on this, why isn't this a ARRAY?
    ARRAY<TV> collision_force;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > collision_force_derivative;
    T stiffness;
    T separation_parameter;
    T self_collision_reciprocity_factor;
    COLLISION_GEOMETRY_ID collision_body_list_id;
    bool perform_self_collision;

    COLLISION_PENALTY_FORCES(DEFORMABLE_PARTICLES<TV>& particles);
    virtual ~COLLISION_PENALTY_FORCES();

    void Set_Perform_Self_Collision(const bool perform_self_collision_input=true)
    {perform_self_collision=perform_self_collision_input;}

    void Set_Stiffness(const T stiffness_input=(T)1e4)
    {stiffness=stiffness_input;}

    void Set_Separation_Parameter(const T separation_parameter_input=(T)1e-4)
    {separation_parameter=separation_parameter_input;}

    void Set_Collision_Body_List(COLLISION_BODY_COLLECTION<TV>& collision_body_list_input);

    template<class T_MESH>
    void Set_Boundary_Only_Collisions(T_MESH& mesh)
    {ARRAY<bool>* old_node_on_boundary=mesh.node_on_boundary;mesh.node_on_boundary=0;
    mesh.Initialize_Node_On_Boundary();for(int p=0;p<particles.Size();p++) if((*mesh.node_on_boundary)(p)) check_collision.Append(p);
    collision_force.Resize(check_collision.m);collision_force_derivative.Resize(check_collision.m);
    delete mesh.node_on_boundary;mesh.node_on_boundary=old_node_on_boundary;}

    void Set_Collision_Body_List_ID(const COLLISION_GEOMETRY_ID id)
    {collision_body_list_id=id;}

    void Set_Self_Collision_Reciprocity_Factor(const T self_collision_reciprocity_factor_input=(T)2)
    {self_collision_reciprocity_factor=self_collision_reciprocity_factor_input;}

    void Check_Collision(const int particle)
    {check_collision.Append_Unique(particle);} // TODO: not very efficient.

    void Omit_Collision(const int particle)
    {int index;if(check_collision.Find(particle,index)) check_collision.Remove_Index(index);} // TODO: not very efficient.

    void Resize_Collision_Arrays_From_Check_Collision()
    {collision_force.Resize(check_collision.m);collision_force_derivative.Resize(check_collision.m);}

    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Update_Forces_And_Derivatives();
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    T CFL_Strain_Rate() const override;
//#####################################################################
};
}
#endif
