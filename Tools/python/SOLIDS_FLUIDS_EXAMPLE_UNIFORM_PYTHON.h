//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLIDS_FLUIDS_EXAMPLE_UNIFORM_PYTHON__
#define __SOLIDS_FLUIDS_EXAMPLE_UNIFORM_PYTHON__
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <boost/python.hpp>
#include <PhysBAM_Solids/PhysBAM_Solids/Fragments/PARTICLE_CONNECTIVITY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS_ND.h>

#define PYTHON_CALL_OVERRIDE(...) \
    if(boost::python::override f=this->get_override(__FUNCTION__)) f(__VA_ARGS__); \
    else PHYSBAM_WARN_IF_NOT_OVERRIDDEN();

namespace PhysBAM{

template<class T_GRID>
class SOLIDS_FLUIDS_EXAMPLE_UNIFORM_PYTHON:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>,public boost::python::wrapper<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> >
{
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> BASE;
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::SCALAR T;
public:
    using BASE::solids_parameters;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM_PYTHON(const STREAM_TYPE stream_type,const int number_of_regions,const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type)
        :BASE(stream_type,number_of_regions,type)
    {}

    void Initialize_Bodies() PHYSBAM_OVERRIDE
    {if(boost::python::override f=this->get_override("Initialize_Bodies"))f();
    else PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
        LOG::cout<<"CFL IS "<<solids_parameters.cfl<<std::endl;
    BASE::Initialize_Bodies();}

    void Set_External_Velocities(RAW_ARRAY<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
    {//BASE::Set_External_Velocities(V,velocity_time,current_position_time);
    PYTHON_CALL_OVERRIDE(boost::ref(V),velocity_time,current_position_time);}

    void Zero_Out_Enslaved_Velocity_Nodes(RAW_ARRAY<TV> V,const T velocity_time,const T current_position_time,const SUPER_FRAGMENT<TV>& super_fragment) PHYSBAM_OVERRIDE
    {//BASE::Zero_Out_Enslaved_Velocity_Nodes(V,velocity_time,current_position_time,super_fragment);
    PYTHON_CALL_OVERRIDE(boost::ref(V),velocity_time,current_position_time,boost::ref(super_fragment));}

    void Set_External_Velocities(RAW_ARRAY<TWIST<TV> > twist,const T velocity_time,const T current_position_time,const SUPER_FRAGMENT<TV>& super_fragment)
    {}//    PYTHON_CALL_OVERRIDE(boost::ref(twist),velocity_time,current_position_time,super_fragment);
/*
    void Zero_Out_Enslaved_Velocity_Nodes(RAW_ARRAY<TWIST<TV> > twist,const T velocity_time,const T current_position_time,const SUPER_FRAGMENT<TV>& super_fragment)
    {BASE::Zero_Out_Enslaved_Velocity_Nodes(twist,velocity_time,current_position_time,super_fragment);
    PYTHON_CALL_OVERRIDE(boost::ref(twist),velocity_time,current_position_time,super_fragment);}
*/
    void Zero_Out_Enslaved_Velocity_Nodes(RAW_ARRAY<TWIST<TV> > twist,const T velocity_time,const T current_position_time,const SUPER_FRAGMENT<TV>& super_fragment) // or zero out components of their velocities
    {}

    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(frame);}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(frame);}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(time);}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(time,substep);}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(time,substep);}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(dt,time);}
    void Add_External_Forces(RAW_ARRAY<TV> F,const T time) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(F),time);}
    void Add_External_Forces(TWIST<TV>& wrench,const T time,const int id) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(wrench),time,id);}
    void Add_External_Forces(RAW_ARRAY<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(wrench),time);}
    void Add_External_Impulses_Before(RAW_ARRAY<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(V),time,dt);}
    void Add_External_Impulses(RAW_ARRAY<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(V),time,dt);}
    void Add_External_Impulse(RAW_ARRAY<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(V),node,time,dt);}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(frame),time,id);}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {
        if(boost::python::override f=this->get_override("Set_Kinematic_Velocities")) return f(boost::ref(twist),time,id);
        else return false;}

    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
    {if(boost::python::override f=this->get_override("Limit_Solids_Dt")){
        boost::python::object o=f(time);
        if(o) dt=min(dt,boost::python::extract<T>(o)());}
    else PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE();}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(time);}
    void Set_External_Positions(RAW_ARRAY<TV> X,const T time) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(X),time);}
    void Zero_Out_Enslaved_Position_Nodes(RAW_ARRAY<TV> X,const T time,const SUPER_FRAGMENT<TV>& super_fragment) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(X),time,boost::ref(super_fragment));}

    void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(time);}

    void Postprocess_Fragments() PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE();}
    void Postprocess_Super_Fragments(const ARRAY<PAIR<SUPER_FRAGMENT_ID,SUPER_FRAGMENT_ID> >& swap_pairs,const ARRAY<SUPER_FRAGMENT_ID>& rebuild,
        SUPER_FRAGMENT_ID old_max) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(swap_pairs),boost::ref(rebuild),old_max);}
    void Update_Super_Fragments(const ARRAY<PAIR<SUPER_FRAGMENT_ID,SUPER_FRAGMENT_ID> >& swap_pairs,const ARRAY<SUPER_FRAGMENT_ID>& rebuild,
        SUPER_FRAGMENT_ID old_max) PHYSBAM_OVERRIDE {PYTHON_CALL_OVERRIDE(boost::ref(swap_pairs),boost::ref(rebuild),old_max);}
    void Add_External_Fragment_Connectivity(PARTICLE_CONNECTIVITY<TV>& particle_connectivity,ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE
    {PYTHON_CALL_OVERRIDE(boost::ref(particle_connectivity),boost::ref(particle_is_simulated));}
    void Add_External_Super_Fragment_Connectivity(SPARSE_UNION_FIND<FRAGMENT_ID>& union_find) PHYSBAM_OVERRIDE
    {PYTHON_CALL_OVERRIDE(boost::ref(union_find));}

    void Collide_All_Bodies()
    {// add structures and rigid bodies to collisions
    solids_parameters.solid_body_collection.deformable_body_collection.collisions.collision_structures.Append_Elements(solids_parameters.solid_body_collection.deformable_body_collection.deformable_geometry.structures);
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.solid_body_collection.rigid_body_collection.rigid_geometry_collection);
    solids_parameters.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(solids_parameters.solid_body_collection.deformable_body_collection.deformable_geometry.structures);}

    void Collide_Rigid_Bodies()
    {solids_parameters.collision_body_list.Add_Bodies(solids_parameters.solid_body_collection.rigid_body_collection.rigid_geometry_collection);}

    void Collide_Segmented_Curve(SEGMENTED_CURVE<TV>& curve)
    {// add structures and rigid bodies to collisions
    solids_parameters.solid_body_collection.deformable_body_collection.collisions.collision_structures.Append(&curve);
    solids_parameters.triangle_repulsions_and_collisions_geometry.structures.Append(&curve);}

    void Unite_All_Fragments(PARTICLE_CONNECTIVITY<TV>& particle_connectivity,ARRAY<bool>& particle_is_simulated)
    {for(int i=2;i<=particle_connectivity.particles_number;i++) particle_connectivity.Union(1,i);}
        
    void Zero_Nodes(const ARRAY<int>& indices,RAW_ARRAY<TV> V)
    {for(int i=0;i<indices.m;i++) V(indices(i))=TV();}

    void Set_Rigid_Velocities(const ARRAY<int>& indices,const ARRAY<TV>& object_positions,const RIGID_BODY<TV>& rigid,RAW_ARRAY<TV> V)
    {for(int i=0;i<indices.m;i++) V(indices(i))=rigid.Pointwise_Object_Velocity(rigid.World_Space_Point(object_positions(i)));}

    void Set_Rigid_Positions(const ARRAY<int>& indices,const ARRAY<TV>& object_positions,const RIGID_BODY<TV>& rigid,RAW_ARRAY<TV> X)
    {for(int i=0;i<indices.m;i++) X(indices(i))=rigid.World_Space_Point(object_positions(i));}

    // TODO: Do other callbacks for zero/set, etc.
    // TODO: Make not overriden case simpler

//#####################################################################
};
}
#endif
