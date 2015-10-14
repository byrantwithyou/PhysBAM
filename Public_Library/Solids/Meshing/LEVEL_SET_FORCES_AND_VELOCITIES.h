//#####################################################################
// Copyright 2002-2007, Christopher Allocco, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVEL_SET_FORCES_AND_VELOCITIES
//#####################################################################
#ifndef __LEVEL_SET_FORCES_AND_VELOCITIES__
#define __LEVEL_SET_FORCES_AND_VELOCITIES__

#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
namespace PhysBAM{

template<class TV>
class LEVEL_SET_FORCES_AND_VELOCITIES:public EXAMPLE_FORCES_AND_VELOCITIES<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_MESH_OBJECT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_BOUNDARY_OBJECT;
    typedef EXAMPLE_FORCES_AND_VELOCITIES<TV> BASE;
    using BASE::Add_External_Forces;using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual
public:
    T_MESH_OBJECT& mesh_object;
    IMPLICIT_OBJECT<TV>& implicit_object;
    T force_attraction_coefficient,velocity_attraction_coefficient;
    bool use_external_forces,use_external_velocities;
    bool use_external_velocities_normal_to_boundary;
    bool allow_tangential_velocity_slip;

    LEVEL_SET_FORCES_AND_VELOCITIES(T_MESH_OBJECT& mesh_object,IMPLICIT_OBJECT<TV>& implicit_object);
    ~LEVEL_SET_FORCES_AND_VELOCITIES();

    void Set_Force_Attraction_Coefficient(const T coefficient=.1)
    {force_attraction_coefficient=coefficient;}

    void Set_Velocity_Attraction_Coefficient(const T coefficient=.1)
    {velocity_attraction_coefficient=coefficient;}

    void Use_External_Forces()
    {use_external_forces=true;use_external_velocities=false;}

    void Use_External_Velocities()
    {use_external_forces=false;use_external_velocities=true;}

    void Use_External_Velocities_Normal_To_Boundary()
    {use_external_velocities_normal_to_boundary=true;}

    void Use_External_Velocities_Towards_Zero_Isocontour()
    {use_external_velocities_normal_to_boundary=false;}

    void Allow_Tangential_Velocity_Slip(const bool allow_tangential_velocity_slip_input=true)
    {allow_tangential_velocity_slip=allow_tangential_velocity_slip_input;}

    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override;
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override;
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override;
};
}
#endif


