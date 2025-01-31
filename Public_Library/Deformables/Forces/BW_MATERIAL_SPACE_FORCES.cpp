//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_MATERIAL_SPACE_FORCES
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/MATRIX_3X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Deformables/Forces/BW_MATERIAL_SPACE_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <cfloat>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> BW_MATERIAL_SPACE_FORCES<TV,d>::
BW_MATERIAL_SPACE_FORCES(DEFORMABLE_PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh_input,const T stiffness_coefficient_input,const T damping_coefficient_input)
    :BW_FORCES<TV,d,3>(particles,triangle_mesh_input,stiffness_coefficient_input,damping_coefficient_input)
{
    // Save initial state
    states.Resize(triangle_mesh.elements.m,no_init);
    material_force_states.Resize(triangle_mesh.elements.m,no_init);
    ARRAY<bool> particle_is_simulated(particles.Size());particle_is_simulated.Fill(true);
    Update_Mpi(particle_is_simulated,0);
    for(int s:force_simplices){
        int node1,node2,node3;triangle_mesh.elements(s).Get(node1,node2,node3);
        typename BASE::STATE& state=states(s);
        state.nodes=VECTOR<int,3>(node1,node2,node3);
        state.stiffness_coefficient=stiffness_coefficient_input;
        state.damping_coefficient=damping_coefficient_input;
        MATERIAL_FORCE_STATE& material_force_state=material_force_states(s);
        // TODO this will not work if the cloth does not start out x-z aligned
        material_force_state.delta_u=VECTOR<T,2>(particles.X(node2).x-particles.X(node1).x,particles.X(node3).x-particles.X(node1).x);
        material_force_state.delta_v=VECTOR<T,2>(particles.X(node2).z-particles.X(node1).z,particles.X(node3).z-particles.X(node1).z);
        material_force_state.rest_state_triangle_area=TRIANGLE_3D<T>(particles.X.Subset(state.nodes)).Area();
        // scale invariant
        material_force_state.rest_state_triangle_area=sqrt(material_force_state.rest_state_triangle_area);
        material_force_state.rest_state_triangle_area=material_force_state.rest_state_triangle_area*sqrt(material_force_state.rest_state_triangle_area);
    }
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> BW_MATERIAL_SPACE_FORCES<TV,d>::
~BW_MATERIAL_SPACE_FORCES()
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV,int d> void BW_MATERIAL_SPACE_FORCES<TV,d>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    triangle_mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV,int d> void BW_MATERIAL_SPACE_FORCES<TV,d>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    Update_Force_Elements(force_simplices,triangle_mesh.elements,particle_is_simulated);
}
//#####################################################################
// Function Compute_UV_Deformation
//#####################################################################
template<class TV,int d> void BW_MATERIAL_SPACE_FORCES<TV,d>::
Compute_UV_Deformation(const int c)
{
    typename BASE::STATE& state=states(c);
    MATERIAL_FORCE_STATE& material_force_state=material_force_states(c);
    material_force_state.denom=(material_force_state.delta_u(0)*material_force_state.delta_v(1)-material_force_state.delta_u(1)*material_force_state.delta_v(0));
    material_force_state.w_u=((particles.X(state.nodes(1))-particles.X(state.nodes(0)))*material_force_state.delta_v(1)-
        (particles.X(state.nodes(2))-particles.X(state.nodes(0)))*material_force_state.delta_v(0))/material_force_state.denom;
    material_force_state.w_u_magnitude=material_force_state.w_u.Normalize();
    material_force_state.w_v=(-(particles.X(state.nodes(1))-particles.X(state.nodes(0)))*material_force_state.delta_u(1)+
        (particles.X(state.nodes(2))-particles.X(state.nodes(0)))*material_force_state.delta_u(0))/material_force_state.denom;
    material_force_state.w_v_magnitude=material_force_state.w_v.Normalize();

    material_force_state.dwu_dx=VECTOR<T,3>((material_force_state.delta_v(0)-material_force_state.delta_v(1))/material_force_state.denom,
        material_force_state.delta_v(1)/material_force_state.denom,-material_force_state.delta_v(0)/material_force_state.denom);
    material_force_state.dwv_dx=VECTOR<T,3>((material_force_state.delta_u(1)-material_force_state.delta_u(0))/material_force_state.denom,
        -material_force_state.delta_u(1)/material_force_state.denom,material_force_state.delta_u(0)/material_force_state.denom);
}
//#####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(T,d) \
    template class BW_MATERIAL_SPACE_FORCES<VECTOR<T,3>,d>;

INSTANTIATION_HELPER(float,1)
INSTANTIATION_HELPER(float,2)
INSTANTIATION_HELPER(double,1)
INSTANTIATION_HELPER(double,2)
}
