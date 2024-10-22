//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_STRETCH_FORCES
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/MATRIX.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Deformables/Forces/BW_STRETCH_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <cfloat>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BW_STRETCH_FORCES<TV>::
BW_STRETCH_FORCES(DEFORMABLE_PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh_input,const T stiffness_coefficient_input,const T damping_coefficient_input)
    :BW_MATERIAL_SPACE_FORCES<TV,2>(particles,triangle_mesh_input,stiffness_coefficient_input,damping_coefficient_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BW_STRETCH_FORCES<TV>::
~BW_STRETCH_FORCES()
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void BW_STRETCH_FORCES<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    for(int s:force_simplices){
        typename BASE::STATE& state=states(s);
        typename BASE::MATERIAL_FORCE_STATE& material_force_state=material_force_states(s);
        Compute_UV_Deformation(s);

        state.C=material_force_state.rest_state_triangle_area*VECTOR<T,2>(material_force_state.w_u_magnitude-material_force_state.b_u,
            material_force_state.w_v_magnitude-material_force_state.b_v);
        state.C_dot=VECTOR<T,2>();
        for(int i=0;i<3;i++){
            state.dC_dx(i)=MATRIX<T,3,2>(material_force_state.rest_state_triangle_area*material_force_state.dwu_dx(i)*material_force_state.w_u,
                material_force_state.rest_state_triangle_area*material_force_state.dwv_dx(i)*material_force_state.w_v);
            state.C_dot+=state.dC_dx(i).Transposed()*particles.V(state.nodes[i]);}
        for(int i=0;i<3;i++) for(int j=0;j<3;j++){
            MATRIX<T,3> dCu_dxi_dxj=(material_force_state.rest_state_triangle_area/material_force_state.w_u_magnitude)*material_force_state.dwu_dx(i)*
                material_force_state.dwu_dx(j)*(MATRIX<T,3>::Identity_Matrix()-MATRIX<T,3>::Outer_Product(material_force_state.w_u,material_force_state.w_u));
            MATRIX<T,3> dCv_dxi_dxj=(material_force_state.rest_state_triangle_area/material_force_state.w_v_magnitude)*material_force_state.dwv_dx(i)*
                material_force_state.dwv_dx(j)*(MATRIX<T,3>::Identity_Matrix()-MATRIX<T,3>::Outer_Product(material_force_state.w_v,material_force_state.w_v));
            state.dC_dxi_dxj_times_C(i,j)=dCu_dxi_dxj*state.C(0)+dCv_dxi_dxj*state.C(1);
            state.dC_dxi_dxj_times_C_dot(i,j)=dCu_dxi_dxj*state.C_dot(0)+dCv_dxi_dxj*state.C_dot(1);}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BW_STRETCH_FORCES<TV>::
Potential_Energy(int s,const T time) const
{
    return (T).5*states(s).stiffness_coefficient*states(s).C.Magnitude_Squared();
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BW_STRETCH_FORCES<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    const_cast<BW_STRETCH_FORCES<TV>* >(this)->Update_Position_Based_State(time,true,true);
    for(int s:force_simplices){
        potential_energy+=Potential_Energy(s,time);}
    return potential_energy;
}
//#####################################################################
// Function Create_BW_Stretch_Force
//#####################################################################
template<class TV> BW_STRETCH_FORCES<TV>* PhysBAM::
Create_BW_Stretch_Force(DEFORMABLE_PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh,const typename TV::SCALAR stiffness_coefficient_input,const typename TV::SCALAR damping_coefficient_input)
{
    BW_STRETCH_FORCES<TV>* sf=new BW_STRETCH_FORCES<TV>(particles,triangle_mesh,stiffness_coefficient_input,damping_coefficient_input);
    return sf;
}
//#####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(T) \
    template class BW_STRETCH_FORCES<VECTOR<T,3> >; \
    template BW_STRETCH_FORCES<VECTOR<T,3> >* Create_BW_Stretch_Force<VECTOR<T,3> >(DEFORMABLE_PARTICLES<VECTOR<T,3> >&,TRIANGLE_MESH&,const T,const T);

INSTANTIATION_HELPER(float)
INSTANTIATION_HELPER(double)
}
