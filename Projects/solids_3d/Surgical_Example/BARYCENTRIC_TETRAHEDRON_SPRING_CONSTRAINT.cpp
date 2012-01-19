//#####################################################################
// Copyright 2006, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT  
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include "BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT.h"
#include <float.h>
using namespace PhysBAM;

//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>::
Initialize_Weights_Outer_Product()
{
    weights_outer_product.Resize(constrained_location_weights.m);
    for(int t=0;t<constrained_location_weights.m;t++){
        T w1=constrained_location_weights(t).x,w2=constrained_location_weights(t).y,w3=constrained_location_weights(t).z,w4=(T)1-(w1+w2+w3);
        weights_outer_product(t).x[0]=w1*w1;weights_outer_product(t).x[1]=w2*w1;weights_outer_product(t).x[2]=w3*w1;weights_outer_product(t).x[3]=w4*w1;
        weights_outer_product(t).x[4]=w1*w2;weights_outer_product(t).x[5]=w2*w2;weights_outer_product(t).x[6]=w3*w2;weights_outer_product(t).x[7]=w4*w2;
        weights_outer_product(t).x[8]=w1*w3;weights_outer_product(t).x[9]=w2*w3;weights_outer_product(t).x[10]=w3*w3;weights_outer_product(t).x[11]=w4*w3;
        weights_outer_product(t).x[12]=w1*w4;weights_outer_product(t).x[13]=w2*w4;weights_outer_product(t).x[14]=w3*w4;weights_outer_product(t).x[15]=w4*w4;}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{

}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>::
Add_Velocity_Independent_Forces(ARRAY<VECTOR<T,3> >& F,const T time) const
{
    ARRAY<VECTOR<T,3> >& X=particles.X.array;
    VECTOR<T,3> weights,force_on_embedded_node;
    if(!youngs_modulus.m){
        for(int t=0;t<constrained_tets.m;t++){
            int i,j,k,l;
            tetrahedron_mesh.elements.Get(constrained_tets(t),i,j,k,l);
            weights=constrained_location_weights(t);
            T w4=((T)1-(weights.x+weights.y+weights.z));
            force_on_embedded_node=-constant_youngs_modulus*(weights.x*X(i)+weights.y*X(j)+weights.z*X(k)+w4*X(l)-constrained_node_locations(t));
            F(i)+=weights.x*force_on_embedded_node;
            F(j)+=weights.y*force_on_embedded_node;
            F(k)+=weights.z*force_on_embedded_node;
            F(l)+=w4*force_on_embedded_node;}}
    else{
        for(int t=0;t<constrained_tets.m;t++){
            int i,j,k,l;
            tetrahedron_mesh.elements.Get(constrained_tets(t),i,j,k,l);
            weights=constrained_location_weights(t);
            T w4=((T)1-(weights.x+weights.y+weights.z));
            force_on_embedded_node=-youngs_modulus(t)*(weights.x*X(i)+weights.y*X(j)+weights.z*X(k)+w4*X(l)-constrained_node_locations(t));
            F(i)+=weights.x*force_on_embedded_node;
            F(j)+=weights.y*force_on_embedded_node;
            F(k)+=weights.z*force_on_embedded_node;
            F(l)+=w4*force_on_embedded_node;}}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>::
Add_Force_Differential(const ARRAY<VECTOR<T,3> >& dX,ARRAY<VECTOR<T,3> >& dF,const T time) const
{
    if(!youngs_modulus.m){
        for(int t=0;t<constrained_tets.m;t++){
            int i,j,k,l;tetrahedron_mesh.elements.Get(constrained_tets(t),i,j,k,l);
            dF(i)+=-constant_youngs_modulus*(weights_outer_product(t).x[0]*dX(i)+weights_outer_product(t).x[1]*dX(j)+weights_outer_product(t).x[2]*dX(k)+weights_outer_product(t).x[3]*dX(l));
            dF(j)+=-constant_youngs_modulus*(weights_outer_product(t).x[4]*dX(i)+weights_outer_product(t).x[5]*dX(j)+weights_outer_product(t).x[6]*dX(k)+weights_outer_product(t).x[7]*dX(l));
            dF(k)+=-constant_youngs_modulus*(weights_outer_product(t).x[8]*dX(i)+weights_outer_product(t).x[9]*dX(j)+weights_outer_product(t).x[10]*dX(k)+weights_outer_product(t).x[11]*dX(l));
            dF(l)+=-constant_youngs_modulus*(weights_outer_product(t).x[12]*dX(i)+weights_outer_product(t).x[13]*dX(j)+weights_outer_product(t).x[14]*dX(k)+weights_outer_product(t).x[15]*dX(l));}}
    else{
        for(int t=0;t<constrained_tets.m;t++){
            int i,j,k,l;tetrahedron_mesh.elements.Get(constrained_tets(t),i,j,k,l);
            dF(i)+=-youngs_modulus(t)*(weights_outer_product(t).x[0]*dX(i)+weights_outer_product(t).x[1]*dX(j)+weights_outer_product(t).x[2]*dX(k)+weights_outer_product(t).x[3]*dX(l));
            dF(j)+=-youngs_modulus(t)*(weights_outer_product(t).x[4]*dX(i)+weights_outer_product(t).x[5]*dX(j)+weights_outer_product(t).x[6]*dX(k)+weights_outer_product(t).x[7]*dX(l));
            dF(k)+=-youngs_modulus(t)*(weights_outer_product(t).x[8]*dX(i)+weights_outer_product(t).x[9]*dX(j)+weights_outer_product(t).x[10]*dX(k)+weights_outer_product(t).x[11]*dX(l));
            dF(l)+=-youngs_modulus(t)*(weights_outer_product(t).x[12]*dX(i)+weights_outer_product(t).x[13]*dX(j)+weights_outer_product(t).x[14]*dX(k)+weights_outer_product(t).x[15]*dX(l));}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>::
Add_Velocity_Dependent_Forces(const ARRAY<VECTOR<T,3> >& V,ARRAY<VECTOR<T,3> >& F,const T time) const
{
    if(!damping.m){
        for(int t=0;t<constrained_tets.m;t++){
            int i,j,k,l;tetrahedron_mesh.elements.Get(constrained_tets(t),i,j,k,l);
            F(i)+=-constant_damping*(weights_outer_product(t).x[0]*V(i)+weights_outer_product(t).x[1]*V(j)+weights_outer_product(t).x[2]*V(k)+weights_outer_product(t).x[3]*V(l));
            F(j)+=-constant_damping*(weights_outer_product(t).x[4]*V(i)+weights_outer_product(t).x[5]*V(j)+weights_outer_product(t).x[6]*V(k)+weights_outer_product(t).x[7]*V(l));
            F(k)+=-constant_damping*(weights_outer_product(t).x[8]*V(i)+weights_outer_product(t).x[9]*V(j)+weights_outer_product(t).x[10]*V(k)+weights_outer_product(t).x[11]*V(l));
            F(l)+=-constant_damping*(weights_outer_product(t).x[12]*V(i)+weights_outer_product(t).x[13]*V(j)+weights_outer_product(t).x[14]*V(k)+weights_outer_product(t).x[15]*V(l));}}
    else{
        for(int t=0;t<constrained_tets.m;t++){
            int i,j,k,l;tetrahedron_mesh.elements.Get(constrained_tets(t),i,j,k,l);
            F(i)+=-damping(t)*(weights_outer_product(t).x[0]*V(i)+weights_outer_product(t).x[1]*V(j)+weights_outer_product(t).x[2]*V(k)+weights_outer_product(t).x[3]*V(l));
            F(j)+=-damping(t)*(weights_outer_product(t).x[4]*V(i)+weights_outer_product(t).x[5]*V(j)+weights_outer_product(t).x[6]*V(k)+weights_outer_product(t).x[7]*V(l));
            F(k)+=-damping(t)*(weights_outer_product(t).x[8]*V(i)+weights_outer_product(t).x[9]*V(j)+weights_outer_product(t).x[10]*V(k)+weights_outer_product(t).x[11]*V(l));
            F(l)+=-damping(t)*(weights_outer_product(t).x[12]*V(i)+weights_outer_product(t).x[13]*V(j)+weights_outer_product(t).x[14]*V(k)+weights_outer_product(t).x[15]*V(l));}}}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>::
Initialize_CFL()
{
    CFL_elastic_time_step_of_fragment.Append(FLT_MAX);CFL_damping_time_step_of_fragment.Append(FLT_MAX);
    CFL_initialized_of_fragment.Append(true);
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<T>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
template class BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BARYCENTRIC_TETRAHEDRON_SPRING_CONSTRAINT<double>;
#endif
