//#####################################################################
// Copyright 2006, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING  
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include "BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING.h"
#include <float.h>
using namespace PhysBAM;

//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>::
Initialize_Weights_Outer_Product()
{
    weights_outer_product.Resize(constrained_location_weights.m);
    for(int t=1;t<=constrained_location_weights.m;t++){
        T wa_1=constrained_location_weights(1,t).x,wa_2=constrained_location_weights(1,t).y,wa_3=constrained_location_weights(1,t).z,wa_4=(T)1-(wa_1+wa_2+wa_3);
        T wb_1=constrained_location_weights(2,t).x,wb_2=constrained_location_weights(2,t).y,wb_3=constrained_location_weights(2,t).z,wb_4=(T)1-(wb_1+wb_2+wb_3);
        //1
        weights_outer_product(1,t).x[0]=wa_1*wa_1;weights_outer_product(1,t).x[1]=wa_2*wa_1;weights_outer_product(1,t).x[2]=wa_3*wa_1;weights_outer_product(1,t).x[3]=wa_4*wa_1;
        weights_outer_product(1,t).x[4]=wa_1*wa_2;weights_outer_product(1,t).x[5]=wa_2*wa_2;weights_outer_product(1,t).x[6]=wa_3*wa_2;weights_outer_product(1,t).x[7]=wa_4*wa_2;
        weights_outer_product(1,t).x[8]=wa_1*wa_3;weights_outer_product(1,t).x[9]=wa_2*wa_3;weights_outer_product(1,t).x[10]=wa_3*wa_3;weights_outer_product(1,t).x[11]=wa_4*wa_3;
        weights_outer_product(1,t).x[12]=wa_1*wa_4;weights_outer_product(1,t).x[13]=wa_2*wa_4;weights_outer_product(1,t).x[14]=wa_3*wa_4;weights_outer_product(1,t).x[15]=wa_4*wa_4;
        //2
        weights_outer_product(2,t).x[0]=wb_1*wa_1;weights_outer_product(2,t).x[1]=wb_2*wa_1;weights_outer_product(1,t).x[2]=wb_3*wa_1;weights_outer_product(1,t).x[3]=wb_4*wa_1;
        weights_outer_product(2,t).x[4]=wb_1*wa_2;weights_outer_product(1,t).x[5]=wb_2*wa_2;weights_outer_product(1,t).x[6]=wb_3*wa_2;weights_outer_product(1,t).x[7]=wb_4*wa_2;
        weights_outer_product(2,t).x[8]=wb_1*wa_3;weights_outer_product(1,t).x[9]=wb_2*wa_3;weights_outer_product(1,t).x[10]=wb_3*wa_3;weights_outer_product(1,t).x[11]=wb_4*wa_3;
        weights_outer_product(2,t).x[12]=wb_1*wa_4;weights_outer_product(1,t).x[13]=wb_2*wa_4;weights_outer_product(1,t).x[14]=wb_3*wa_4;weights_outer_product(1,t).x[15]=wb_4*wa_4;
        //3
        weights_outer_product(3,t).x[0]=wa_1*wb_1;weights_outer_product(2,t).x[1]=wa_2*wb_1;weights_outer_product(1,t).x[2]=wa_3*wb_1;weights_outer_product(1,t).x[3]=wa_4*wb_1;
        weights_outer_product(3,t).x[4]=wa_1*wb_2;weights_outer_product(1,t).x[5]=wa_2*wb_2;weights_outer_product(1,t).x[6]=wa_3*wb_2;weights_outer_product(1,t).x[7]=wa_4*wb_2;
        weights_outer_product(3,t).x[8]=wa_1*wb_3;weights_outer_product(1,t).x[9]=wa_2*wb_3;weights_outer_product(1,t).x[10]=wa_3*wb_3;weights_outer_product(1,t).x[11]=wa_4*wb_3;
        weights_outer_product(3,t).x[12]=wa_1*wb_4;weights_outer_product(1,t).x[13]=wa_2*wb_4;weights_outer_product(1,t).x[14]=wa_3*wb_4;weights_outer_product(1,t).x[15]=wa_4*wb_4;
        //4
        weights_outer_product(4,t).x[0]=wb_1*wb_1;weights_outer_product(4,t).x[1]=wb_2*wb_1;weights_outer_product(4,t).x[2]=wb_3*wb_1;weights_outer_product(4,t).x[3]=wb_4*wb_1;
        weights_outer_product(4,t).x[4]=wb_1*wb_2;weights_outer_product(4,t).x[5]=wb_2*wb_2;weights_outer_product(4,t).x[6]=wb_3*wb_2;weights_outer_product(4,t).x[7]=wb_4*wb_2;
        weights_outer_product(4,t).x[8]=wb_1*wb_3;weights_outer_product(4,t).x[9]=wb_2*wb_3;weights_outer_product(4,t).x[10]=wb_3*wb_3;weights_outer_product(4,t).x[11]=wb_4*wb_3;
        weights_outer_product(4,t).x[12]=wb_1*wb_4;weights_outer_product(4,t).x[13]=wb_2*wb_4;weights_outer_product(4,t).x[14]=wb_3*wb_4;weights_outer_product(4,t).x[15]=wb_4*wb_4;}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{

}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>::
Add_Velocity_Independent_Forces(ARRAY<VECTOR<T,3> >& F,const T time) const
{
    ARRAY<VECTOR<T,3> >& X=particles.X.array;
    VECTOR<T,3> weights_a,weights_b,force_on_embedded_node_a,embedded_location_a,embedded_location_b;
    if(!youngs_modulus.m){
        for(int t=1;t<=constrained_tets.m;t++){
            int ia,ja,ka,la,ib,jb,kb,lb;
            tetrahedron_mesh.elements.Get(constrained_tets(1,t),ia,ja,ka,la);
            tetrahedron_mesh.elements.Get(constrained_tets(2,t),ib,jb,kb,lb);
            weights_a=constrained_location_weights(1,t);
            weights_b=constrained_location_weights(2,t);
            T wa_4=((T)1-(weights_a.x+weights_a.y+weights_a.z));
            T wb_4=((T)1-(weights_b.x+weights_b.y+weights_b.z));
            embedded_location_a=weights_a.x*X(ia)+weights_a.y*X(ja)+weights_a.z*X(ka)+wa_4*X(la);
            embedded_location_b=weights_b.x*X(ib)+weights_b.y*X(jb)+weights_b.z*X(kb)+wb_4*X(lb);
            force_on_embedded_node_a=-constant_youngs_modulus*(embedded_location_a-embedded_location_b);
            F(ia)+=weights_a.x*force_on_embedded_node_a;
            F(ja)+=weights_a.y*force_on_embedded_node_a;
            F(ka)+=weights_a.z*force_on_embedded_node_a;
            F(la)+=wa_4*force_on_embedded_node_a;
            F(ib)-=weights_b.x*force_on_embedded_node_a;
            F(jb)-=weights_b.y*force_on_embedded_node_a;
            F(kb)-=weights_b.z*force_on_embedded_node_a;
            F(lb)-=wb_4*force_on_embedded_node_a;}}
    else{
        for(int t=1;t<=constrained_tets.m;t++){
            int ia,ja,ka,la,ib,jb,kb,lb;
            tetrahedron_mesh.elements.Get(constrained_tets(1,t),ia,ja,ka,la);
            tetrahedron_mesh.elements.Get(constrained_tets(2,t),ib,jb,kb,lb);
            weights_a=constrained_location_weights(1,t);
            weights_b=constrained_location_weights(2,t);
            T wa_4=((T)1-(weights_a.x+weights_a.y+weights_a.z));
            T wb_4=((T)1-(weights_b.x+weights_b.y+weights_b.z));
            embedded_location_a=weights_a.x*X(ia)+weights_a.y*X(ja)+weights_a.z*X(ka)+wa_4*X(la);
            embedded_location_b=weights_b.x*X(ib)+weights_b.y*X(jb)+weights_b.z*X(kb)+wb_4*X(lb);
            force_on_embedded_node_a=-youngs_modulus(t)*(embedded_location_a-embedded_location_b);
            F(ia)+=weights_a.x*force_on_embedded_node_a;
            F(ja)+=weights_a.y*force_on_embedded_node_a;
            F(ka)+=weights_a.z*force_on_embedded_node_a;
            F(la)+=wa_4*force_on_embedded_node_a;
            F(ib)-=weights_b.x*force_on_embedded_node_a;
            F(jb)-=weights_b.y*force_on_embedded_node_a;
            F(kb)-=weights_b.z*force_on_embedded_node_a;
            F(lb)-=wb_4*force_on_embedded_node_a;}}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>::
Add_Force_Differential(const ARRAY<VECTOR<T,3> >& dX,ARRAY<VECTOR<T,3> >& dF,const T time) const
{
    if(!youngs_modulus.m){
        for(int t=1;t<=constrained_tets.m;t++){
            int ia,ja,ka,la;tetrahedron_mesh.elements.Get(constrained_tets(1,t),ia,ja,ka,la);
            int ib,jb,kb,lb;tetrahedron_mesh.elements.Get(constrained_tets(2,t),ib,jb,kb,lb);
            //1
            dF(ia)+=-constant_youngs_modulus*(weights_outer_product(1,t).x[0]*dX(ia)+weights_outer_product(1,t).x[1]*dX(ja)+weights_outer_product(1,t).x[2]*dX(ka)+weights_outer_product(1,t).x[3]*dX(la));
            dF(ja)+=-constant_youngs_modulus*(weights_outer_product(1,t).x[4]*dX(ia)+weights_outer_product(1,t).x[5]*dX(ja)+weights_outer_product(1,t).x[6]*dX(ka)+weights_outer_product(1,t).x[7]*dX(la));
            dF(ka)+=-constant_youngs_modulus*(weights_outer_product(1,t).x[8]*dX(ia)+weights_outer_product(1,t).x[9]*dX(ja)+weights_outer_product(1,t).x[10]*dX(ka)+weights_outer_product(1,t).x[11]*dX(la));
            dF(la)+=-constant_youngs_modulus*(weights_outer_product(1,t).x[12]*dX(ia)+weights_outer_product(1,t).x[13]*dX(ja)+weights_outer_product(1,t).x[14]*dX(ka)+weights_outer_product(1,t).x[15]*dX(la));
            //3
            dF(ia)+=constant_youngs_modulus*(weights_outer_product(3,t).x[0]*dX(ib)+weights_outer_product(3,t).x[1]*dX(jb)+weights_outer_product(3,t).x[2]*dX(kb)+weights_outer_product(3,t).x[3]*dX(lb));
            dF(ja)+=constant_youngs_modulus*(weights_outer_product(3,t).x[4]*dX(ib)+weights_outer_product(3,t).x[5]*dX(jb)+weights_outer_product(3,t).x[6]*dX(kb)+weights_outer_product(3,t).x[7]*dX(lb));
            dF(ka)+=constant_youngs_modulus*(weights_outer_product(3,t).x[8]*dX(ib)+weights_outer_product(3,t).x[9]*dX(jb)+weights_outer_product(3,t).x[10]*dX(kb)+weights_outer_product(3,t).x[11]*dX(lb));
            dF(la)+=constant_youngs_modulus*(weights_outer_product(3,t).x[12]*dX(ib)+weights_outer_product(3,t).x[13]*dX(jb)+weights_outer_product(3,t).x[14]*dX(kb)+weights_outer_product(3,t).x[15]*dX(lb));
            //2
            dF(ib)+=constant_youngs_modulus*(weights_outer_product(2,t).x[0]*dX(ia)+weights_outer_product(2,t).x[1]*dX(ja)+weights_outer_product(2,t).x[2]*dX(ka)+weights_outer_product(2,t).x[3]*dX(la));
            dF(jb)+=constant_youngs_modulus*(weights_outer_product(2,t).x[4]*dX(ia)+weights_outer_product(2,t).x[5]*dX(ja)+weights_outer_product(2,t).x[6]*dX(ka)+weights_outer_product(2,t).x[7]*dX(la));
            dF(kb)+=constant_youngs_modulus*(weights_outer_product(2,t).x[8]*dX(ia)+weights_outer_product(2,t).x[9]*dX(ja)+weights_outer_product(2,t).x[10]*dX(ka)+weights_outer_product(2,t).x[11]*dX(la));
            dF(lb)+=constant_youngs_modulus*(weights_outer_product(2,t).x[12]*dX(ia)+weights_outer_product(2,t).x[13]*dX(ja)+weights_outer_product(2,t).x[14]*dX(ka)+weights_outer_product(2,t).x[15]*dX(la));
            //4
            dF(ib)+=-constant_youngs_modulus*(weights_outer_product(4,t).x[0]*dX(ib)+weights_outer_product(4,t).x[1]*dX(jb)+weights_outer_product(4,t).x[2]*dX(kb)+weights_outer_product(4,t).x[3]*dX(lb));
            dF(jb)+=-constant_youngs_modulus*(weights_outer_product(4,t).x[4]*dX(ib)+weights_outer_product(4,t).x[5]*dX(jb)+weights_outer_product(4,t).x[6]*dX(kb)+weights_outer_product(4,t).x[7]*dX(lb));
            dF(kb)+=-constant_youngs_modulus*(weights_outer_product(4,t).x[8]*dX(ib)+weights_outer_product(4,t).x[9]*dX(jb)+weights_outer_product(4,t).x[10]*dX(kb)+weights_outer_product(4,t).x[11]*dX(lb));
            dF(lb)+=-constant_youngs_modulus*(weights_outer_product(4,t).x[12]*dX(ib)+weights_outer_product(4,t).x[13]*dX(jb)+weights_outer_product(4,t).x[14]*dX(kb)+weights_outer_product(4,t).x[15]*dX(lb));}}
    else{
        for(int t=1;t<=constrained_tets.m;t++){
            int ia,ja,ka,la;tetrahedron_mesh.elements.Get(constrained_tets(1,t),ia,ja,ka,la);
            int ib,jb,kb,lb;tetrahedron_mesh.elements.Get(constrained_tets(2,t),ib,jb,kb,lb);
            //1
            dF(ia)+=-youngs_modulus(t)*(weights_outer_product(1,t).x[0]*dX(ia)+weights_outer_product(1,t).x[1]*dX(ja)+weights_outer_product(1,t).x[2]*dX(ka)+weights_outer_product(1,t).x[3]*dX(la));
            dF(ja)+=-youngs_modulus(t)*(weights_outer_product(1,t).x[4]*dX(ia)+weights_outer_product(1,t).x[5]*dX(ja)+weights_outer_product(1,t).x[6]*dX(ka)+weights_outer_product(1,t).x[7]*dX(la));
            dF(ka)+=-youngs_modulus(t)*(weights_outer_product(1,t).x[8]*dX(ia)+weights_outer_product(1,t).x[9]*dX(ja)+weights_outer_product(1,t).x[10]*dX(ka)+weights_outer_product(1,t).x[11]*dX(la));
            dF(la)+=-youngs_modulus(t)*(weights_outer_product(1,t).x[12]*dX(ia)+weights_outer_product(1,t).x[13]*dX(ja)+weights_outer_product(1,t).x[14]*dX(ka)+weights_outer_product(1,t).x[15]*dX(la));
            //3
            dF(ia)+=youngs_modulus(t)*(weights_outer_product(3,t).x[0]*dX(ib)+weights_outer_product(3,t).x[1]*dX(jb)+weights_outer_product(3,t).x[2]*dX(kb)+weights_outer_product(3,t).x[3]*dX(lb));
            dF(ja)+=youngs_modulus(t)*(weights_outer_product(3,t).x[4]*dX(ib)+weights_outer_product(3,t).x[5]*dX(jb)+weights_outer_product(3,t).x[6]*dX(kb)+weights_outer_product(3,t).x[7]*dX(lb));
            dF(ka)+=youngs_modulus(t)*(weights_outer_product(3,t).x[8]*dX(ib)+weights_outer_product(3,t).x[9]*dX(jb)+weights_outer_product(3,t).x[10]*dX(kb)+weights_outer_product(3,t).x[11]*dX(lb));
            dF(la)+=youngs_modulus(t)*(weights_outer_product(3,t).x[12]*dX(ib)+weights_outer_product(3,t).x[13]*dX(jb)+weights_outer_product(3,t).x[14]*dX(kb)+weights_outer_product(3,t).x[15]*dX(lb));
            //2
            dF(ib)+=youngs_modulus(t)*(weights_outer_product(2,t).x[0]*dX(ia)+weights_outer_product(2,t).x[1]*dX(ja)+weights_outer_product(2,t).x[2]*dX(ka)+weights_outer_product(2,t).x[3]*dX(la));
            dF(jb)+=youngs_modulus(t)*(weights_outer_product(2,t).x[4]*dX(ia)+weights_outer_product(2,t).x[5]*dX(ja)+weights_outer_product(2,t).x[6]*dX(ka)+weights_outer_product(2,t).x[7]*dX(la));
            dF(kb)+=youngs_modulus(t)*(weights_outer_product(2,t).x[8]*dX(ia)+weights_outer_product(2,t).x[9]*dX(ja)+weights_outer_product(2,t).x[10]*dX(ka)+weights_outer_product(2,t).x[11]*dX(la));
            dF(lb)+=youngs_modulus(t)*(weights_outer_product(2,t).x[12]*dX(ia)+weights_outer_product(2,t).x[13]*dX(ja)+weights_outer_product(2,t).x[14]*dX(ka)+weights_outer_product(2,t).x[15]*dX(la));
            //4
            dF(ib)+=-youngs_modulus(t)*(weights_outer_product(4,t).x[0]*dX(ib)+weights_outer_product(4,t).x[1]*dX(jb)+weights_outer_product(4,t).x[2]*dX(kb)+weights_outer_product(4,t).x[3]*dX(lb));
            dF(jb)+=-youngs_modulus(t)*(weights_outer_product(4,t).x[4]*dX(ib)+weights_outer_product(4,t).x[5]*dX(jb)+weights_outer_product(4,t).x[6]*dX(kb)+weights_outer_product(4,t).x[7]*dX(lb));
            dF(kb)+=-youngs_modulus(t)*(weights_outer_product(4,t).x[8]*dX(ib)+weights_outer_product(4,t).x[9]*dX(jb)+weights_outer_product(4,t).x[10]*dX(kb)+weights_outer_product(4,t).x[11]*dX(lb));
            dF(lb)+=-youngs_modulus(t)*(weights_outer_product(4,t).x[12]*dX(ib)+weights_outer_product(4,t).x[13]*dX(jb)+weights_outer_product(4,t).x[14]*dX(kb)+weights_outer_product(4,t).x[15]*dX(lb));}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>::
Add_Velocity_Dependent_Forces(const ARRAY<VECTOR<T,3> >& V,ARRAY<VECTOR<T,3> >& F,const T time) const
{

}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>::
Initialize_CFL()
{
    CFL_elastic_time_step_of_fragment.Append(FLT_MAX);CFL_damping_time_step_of_fragment.Append(FLT_MAX);
    CFL_initialized_of_fragment.Append(true);
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class T> void BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<T>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
template class BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING<double>;
#endif
