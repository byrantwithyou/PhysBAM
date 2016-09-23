//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Deformables/Forces/OPENSUBDIV_SURFACE_CURVATURE_FORCE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{
namespace{
double quadrature_weight[8][7]=
{
    {},
    {1},
    {0.5,0.5},
    {0.444444444444444444,0.277777777777777778,0.277777777777777778},
    {0.173927422568726929,0.173927422568726929,0.326072577431273071,0.326072577431273071},
    {0.284444444444444444,0.239314335249683234,0.239314335249683234,0.118463442528094544,0.118463442528094544},
    {0.0856622461895851732,0.0856622461895851732,0.233956967286345525,0.233956967286345525,0.180380786524069304,0.180380786524069304},
    {0.208979591836734694,0.0647424830844348507,0.0647424830844348507,0.190915025252559459,0.190915025252559459,0.139852695744638342,0.139852695744638342}
};
}
//#####################################################################
//Constructor
//#####################################################################
template<class T,int gauss_order> OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
OPENSUBDIV_SURFACE_CURVATURE_FORCE(DEFORMABLE_PARTICLES<TV>& particles,const OPENSUBDIV_SURFACE<TV>& surf_input,const MOONEY_RIVLIN_CURVATURE<T>& model_input)
    :DEFORMABLES_FORCES<TV>(particles),surf(surf_input),model(model_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int gauss_order> OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
~OPENSUBDIV_SURFACE_CURVATURE_FORCE()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T,int gauss_order> void OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int e=0;e<surf.m;e++){
        const ARRAY<int> a(surf.control_points.Subset(surf.face_data(e).nodes));
        for(int i=0;i<a.m;i++)
            F(a(i))-=data(e).ge(i);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T,int gauss_order> void OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class T,int gauss_order> void OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
#pragma omp parallel for 
    for(int face=0;face<surf.m;face++){
        const ARRAY<int> nodes(surf.control_points.Subset(surf.face_data(face).nodes));
        const ARRAY<MATRIX<VECTOR<T,5>,gauss_order> >& A=surf.face_data(face).A;
        ARRAY<TV> f(nodes.m);
        ARRAY<TV> v(V.Subset(nodes));
        MATRIX<VECTOR<TV,5>,gauss_order> Av;
        MATRIX<VECTOR<TV,5>,gauss_order> heAv;
        for(int ga=0;ga<nodes.m;ga++){
            const TV& vga=v(ga);
            for(int i=0;i<gauss_order;i++)
                for(int j=0;j<gauss_order;j++){
                    VECTOR<TV,5>& Avij=Av(i,j);
                    const VECTOR<T,5>& Aga=A(ga)(i,j);
                    for(int five=0;five<5;five++)
                        Avij(five)+=Aga(five)*vga;}}
        
        for(int i=0;i<gauss_order;i++)
            for(int j=0;j<gauss_order;j++){
                const VECTOR<TV,5>& Avij=Av(i,j);
                VECTOR<TV,5>& heAvij=heAv(i,j);
                const TT& he=data(face).he(i,j);
                for(int five2=0;five2<5;five2++)
                    for(int five1=0;five1<5;five1++)
                        heAvij(five1)+=he(five1,five2)*Avij(five2);}

        for(int al=0;al<nodes.m;al++){
            TV& fal=f(al);
            for(int i=0;i<gauss_order;i++)
                for(int j=0;j<gauss_order;j++){
                    const VECTOR<TV,5>& heAvij=heAv(i,j);
                    const VECTOR<T,5>& Aal=A(al)(i,j);
                    for(int five1=0;five1<5;five1++)
                        fal+=Aal(five1)*heAvij(five1);}}
#pragma omp critical
                F.Subset(nodes)-=f;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class T,int gauss_order> void OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T,int gauss_order> T OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T,int gauss_order> void OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int i=1;i<surf.control_points.m;i++)
        dependency_mesh.elements.Append(VECTOR<int,2>(surf.control_points(0),surf.control_points(i)));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int gauss_order> void OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    T local_pe=0;
    if(data.m!=surf.m){
        data.Resize(surf.m);
        for(int face=0;face<surf.m;face++){
            data(face).ge.Resize(surf.face_data(face).nodes.m);}}

#pragma omp parallel for reduction(+:local_pe)
    for(int face=0;face<surf.m;face++){
        TM2 ge;
        const ARRAY<int> nodes(surf.control_points.Subset(surf.face_data(face).nodes));
        const ARRAY<MATRIX<VECTOR<T,5>,gauss_order> >& A=surf.face_data(face).A;
        for(int i=0;i<nodes.m;i++) data(face).ge(i)=TV();

        for(int i=0;i<gauss_order;i++){
            for(int j=0;j<gauss_order;j++){
                VECTOR<TV,5> a;
                for(int al=0;al<nodes.m;al++)
                    for(int five=0;five<5;five++)
                        a(five)+=A(al)(i,j)(five)*particles.X(nodes(al));

                T weight=quadrature_weight[gauss_order][i]*quadrature_weight[gauss_order][j];
                if(update_hessian)
                    local_pe+=model.Potential_Energy(a(0),a(1),a(2),a(3),a(4),surf.face_data(face).G0_inv(i,j),surf.face_data(face).G0_det(i,j),ge,data(face).he(i,j),weight);
                else
                    local_pe+=model.Potential_Energy(a(0),a(1),a(2),a(3),a(4),surf.face_data(face).G0_inv(i,j),surf.face_data(face).G0_det(i,j),ge,weight);

                for(int al=0;al<nodes.m;al++)
                    for(int five=0;five<5;five++)
                        data(face).ge(al)+=ge(five)*A(al)(i,j)(five);}}}
    pe=local_pe;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T,int gauss_order> T OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
Potential_Energy(const T time) const
{
    return pe;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T,int gauss_order> void OPENSUBDIV_SURFACE_CURVATURE_FORCE<T,gauss_order>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<float,1>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<double,1>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<float,2>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<double,2>;
template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<float,3>;
template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<double,3>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<float,4>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<double,4>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<float,5>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<double,5>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<float,6>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<double,6>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<float,7>;
//template class OPENSUBDIV_SURFACE_CURVATURE_FORCE<double,7>;
}
