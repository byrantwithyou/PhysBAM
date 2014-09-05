//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Deformables/Forces/BEZIER_CURVATURE_FORCE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BEZIER_CURVATURE_FORCE<TV>::
BEZIER_CURVATURE_FORCE(DEFORMABLE_PARTICLES<TV>& particles,const BEZIER_SPLINE<TV,3>& spline_input,
    T curvature_stiffness_input,T length_stiffness_input)
    :DEFORMABLES_FORCES<TV>(particles),spline(spline_input),curvature_stiffness(curvature_stiffness_input),
    length_stiffness(length_stiffness_input)
{
    X0=spline.particles.X;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BEZIER_CURVATURE_FORCE<TV>::
~BEZIER_CURVATURE_FORCE()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void BEZIER_CURVATURE_FORCE<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=0;i<spline.control_points.m;i++)
        F.Subset(spline.control_points(i))-=data(i).ge;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void BEZIER_CURVATURE_FORCE<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void BEZIER_CURVATURE_FORCE<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const
{
    for(int i=0;i<spline.control_points.m;i++){
        const VECTOR<int,4>& nodes=spline.control_points(i);
        const MATRIX<MATRIX<T,TV::m>,4>& he=data(i).he;
        VECTOR<TV,4> v(V.Subset(nodes)*scale),f;
        for(int j=0;j<4;j++)
            for(int k=0;k<4;k++)
                f(j)+=he(j,k)*v(k);
        F.Subset(nodes)-=f;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void BEZIER_CURVATURE_FORCE<TV>::
Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR BEZIER_CURVATURE_FORCE<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void BEZIER_CURVATURE_FORCE<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int i=0;i<spline.control_points.m;i++){
        const VECTOR<int,4>& nodes=spline.control_points(i);
        for(int j=1;j<4;j++)
            dependency_mesh.elements.Append(VECTOR<int,2>(nodes(0),nodes(j)));}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void BEZIER_CURVATURE_FORCE<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    MATRIX<T,2> rot(0,1,-1,0);
    pe=0;

    if(data.m!=spline.control_points.m){
        data.Resize(spline.control_points.m);
        for(int i=0;i<spline.control_points.m;i++){
            const VECTOR<int,4>& nodes=spline.control_points(i);
            VECTOR<TV,4> XX(X0.Subset(nodes));
            TV phi=(XX(0)+XX(1)*3+XX(2)*3+XX(3))*(T).125;
            TV phi_p=(-XX(0)-XX(1)+XX(2)+XX(3))*(T).75;
            TV phi_pp=(XX(0)-XX(1)-XX(2)+XX(3))*3;
            T length=phi_p.Magnitude();
            T m_bar=1/length;
            T m_hat=sqr(m_bar);
            TV TT=m_bar*phi_p;
            TV N=rot*TT;
            T xi=N.Dot(phi_pp);
            T kappa=m_hat*xi;
            data(i).kappa0=kappa;
            data(i).length0=length;}}

    for(int i=0;i<spline.control_points.m;i++){
        const VECTOR<int,4>& nodes=spline.control_points(i);
        DATA& dat=data(i);
        VECTOR<TV,4> XX(particles.X.Subset(nodes));
        auto Z0=From_Var<4,0>(XX(0));
        auto Z1=From_Var<4,1>(XX(1));
        auto Z2=From_Var<4,2>(XX(2));
        auto Z3=From_Var<4,3>(XX(3));
        auto phi=(Z0+3*Z1+3*Z2+Z3)*(T).125;
        auto phi_p=(-Z0-Z1+Z2+Z3)*(T).75;
        auto phi_pp=(Z0-Z1-Z2+Z3)*3;
        auto length=phi_p.Magnitude();
        auto m_bar=1/length;
        auto m_hat=sqr(m_bar);
        auto TT=m_bar*phi_p;
        auto N=rot*TT;
        auto xi=N.Dot(phi_pp);
        auto kappa=m_hat*xi;
        auto tau=kappa-dat.kappa0;
        auto curvature_pe=curvature_stiffness/2*sqr(tau);
        auto dl=length-dat.length0;
        auto length_pe=length_stiffness/2*sqr(dl);
        auto new_pe=curvature_pe+length_pe;
        pe+=new_pe.x;
        Extract(dat.ge,new_pe.dx);
        Extract(dat.he,new_pe.ddx);}
    printf("PE: %g\n",pe);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BEZIER_CURVATURE_FORCE<TV>::
Potential_Energy(const T time) const
{
    return pe;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void BEZIER_CURVATURE_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
template class BEZIER_CURVATURE_FORCE<VECTOR<float,2> >;
template class BEZIER_CURVATURE_FORCE<VECTOR<double,2> >;
}
