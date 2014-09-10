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
namespace{
double quadrature_loc[8][7]=
{
    {},
    {0.5},
    {0.211324865405187118,1-0.211324865405187118},
    {0.5,0.112701665379258311,1-0.112701665379258311},
    {0.069431844202973713,0.930568155797026287,0.669990521792428133,0.330009478207571867},
    {0.500000000000000000,0.230765344947158455,0.769234655052841545,0.953089922969331997,0.046910077030668003},
    {0.966234757101576014,0.033765242898423986,0.619309593041598456,0.380690406958401544,0.169395306766867742,0.830604693233132258},
    {0.500000000000000000,0.974553956171379264,0.0254460438286207357,0.702922575688698588,0.297077424311301413,0.870765592799697222,0.129234407200302778}
};
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
            for(int j=0;j<gauss_order;j++){
                T t=quadrature_loc[gauss_order][j];
                TV phi_p=-3*sqr(t-1)*XX(0)+3*(t-1)*(3*t-1)*XX(1)-3*t*(3*t-2)*XX(2)+3*sqr(t)*XX(3);
                TV phi_pp=-6*(t-1)*XX(0)+6*(3*t-2)*XX(1)-6*(3*t-1)*XX(2)+6*t*XX(3);
                T length=phi_p.Magnitude();
                T m_bar=1/length;
                T m_hat=cube(m_bar);
                TV N=rot*phi_p;
                T xi=N.Dot(phi_pp);
                T kappa=m_hat*xi;
                data(i).kappa0(j)=kappa;
                data(i).length0(j)=length;}}}

    for(int i=0;i<spline.control_points.m;i++){
        const VECTOR<int,4>& nodes=spline.control_points(i);
        DATA& dat=data(i);
        VECTOR<TV,4> XX(particles.X.Subset(nodes));
        VECTOR<TV,4> ge;
        MATRIX<MATRIX<T,TV::m>,4> he;
        dat.ge=ge;
        dat.he=he;
        for(int j=0;j<gauss_order;j++){
            T t=quadrature_loc[gauss_order][j];
            auto Z0=From_Var<4,0>(XX(0));
            auto Z1=From_Var<4,1>(XX(1));
            auto Z2=From_Var<4,2>(XX(2));
            auto Z3=From_Var<4,3>(XX(3));
            auto phi_p=-3*sqr(t-1)*Z0+3*(t-1)*(3*t-1)*Z1-3*t*(3*t-2)*Z2+3*sqr(t)*Z3;
            auto phi_pp=-6*(t-1)*Z0+6*(3*t-2)*Z1-6*(3*t-1)*Z2+6*t*Z3;
            auto length=phi_p.Magnitude();
            auto m_bar=1/length;
            auto m_hat=cube(m_bar);
            auto N=rot*phi_p;
            auto xi=N.Dot(phi_pp);
            auto kappa=m_hat*xi;
            auto tau=kappa-dat.kappa0(j);
            auto curvature_pe=curvature_stiffness/2*sqr(tau);
            auto dl=length/dat.length0(j)-1;
            auto length_pe=length_stiffness/2*sqr(dl);
            auto new_pe=quadrature_weight[gauss_order][j]/spline.control_points.m*(curvature_pe+length_pe);
            pe+=new_pe.x;
            Extract(ge,new_pe.dx);
            Extract(he,new_pe.ddx);
            dat.ge+=ge;
            dat.he+=he;}}
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
