//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Deformables/Forces/B_SPLINE_CURVATURE_FORCE.h>
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
template<class TV> B_SPLINE_CURVATURE_FORCE<TV>::
B_SPLINE_CURVATURE_FORCE(DEFORMABLE_PARTICLES<TV>& particles,const B_SPLINE<TV,3>& spline_input,
    T curvature_stiffness_input,T length_stiffness_input)
    :DEFORMABLES_FORCES<TV>(particles),spline(spline_input),curvature_stiffness(curvature_stiffness_input),
    length_stiffness(length_stiffness_input)
{
    X0=spline.particles.X;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> B_SPLINE_CURVATURE_FORCE<TV>::
~B_SPLINE_CURVATURE_FORCE()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void B_SPLINE_CURVATURE_FORCE<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=2;i<spline.knots.m-3;i++)
        for(int j=0;j<4;j++)
            F(spline.control_points(i-2+j))-=data(i-2).ge(j);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void B_SPLINE_CURVATURE_FORCE<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void B_SPLINE_CURVATURE_FORCE<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=2;i<spline.knots.m-3;i++){
        VECTOR<int,4> nodes(spline.control_points.Subset(VECTOR<int,4>(i-2,i-1,i,i+1)));
        const MATRIX<MATRIX<T,TV::m>,4>& he=data(i-2).he;
        VECTOR<TV,4> v(V.Subset(nodes)),f;
        for(int j=0;j<4;j++)
            for(int k=0;k<4;k++)
                f(j)+=he(j,k)*v(k);
        F.Subset(nodes)-=f;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void B_SPLINE_CURVATURE_FORCE<TV>::
Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR B_SPLINE_CURVATURE_FORCE<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void B_SPLINE_CURVATURE_FORCE<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int i=1; i<spline.control_points.m;i++)
        dependency_mesh.elements.Append(VECTOR<int,2>(spline.control_points(0),spline.control_points(i)));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void B_SPLINE_CURVATURE_FORCE<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    pe=0;

    MATRIX<T,2> rot(0,1,-1,0);
    VECTOR<T,6> k;
    VECTOR<TV,4> XX;

    if(data.m!=spline.knots.m-5){
        data.Resize(spline.knots.m-5);
        for(int i=2;i<spline.knots.m-3;i++){
            for(int ctr=0;ctr<4;ctr++) XX(ctr)=X0(spline.control_points(i-2+ctr));
            for(int ctr=0;ctr<6;ctr++) k[ctr]=spline.knots(i-2+ctr);
            for(int j=0;j<gauss_order;j++){
                T t=quadrature_loc[gauss_order][j];
                t=k[2]+(k[3]-k[2])*t;

                TV x10=XX(0)+(XX(1)-XX(0))*(t-k[0])/(k[3]-k[0]);
                TV x11=XX(1)+(XX(2)-XX(1))*(t-k[1])/(k[4]-k[1]);
                TV x12=XX(2)+(XX(3)-XX(2))*(t-k[2])/(k[5]-k[2]);
                TV x20=x10+(x11-x10)*(t-k[1])/(k[3]-k[1]);
                TV x21=x11+(x12-x11)*(t-k[2])/(k[4]-k[2]);
                
                TV x10_p=(XX(1)-XX(0))/(k[3]-k[0]);
                TV x11_p=(XX(2)-XX(1))/(k[4]-k[1]);
                TV x12_p=(XX(3)-XX(2))/(k[5]-k[2]);
                TV x20_p=x10_p+(x11_p-x10_p)*(t-k[1])/(k[3]-k[1])+(x11-x10)/(k[3]-k[1]);
                TV x21_p=x11_p+(x12_p-x11_p)*(t-k[2])/(k[4]-k[2])+(x12-x11)/(k[4]-k[2]);
                TV phi_p=x20_p+(x21_p-x20_p)*(t-k[2])/(k[3]-k[2])+(x21-x20)/(k[3]-k[2]);
                
                TV x20_pp=(x11_p-x10_p)/(k[3]-k[1]);
                TV x21_pp=(x12_p-x11_p)/(k[4]-k[2]);
                TV phi_pp=x20_pp+(x21_pp-x20_pp)*(t-k[2])/(k[3]-k[2])+(x21_p-x20_p)/(k[3]-k[2])+(x21_p-x20_p)/(k[3]-k[2]);

                T length=phi_p.Magnitude();
                T m_bar=1/length;
                T m_hat=cube(m_bar);
                TV N=rot*phi_p;
                T xi=N.Dot(phi_pp);
                T kappa=m_hat*xi;
                data(i-2).kappa0(j)=kappa;
                data(i-2).length0(j)=length;}}}

    for(int i=2;i<spline.knots.m-3;i++){
        DATA& dat=data(i-2);
        for(int ctr=0;ctr<4;ctr++) XX(ctr)=particles.X(spline.control_points(i-2+ctr));
        for(int ctr=0;ctr<6;ctr++) k[ctr]=spline.knots(i-2+ctr);
        VECTOR<TV,4> ge;
        MATRIX<MATRIX<T,TV::m>,4> he;
        dat.ge=ge;
        dat.he=he;
        for(int j=0;j<gauss_order;j++){
            T t=quadrature_loc[gauss_order][j];
            t=k[2]+(k[3]-k[2])*t;
            
            typedef DIFF_LAYOUT<T,TV::m,TV::m,TV::m,TV::m> LAYOUT;
            auto x00=Hess_From_Var<LAYOUT,0>(XX(0));
            auto x01=Hess_From_Var<LAYOUT,1>(XX(1));
            auto x02=Hess_From_Var<LAYOUT,2>(XX(2));
            auto x03=Hess_From_Var<LAYOUT,3>(XX(3));

            auto x10=x00+(x01-x00)*(t-k[0])/(k[3]-k[0]);
            auto x11=x01+(x02-x01)*(t-k[1])/(k[4]-k[1]);
            auto x12=x02+(x03-x02)*(t-k[2])/(k[5]-k[2]);
            auto x20=x10+(x11-x10)*(t-k[1])/(k[3]-k[1]);
            auto x21=x11+(x12-x11)*(t-k[2])/(k[4]-k[2]);
            
            auto x10_p=(x01-x00)/(k[3]-k[0]);
            auto x11_p=(x02-x01)/(k[4]-k[1]);
            auto x12_p=(x03-x02)/(k[5]-k[2]);
            auto x20_p=x10_p+(x11_p-x10_p)*(t-k[1])/(k[3]-k[1])+(x11-x10)/(k[3]-k[1]);
            auto x21_p=x11_p+(x12_p-x11_p)*(t-k[2])/(k[4]-k[2])+(x12-x11)/(k[4]-k[2]);
            auto phi_p=x20_p+(x21_p-x20_p)*(t-k[2])/(k[3]-k[2])+(x21-x20)/(k[3]-k[2]);
            
            auto x20_pp=(x11_p-x10_p)/(k[3]-k[1]);
            auto x21_pp=(x12_p-x11_p)/(k[4]-k[2]);
            auto phi_pp=x20_pp+(x21_pp-x20_pp)*(t-k[2])/(k[3]-k[2])+(x21_p-x20_p)/(k[3]-k[2])+(x21_p-x20_p)/(k[3]-k[2]);

            auto length=phi_p.Magnitude();
            auto m_bar=(T)1/length;
            auto m_hat=cube(m_bar);
            auto N=rot*phi_p;
            auto xi=N.Dot(phi_pp);
            auto kappa=m_hat*xi;
//            auto tau=kappa-dat.kappa0(j);
            auto tau=atan2(kappa,(T)1)-atan2(dat.kappa0(j),(T)1);
            auto curvature_pe=curvature_stiffness/2*sqr(tau);
            auto dl=length/dat.length0(j)-1;
            auto length_pe=length_stiffness/2*sqr(dl);
            auto new_pe=(T)quadrature_weight[gauss_order][j]/spline.control_points.m*(curvature_pe+length_pe);
            pe+=new_pe.x;
            Extract(ge,new_pe.dx);
            Extract(he,new_pe.ddx);
            dat.ge+=ge;
            dat.he+=he;}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR B_SPLINE_CURVATURE_FORCE<TV>::
Potential_Energy(const T time) const
{
    return pe;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void B_SPLINE_CURVATURE_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
template class B_SPLINE_CURVATURE_FORCE<VECTOR<float,2> >;
template class B_SPLINE_CURVATURE_FORCE<VECTOR<double,2> >;
}
