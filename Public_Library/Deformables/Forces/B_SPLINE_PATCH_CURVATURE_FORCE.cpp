//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Deformables/Forces/B_SPLINE_PATCH_CURVATURE_FORCE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
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
//Constructor
//#####################################################################
template<class T,int gauss_order> B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
B_SPLINE_PATCH_CURVATURE_FORCE(DEFORMABLE_PARTICLES<TV>& particles,const B_SPLINE_PATCH<TV,3>& spline_input,const MOONEY_RIVLIN_CURVATURE<T>& model_input,T density)
    :LAZY_HESSIAN_FORCE<TV>(particles),spline(spline_input),model(model_input),density(density),recompute_hessian(true)
{
    X0=spline.particles.X;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int gauss_order> B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
~B_SPLINE_PATCH_CURVATURE_FORCE()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T,int gauss_order> void B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int e=0;e<spline.m;e++){
        VECTOR<int,16> a = spline.Control_Points_For_Element(e);
        for(int i=0;i<a.m;i++)
            F(a(i))-=data(e).ge.Column(i);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T,int gauss_order> void B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class T,int gauss_order> void B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const
{
#pragma omp parallel for 
    for(int e=0;e<spline.m;e++){
        VECTOR<int,16> nodes = spline.Control_Points_For_Element(e);
        RANGE<VECTOR<T,2>> r=spline.Range_For_Element(e);
        const DATA& dat=data(e);
        VECTOR<TV,16> v(V.Subset(nodes)*scale),f;
        for(int i=0;i<gauss_order;i++)
            for(int j=0;j<gauss_order;j++){
                const VECTOR<VECTOR<T,16>,5>& Aij=dat.A(i,j);
                const TT &he=dat.he(i,j);
                T weight=quadrature_weight[gauss_order][i]*quadrature_weight[gauss_order][j]*r.Size();
                VECTOR<TV,5> Av;
                for(int five=0;five<5;five++)
                    for(int ga=0;ga<nodes.m;ga++)
                        Av(five)+=weight*Aij(five)(ga)*v(ga);

                VECTOR<TV,5> heAv;
                for(int five1=0;five1<5;five1++)
                    for(int five2=0;five2<5;five2++)
                        heAv(five1)+=he(five1,five2)*Av(five2);

                for(int al=0;al<nodes.m;al++)
                    for(int five1=0;five1<5;five1++)
                        f(al)+=Aij(five1)(al)*heAv(five1);}
#pragma omp critical
        F.Subset(nodes)-=f;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class T,int gauss_order> void B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T,int gauss_order> T B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T,int gauss_order> void B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(RANGE_ITERATOR<IV::m> it(spline.control_points.domain);it.Valid();it.Next())
        dependency_mesh.elements.Append(VECTOR<int,2>(spline.control_points(0,0),spline.control_points(it.index)));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T,int gauss_order> void B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    T local_pe=0;
    if(data.m!=spline.m){
        VECTOR<T,16> w;
        VECTOR<VECTOR<T,16>,2> dw;
        VECTOR<VECTOR<T,16>,3> ddw;
        data.Resize(spline.m);
        for(int e=0;e<spline.m;e++){
            VECTOR<int,16> nodes=spline.Control_Points_For_Element(e);
            VECTOR<TV,16> x0(X0.Subset(nodes));
            RANGE<VECTOR<T,2>> r=spline.Range_For_Element(e);
            VECTOR<T,2> dr=r.Edge_Lengths();
            DATA& dat=data(e);
            for(int i=0;i<gauss_order;i++){
                T s=quadrature_loc[gauss_order][i];
                s=r.min_corner(0)+dr(0)*s;
                for(int j=0;j<gauss_order;j++){
                    T t=quadrature_loc[gauss_order][j];
                    t=r.min_corner(1)+dr(1)*t;
                    VECTOR<VECTOR<T,16>,5>& Aij=dat.A(i,j);
                    VECTOR<TM,3>& G0_inv=dat.G0_inv(i,j);
                    VECTOR<T,3>& G0_det=dat.G0_det(i,j);
                    spline.Calculate_Weights(s,t,w,Aij(0),Aij(1),Aij(2),Aij(3),Aij(4));
                    w*=w.Dot(VECTOR<T,16>::All_Ones_Vector());

                    VECTOR<TV,5> a;
                    for(int five=0;five<5;five++)
                        for(int al=0;al<nodes.m;al++)
                            a(five)+=Aij(five)(al)*x0(al);

                    TV a3=a(0).Cross(a(1));
                    TV da31=a(0).Cross(a(3))+a(2).Cross(a(1));
                    TV da32=a(0).Cross(a(4))+a(3).Cross(a(1));

                    T m_sqr=a3.Magnitude_Squared();
                    T dm_sqr1=2*da31.Dot(a3);
                    T dm_sqr2=2*da32.Dot(a3);

                    T J0=sqrt(m_sqr);
                    T beta=J0/m_sqr;

                    TV lambda_n=beta*a3;
                    TV dlambda_n1=beta*(da31-(dm_sqr1/m_sqr)*a3);
                    TV dlambda_n2=beta*(da31-(dm_sqr2/m_sqr)*a3);

                    T mass=0;
                    const VECTOR<T,3> simpson_weights=TV(1,4,1);
                    for(int k=0;k<3;k++){
                        G0_inv(k)=TM(a(0)+(k-1)*model.thickness/2*dlambda_n1,a(1)+(k-1)*model.thickness/2*dlambda_n2,lambda_n);
                        if(G0_inv(k).Determinant()<0) LOG::cout<<"det="<<G0_inv(k).Determinant()<<std::endl;
                        G0_det(k)=abs(G0_inv(k).Determinant());
                        mass+=simpson_weights(k)*G0_det(k);
                        G0_inv(k).Invert();}

                    mass*=density*model.thickness*quadrature_weight[gauss_order][i]*quadrature_weight[gauss_order][j]*r.Size()/6;
                    PHYSBAM_ASSERT(mass>0);
                    for(int a=0;a<nodes.m;a++)
                        particles.mass(nodes(a))+=mass*w(a);}}}}

#pragma omp parallel for reduction(+:local_pe)
    for(int e=0;e<spline.m;e++){
        TM2 ge;
        VECTOR<int,16> nodes=spline.Control_Points_For_Element(e);
        VECTOR<TV,16> x(particles.X.Subset(nodes));
        RANGE<VECTOR<T,2>> r=spline.Range_For_Element(e);
        DATA& dat=data(e);
        dat.ge=MATRIX<T,3,16>();
        for(int i=0;i<gauss_order;i++){
            for(int j=0;j<gauss_order;j++){
                VECTOR<VECTOR<T,16>,5>& Aij=dat.A(i,j);

                VECTOR<TV,5> a;
                for(int five=0;five<5;five++)
                    for(int al=0;al<nodes.m;al++)
                        a(five)+=Aij(five)(al)*x(al);

                T weight=quadrature_weight[gauss_order][i]*quadrature_weight[gauss_order][j]*r.Size();
                if(recompute_hessian)
                    local_pe+=weight*model.Potential_Energy(a(0),a(1),a(2),a(3),a(4),dat.G0_inv(i,j),dat.G0_det(i,j),ge,dat.he(i,j));
                else
                    local_pe+=weight*model.Potential_Energy(a(0),a(1),a(2),a(3),a(4),dat.G0_inv(i,j),dat.G0_det(i,j),ge);
                ge*=weight;

                for(int al=0;al<nodes.m;al++)
                    for(int dim=0;dim<3;dim++)
                        for(int five=0;five<5;five++)
                            dat.ge(dim,al)+=ge(dim,five)*Aij(five)(al);}}}
    pe=local_pe;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T,int gauss_order> T B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
Potential_Energy(const T time) const
{
    return pe;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T,int gauss_order> void B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Need_To_Recompute_Hessian
//#####################################################################
template<class T,int gauss_order> void B_SPLINE_PATCH_CURVATURE_FORCE<T,gauss_order>::
Need_To_Recompute_Hessian(bool h)
{
    recompute_hessian=h;
}
template class B_SPLINE_PATCH_CURVATURE_FORCE<float,1>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<double,1>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<float,2>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<double,2>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<float,3>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<double,3>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<float,4>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<double,4>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<float,5>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<double,5>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<float,6>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<double,6>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<float,7>;
template class B_SPLINE_PATCH_CURVATURE_FORCE<double,7>;
}
