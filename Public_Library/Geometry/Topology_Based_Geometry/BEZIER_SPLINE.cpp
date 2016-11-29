//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BEZIER_SPLINE  
//##################################################################### 
#include <Core/Arrays/IDENTITY_ARRAY.h>
#include <Core/Matrices/BANDED_MATRIX.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> BEZIER_SPLINE<TV,d>::
BEZIER_SPLINE()
    :need_destroy_particles(true),particles(*new GEOMETRY_PARTICLES<TV>)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> BEZIER_SPLINE<TV,d>::
BEZIER_SPLINE(GEOMETRY_PARTICLES<TV>& particles)
    :need_destroy_particles(false),particles(particles)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> BEZIER_SPLINE<TV,d>::
~BEZIER_SPLINE()
{
    if(need_destroy_particles) delete &particles;
}
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV,int d> BEZIER_SPLINE<TV,d>* BEZIER_SPLINE<TV,d>::
Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const
{
    BEZIER_SPLINE* bs=new BEZIER_SPLINE(new_particles);
    int offset=new_particles.Size();
    new_particles.Append(particles);
    if(particle_indices) particle_indices->Append_Elements(IDENTITY_ARRAY<>(particles.Size())+offset);
    bs->control_points=control_points+offset;
    return bs;
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV,int d> void BEZIER_SPLINE<TV,d>::
Read(TYPED_ISTREAM& input)
{
    int size;
    Read_Binary(input,control_points,size);
    particles.Clean_Memory(); // strip everything away except for position
    particles.Resize(size);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,particles.X.Get_Array_Pointer(),size);
    else Read_Binary_Array<float>(input.stream,particles.X.Get_Array_Pointer(),size);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV,int d> void BEZIER_SPLINE<TV,d>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,control_points,particles.X);
}
//#####################################################################
// Function Smooth_Fit
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit(BEZIER_SPLINE<TV,3>& bs,ARRAY_VIEW<TV> X)
{
    typedef typename TV::SCALAR T;
    ARRAY<TV> rhs(X.m*2-2);
    BANDED_MATRIX<T,4> matrix(X.m*2-2);
    matrix.diagonal_column=1;
    matrix.A(0)=VECTOR<T,4>(0,2,-1,0);
    rhs(0)=X(0);
    for(int i=1;i<X.m-1;i++){
        matrix.A(2*i-1)=VECTOR<T,4>(-1,2,-2,1);
        matrix.A(2*i)=VECTOR<T,4>(1,1,0,0);
        rhs(2*i)=X(i)*2;}
    matrix.A(X.m*2-3)=VECTOR<T,4>(-1,2,0,0);
    rhs(X.m*2-3)=X(X.m-1);
    matrix.QR_Solve(rhs);

    bs.particles.Add_Elements(X.m*3-2);
    VECTOR<int,4> base(0,1,2,3);
    for(int i=0;i<X.m-1;i++){
        bs.particles.X(3*i)=X(i);
        bs.particles.X(3*i+1)=rhs(2*i);
        bs.particles.X(3*i+2)=rhs(2*i+1);
        bs.control_points.Append(base+3*i);}
    bs.particles.X.Last()=X.Last();
}
//#####################################################################
// Function Smooth_Fit_Loop
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit_Loop(BEZIER_SPLINE<TV,3>& bs,ARRAY_VIEW<TV> X)
{
    typedef typename TV::SCALAR T;
    T a0=1,a1=0,b0=0,b1=1;
    TV z0,z1;
    for(int i=0;i<X.m;i++)
    {
        T c0=-a1,d0=-b1,c1=a0-2*a1+2*b0,d1=b0-2*b1+2*d0;
        TV y0=2*X(i)-z1,y1=z0-2*z1+2*y0;
        a0=c0;
        a1=c1;
        b0=d0;
        b1=d1;
        z0=y0;
        z1=y1;
    }
    T r=a0*b1-a0-b1+1-b0*a1;
    TV P=((1-b1)*z0+b0*z1)/r;
    TV Q=((1-a0)*z1+a1*z0)/r;
    bs.particles.Add_Elements(X.m*3-2);
    bs.particles.X(0)=X(0);
    bs.particles.X(1)=P;
    bs.particles.X(2)=Q;
    for(int i=0;i<X.m-1;i++)
    {
        TV p=-Q+2*X(i);
        TV q=P-2*Q+2*p;
        bs.particles.X(3*i)=X(i);
        bs.particles.X(3*i+1)=p;
        bs.particles.X(3*i+2)=q;
        P=p;
        Q=q;
    }
    bs.particles.X.Last()=X.Last();
}
//#####################################################################
// Function Create_Segmented_Curve
//#####################################################################
template<class TV,int d> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT* PhysBAM::
Create_Segmented_Curve(const BEZIER_SPLINE<TV,d>& spline,bool same_particles)
{
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT T_SEGMENTED_CURVE;
    T_SEGMENTED_CURVE* segmented_curve=0;
    if(same_particles) segmented_curve=T_SEGMENTED_CURVE::Create(spline.particles);
    else segmented_curve=T_SEGMENTED_CURVE::Create();
    for(int i=0;i<spline.control_points.m;i++){
        const VECTOR<int,d+1>& elem=spline.control_points(i);
        segmented_curve->mesh.elements.Append(VECTOR<int,2>(elem(0),elem.Last()));}

    if(!same_particles){
        ARRAY<int> map(spline.particles.X.m,true,-1);
        ARRAY_VIEW<int> av=segmented_curve->mesh.elements.Flattened();
        int next=0;
        for(int i=0;i<av.m;i++)
            if(map(av(i))<0){
                map(av(i))=next++;
                segmented_curve->particles.Add_Element();
                segmented_curve->particles.X.Last()=spline.particles.X(av(i));}
        av=map.Subset(av);}
    segmented_curve->Update_Number_Nodes();
    return segmented_curve;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> BEZIER_SPLINE<TV,d>* BEZIER_SPLINE<TV,d>::
Create()
{
    return new BEZIER_SPLINE;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> BEZIER_SPLINE<TV,d>* BEZIER_SPLINE<TV,d>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    return new BEZIER_SPLINE(particles);
}
namespace{
const int comb[9][9]=
{
    {1},
    {1,1},
    {1,2,1},
    {1,3,3,1},
    {1,4,6,4,1},
    {1,5,10,10,5,1},
    {1,6,15,20,15,6,1},
    {1,7,21,35,35,21,7,1},
    {1,8,28,56,70,56,28,8,1}
};
}
//#####################################################################
// Function Evaluate
//#####################################################################
template<class TV,int d> TV BEZIER_SPLINE<TV,d>::
Evaluate(int id,T t) const
{
    VECTOR<TV,d+1> X(particles.X.Subset(control_points(id)));
    T pow_t[d+1]={1};
    T pow_u[d+1]={1};
    T u=1-t;
    for(int i=0;i<d;i++){
        pow_t[i+1]=pow_t[i]*t;
        pow_u[i+1]=pow_u[i]*u;}
    TV z;
    for(int i=0;i<d+1;i++)
        z+=comb[d][i]*pow_t[i]*pow_u[d-i]*X(i);
    return z;
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV, int d> std::string BEZIER_SPLINE<TV,d>::
Name() const
{
    return Static_Name();
}
//#####################################################################
// Function Static_Name
//#####################################################################
template<class TV, int d> std::string BEZIER_SPLINE<TV,d>::
Static_Name()
{
    return LOG::sprintf("BEZIER_SPLINE<VECTOR<T,%d> ,%d>",TV::m,d);
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV, int d> std::string BEZIER_SPLINE<TV,d>::
Extension() const
{
    return Static_Extension();
}
//#####################################################################
// Function Static_Extension
//#####################################################################
template<class TV, int d> std::string BEZIER_SPLINE<TV,d>::
Static_Extension()
{
    return "";
}
namespace PhysBAM{
template class BEZIER_SPLINE<VECTOR<float,2>,3>;
template class BEZIER_SPLINE<VECTOR<float,3>,3>;
template class BEZIER_SPLINE<VECTOR<double,2>,3>;
template class BEZIER_SPLINE<VECTOR<double,3>,3>;
template TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<float,2>,1>::OBJECT* Create_Segmented_Curve<VECTOR<float,2>,3>(BEZIER_SPLINE<VECTOR<float,2>,3> const&,bool);
template TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<double,2>,1>::OBJECT* Create_Segmented_Curve<VECTOR<double,2>,3>(BEZIER_SPLINE<VECTOR<double,2>,3> const&,bool);
template void Smooth_Fit<VECTOR<float,2> >(BEZIER_SPLINE<VECTOR<float,2>,3>&,ARRAY_VIEW<VECTOR<float,2>,int>);
template void Smooth_Fit<VECTOR<double,2> >(BEZIER_SPLINE<VECTOR<double,2>,3>&,ARRAY_VIEW<VECTOR<double,2>,int>);
}
