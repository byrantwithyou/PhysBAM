//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class B_SPLINE
//##################################################################### 
#include <Tools/Arrays/IDENTITY_ARRAY.h>
#include <Tools/Matrices/BANDED_MATRIX.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> B_SPLINE<TV,d>::
B_SPLINE()
    :need_destroy_particles(true),particles(*new GEOMETRY_PARTICLES<TV>)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> B_SPLINE<TV,d>::
B_SPLINE(GEOMETRY_PARTICLES<TV>& particles)
    :need_destroy_particles(false),particles(particles)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> B_SPLINE<TV,d>::
~B_SPLINE()
{
    if(need_destroy_particles) delete &particles;
}
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV,int d> B_SPLINE<TV,d>* B_SPLINE<TV,d>::
Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const
{
    B_SPLINE* bs=new B_SPLINE(new_particles);
    int offset=new_particles.Size();
    new_particles.Append(particles);
    if(particle_indices) particle_indices->Append_Elements(IDENTITY_ARRAY<>(particles.Size())+offset);
    bs->control_points=control_points+offset;
    bs->knots=knots;
    return bs;
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV,int d> void B_SPLINE<TV,d>::
Read(TYPED_ISTREAM& input)
{
    int size;
    Read_Binary(input,control_points,knots,size);
    particles.Clean_Memory(); // strip everything away except for position
    particles.Resize(size);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,particles.X.Get_Array_Pointer(),size);
    else Read_Binary_Array<float>(input.stream,particles.X.Get_Array_Pointer(),size);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV,int d> void B_SPLINE<TV,d>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,control_points,knots,particles.X);
}
//#####################################################################
// Function Smooth_Fit
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit(B_SPLINE<TV,3>& bs,ARRAY_VIEW<TV> X)
{
    // typedef typename TV::SCALAR T;
    // ARRAY<TV> rhs(X.m*2-2);
    // BANDED_MATRIX<T,4> matrix(X.m*2-2);
    // matrix.diagonal_column=1;
    // matrix.A(0)=VECTOR<T,4>(0,2,-1,0);
    // rhs(0)=X(0);
    // for(int i=1;i<X.m-1;i++){
    //     matrix.A(2*i-1)=VECTOR<T,4>(-1,2,-2,1);
    //     matrix.A(2*i)=VECTOR<T,4>(1,1,0,0);
    //     rhs(2*i)=X(i)*2;}
    // matrix.A(X.m*2-3)=VECTOR<T,4>(-1,2,0,0);
    // rhs(X.m*2-3)=X(X.m-1);
    // matrix.QR_Solve(rhs);

    // bs.particles.Add_Elements(X.m*3-2);
    // VECTOR<int,4> base(0,1,2,3);
    // for(int i=0;i<X.m-1;i++){
    //     bs.particles.X(3*i)=X(i);
    //     bs.particles.X(3*i+1)=rhs(2*i);
    //     bs.particles.X(3*i+2)=rhs(2*i+1);
    //     bs.control_points.Append(base+3*i);}
    // bs.particles.X.Last()=X.Last();
}
//#####################################################################
// Function Smooth_Fit
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit_Loop(B_SPLINE<TV,3>& bs,ARRAY_VIEW<TV> X)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Create_Segmented_Curve
//#####################################################################
template<class TV,int d> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT* PhysBAM::
Create_Segmented_Curve(const B_SPLINE<TV,d>& spline,bool same_particles)
{
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT T_SEGMENTED_CURVE;
    T_SEGMENTED_CURVE* segmented_curve=0;
    if(same_particles) segmented_curve=T_SEGMENTED_CURVE::Create(spline.particles);
    else segmented_curve=T_SEGMENTED_CURVE::Create();
    for(int i=0;i<spline.control_points.m-1;i++){
        int a=spline.control_points(i),b=spline.control_points(i+1);
        if(a!=b)
            segmented_curve->mesh.elements.Append(VECTOR<int,2>(a,b));}

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
template<class TV,int d> B_SPLINE<TV,d>* B_SPLINE<TV,d>::
Create()
{
    return new B_SPLINE;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> B_SPLINE<TV,d>* B_SPLINE<TV,d>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    return new B_SPLINE(particles);
}
//#####################################################################
// Function Evaluate
//#####################################################################
template<class TV,int d> TV B_SPLINE<TV,d>::
Evaluate(T t) const
{
    t=clamp(t,knots(d),knots.Last());
    int id=std::upper_bound(knots.begin(),knots.end(),t)-knots.begin()-1;
    assert(id>=d);
    TV x[d+1][d+1];
    for(int i=0;i<=d;i++) x[0][i]=particles.X(control_points(i-d+id));

    for(int k=0;k<d;k++)
        for(int i=k+1;i<=d;i++){
            T u0=knots(i-d+id),u1=knots(i-k+id),a=(t-u0)/(u1-u0);
            x[k+1][i]=(1-a)*x[k][i-1]+a*x[k][i];}

    return x[d][d];
}
namespace PhysBAM{
template class B_SPLINE<VECTOR<float,2>,3>;
template class B_SPLINE<VECTOR<float,3>,3>;
template class B_SPLINE<VECTOR<double,2>,3>;
template class B_SPLINE<VECTOR<double,3>,3>;
template TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<float,2>,1>::OBJECT* Create_Segmented_Curve<VECTOR<float,2>,3>(B_SPLINE<VECTOR<float,2>,3> const&,bool);
template TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<double,2>,1>::OBJECT* Create_Segmented_Curve<VECTOR<double,2>,3>(B_SPLINE<VECTOR<double,2>,3> const&,bool);
template void Smooth_Fit<VECTOR<float,2> >(B_SPLINE<VECTOR<float,2>,3>&,ARRAY_VIEW<VECTOR<float,2>,int>);
template void Smooth_Fit<VECTOR<double,2> >(B_SPLINE<VECTOR<double,2>,3>&,ARRAY_VIEW<VECTOR<double,2>,int>);
}
