//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class B_SPLINE
//##################################################################### 
#include <Tools/Arrays/IDENTITY_ARRAY.h>
#include <Tools/Matrices/BANDED_MATRIX.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE.h>
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
    // This always uses a uniform knot distribution in parameter space.
    // Maybe in future add a version that uses some other structure?
    typedef typename TV::SCALAR T;
    ARRAY<TV> rhs(X.m); // We have X.m distinct knot points, so X.m+2 control points need to be found. First and last will be automatic though.
    BANDED_MATRIX<T,3> matrix(X.m);
    matrix.diagonal_column=1;
    for(int i=0;i<X.m;i++) rhs(i)=X(i)*12;

    matrix.A(0)=VECTOR<T,3>(0,18,-6);
    if(X.m<=3){
        PHYSBAM_ASSERT(X.m>=2);
        if(X.m==3) matrix.A(1)=VECTOR<T,3>(3,6,3);}
    else{
        matrix.A(1)=VECTOR<T,3>(3,7,2);
        for(int i=2;i<X.m-2;i++) matrix.A(i)=VECTOR<T,3>(2,8,2);
        matrix.A(X.m-2)=VECTOR<T,3>(2,7,3);}
    matrix.A(X.m-1)=VECTOR<T,3>(-6,18,0);

    matrix.QR_Solve(rhs);
    bs.knots.Append(0);
    bs.knots.Append(0);
    bs.knots.Append(0);
    bs.particles.Add_Elements(X.m+2);
    bs.control_points=IDENTITY_ARRAY<>(X.m+2);
    bs.particles.X(0)=X(0);
    for(int i=1;i<=X.m;i++){
        bs.particles.X(i)=rhs(i-1);
        bs.knots.Append((T)i/(X.m-1));}
    bs.particles.X.Last()=X.Last();
    bs.knots(X.m+2)=1;
    bs.knots.Append(1);
    
    // Final knots list: 0, 0, 0, 1/n, 2/n, ..., (n-2)/(n-1), 1, 1, 1.
    // Final control points list: 0,1,2,...,n+1
}
//#####################################################################
// Function Smooth_Fit_Loop
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit_Loop(B_SPLINE<TV,3>& bs,ARRAY_VIEW<TV> X)
{
    // X.m knots, equally spaced
    // matrix with 2/3 on diagonal, 1/6 above&below and in corners
    // How do I solve a matrix like that in PhysBAM?
    // Also: the current Evaluate function will cry if we give it that knot list.
    typedef typename TV::SCALAR T;
    ARRAY<TV> rhs(X.m);
    BANDED_MATRIX<T,3> matrix(X.m);
    matrix.diagonal_column=1;
    for(int i=0;i<X.m;i++){matrix.A(i)=VECTOR<T,3>(1,4,1);rhs(i)=X(i)*6;}
    matrix.QR_Solve(rhs);

    for(int k=-2;k<X.m+3;k++) bs.knots.Append((T)k/(X.m));
    bs.particles.Add_Elements(rhs.m);
    bs.particles.X=rhs;
    bs.control_points=IDENTITY_ARRAY<>(X.m);
    for(int j=0;j<3;j++) bs.control_points.Append(j);
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
    t=clamp(t,knots(d-1),knots(knots.m-d));
    
    int id=std::upper_bound(knots.begin(),knots.end()-d,t)-knots.begin();
//    PHYSBAM_ASSERT(id>=d);
    TV x[d+1][d+1];
    for(int i=0;i<=d;i++){
        x[0][i]=particles.X(control_points(i-d+id));}
    id--;
    for(int k=0;k<d;k++)
        for(int i=k+1;i<=d;i++){
            T u0=knots(i-d+id),u1=knots(i-k+id),a=u1!=u0?(t-u0)/(u1-u0):0;
            x[k+1][i]=(1-a)*x[k][i-1]+a*x[k][i];}

    return x[d][d];
}
//#####################################################################
// Function Fill_Bezier
//#####################################################################
template<class TV> void PhysBAM::
Fill_Bezier(BEZIER_SPLINE<TV,3>& bez,const B_SPLINE<TV,3>& bs)
{
    typedef typename TV::SCALAR T;
    int m=bs.knots.m-4; // m = number of Bezier segments + 1. (i.e. number of `physical' Bezier points)
    bez.particles.Resize(m*3-2);
    if(bez.control_points.m<m-1){
        for(int i=bez.control_points.m;i<m-1;i++)
            bez.control_points.Append(VECTOR<int,4>(3*i,3*i+1,3*i+2,3*i+3));}
    else bez.control_points.Resize(m-1);
    for(int i=0;i<m-1;i++)
    {
        T k23=bs.knots(i+2)-bs.knots(i+1);
        T k14=bs.knots(i+3)-bs.knots(i);
        T k24=bs.knots(i+3)-bs.knots(i+1);
        T k34=bs.knots(i+3)-bs.knots(i+2);
        T k15=bs.knots(i+4)-bs.knots(i);
        T k25=bs.knots(i+4)-bs.knots(i+1);
        T k35=bs.knots(i+4)-bs.knots(i+2);
        T k45=bs.knots(i+4)-bs.knots(i+3);

        T a11=k34*k34/(k14*k24);
        T a13=k23*k23/(k25*k24);
        T a12=1-a11-a13;
        T a22=k35/k25;
        T a23=k23/k25;
        T a32=k45/k25;
        T a33=k24/k25;
        bez.particles.X(3*i)=a11*bs.particles.X(bs.control_points(i))
            +a12*bs.particles.X(bs.control_points(i+1))+a13*bs.particles.X(bs.control_points(i+2));
        bez.particles.X(3*i+1)=a22*bs.particles.X(bs.control_points(i+1))+a23*bs.particles.X(bs.control_points(i+2));
        bez.particles.X(3*i+2)=a32*bs.particles.X(bs.control_points(i+1))+a33*bs.particles.X(bs.control_points(i+2));
    }
    T k25=bs.knots(m+2)-bs.knots(m-1);
    T k35=bs.knots(m+2)-bs.knots(m);
    T k45=bs.knots(m+2)-bs.knots(m+1);
    T k34=bs.knots(m+1)-bs.knots(m);
    T k36=bs.knots(m+3)-bs.knots(m);
    T a42=k45*k45/(k35*k25);
    T a44=k34*k34/(k36*k35);
    T a43=1-a42-a44;
    bez.particles.X.Last()=a42*bs.particles.X(bs.control_points(m-1))
        +a43*bs.particles.X(bs.control_points(m))+a44*bs.particles.X(bs.control_points(m+1));
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV, int d> std::string B_SPLINE<TV,d>::
Name() const
{
    return Static_Name();
}
//#####################################################################
// Function Static_Name
//#####################################################################
template<class TV, int d> std::string B_SPLINE<TV,d>::
Static_Name()
{
    return STRING_UTILITIES::string_sprintf("B_SPLINE<VECTOR<T,%d> ,%d>",TV::dimension,d);
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV, int d> std::string B_SPLINE<TV,d>::
Extension() const
{
    return Static_Extension();
}
//#####################################################################
// Function Static_Extension
//#####################################################################
template<class TV, int d> std::string B_SPLINE<TV,d>::
Static_Extension()
{
    return "";
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
template void Smooth_Fit_Loop<VECTOR<float,2> >(B_SPLINE<VECTOR<float,2>,3>&,ARRAY_VIEW<VECTOR<float,2>,int>);
template void Smooth_Fit_Loop<VECTOR<double,2> >(B_SPLINE<VECTOR<double,2>,3>&,ARRAY_VIEW<VECTOR<double,2>,int>);
template void Fill_Bezier<VECTOR<float,2> >(BEZIER_SPLINE<VECTOR<float,2>,3>&,B_SPLINE<VECTOR<float,2>,3> const&);
template void Fill_Bezier<VECTOR<double,2> >(BEZIER_SPLINE<VECTOR<double,2>,3>&,B_SPLINE<VECTOR<double,2>,3> const&);
}
