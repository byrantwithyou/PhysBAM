//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class B_SPLINE_PATCH
//##################################################################### 
#include <Tools/Arrays/IDENTITY_ARRAY.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Matrices/BANDED_MATRIX.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE_PATCH.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> B_SPLINE_PATCH<TV,d>::
B_SPLINE_PATCH()
    :need_destroy_particles(true),particles(*new GEOMETRY_PARTICLES<TV>)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> B_SPLINE_PATCH<TV,d>::
B_SPLINE_PATCH(GEOMETRY_PARTICLES<TV>& particles)
    :need_destroy_particles(false),particles(particles)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> B_SPLINE_PATCH<TV,d>::
~B_SPLINE_PATCH()
{
    if(need_destroy_particles) delete &particles;
}
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV,int d> B_SPLINE_PATCH<TV,d>* B_SPLINE_PATCH<TV,d>::
Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const
{
    B_SPLINE_PATCH* bs=new B_SPLINE_PATCH(new_particles);
    int offset=new_particles.Size();
    new_particles.Append(particles);
    if(particle_indices) particle_indices->Append_Elements(IDENTITY_ARRAY<>(particles.Size())+offset);
    bs->control_points=control_points+offset;
    bs->knots_t=knots_t;
    bs->knots_s=knots_s;
    return bs;
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV,int d> void B_SPLINE_PATCH<TV,d>::
Read(TYPED_ISTREAM& input)
{
    int size;
    Read_Binary(input,control_points,knots_t,knots_s,size);
    particles.Clean_Memory(); // strip everything away except for position
    particles.Resize(size);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,particles.X.Get_Array_Pointer(),size);
    else Read_Binary_Array<float>(input.stream,particles.X.Get_Array_Pointer(),size);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV,int d> void B_SPLINE_PATCH<TV,d>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,control_points,knots_t,knots_s,particles.X);
}
//#####################################################################
// Function Smooth_Fit
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit(B_SPLINE_PATCH<TV,3>& bs,const ARRAY<ARRAY<TV> >& X)
{
    typedef typename TV::SCALAR T;
    int M=X.m; int N=X(0).m;

    bs.control_points.Resize(M+2);
    for(int i=0;i<M+2;i++) bs.control_points(i)=IDENTITY_ARRAY<>(N+2)+i*(N+2);

    bs.knots_s.Append(0); bs.knots_s.Append(0);
    for(int i=0;i<M;i++) bs.knots_s.Append((T)i/(M-1));
    bs.knots_s.Append(1); bs.knots_s.Append(1);
    bs.knots_t.Append(0); bs.knots_t.Append(0);
    for(int i=0;i<N;i++) bs.knots_t.Append((T)i/(N-1));
    bs.knots_t.Append(1); bs.knots_t.Append(1);

    bs.particles.Add_Elements((M+2)*(N+2));

    ARRAY<ARRAY<TV> > rhs1(M);

    for(int i=0;i<M;i++) rhs1(i).Resize(N);

    BANDED_MATRIX<T,3> mat(M);
    mat.diagonal_column=1;
    for(int i=0;i<M;i++)
        for(int j=0;j<N;j++)
            rhs1(i)(j)=X(i)(j)*12;

    mat.A(0)=VECTOR<T,3>(0,18,-6);
    if(M<=3){
        PHYSBAM_ASSERT(M>=2);
        if(M==3) mat.A(1)=VECTOR<T,3>(3,6,3);}
    else{
        mat.A(1)=VECTOR<T,3>(3,7,2);
        for(int i=2;i<M-2;i++) mat.A(i)=VECTOR<T,3>(2,8,2);
        mat.A(M-2)=VECTOR<T,3>(2,7,3);}
    mat.A(M-1)=VECTOR<T,3>(-6,18,0);

    mat.QR_Solve(rhs1); // Vertical solve to get intermediate points

    ARRAY<ARRAY<TV> > rhs2(N);
    for(int i=0;i<N;i++) rhs2(i).Resize(M+2);

    for(int j=0;j<N;j++){
        rhs2(j)(0)=X(0)(j)*12;
        for(int i=1;i<=M;i++) rhs2(j)(i)=rhs1(i-1)(j)*12;
        rhs2(j)(M+1)=X(M-1)(j)*12;}

    mat.Resize(N);
    mat.diagonal_column=1;

    mat.A(0)=VECTOR<T,3>(0,18,-6);
    if(N<=3){
        PHYSBAM_ASSERT(N>=2);
        if(N==3) mat.A(1)=VECTOR<T,3>(3,6,3);}
    else{
        mat.A(1)=VECTOR<T,3>(3,7,2);
        for(int i=2;i<N-2;i++) mat.A(i)=VECTOR<T,3>(2,8,2);
        mat.A(N-2)=VECTOR<T,3>(2,7,3);}
    mat.A(N-1)=VECTOR<T,3>(-6,18,0);

    for(int i=0;i<M+2;i++){
        bs.particles.X(bs.control_points(i)(0))=rhs2(0)(i)/12;
        bs.particles.X(bs.control_points(i)(N+1))=rhs2(N-1)(i)/12;}
 
    mat.QR_Solve(rhs2); // Horizontal solve to get final control points

    for(int i=0;i<M+2;i++)
        for(int j=1;j<=N;j++)
            bs.particles.X(bs.control_points(i)(j))=rhs2(j-1)(i);

// knots_t: 0, 0, 0, 1/(N-1), 2/(N-1), ... , (N-2)/(N-1), 1, 1, 1
// knots_s: 0, 0, 0, 1/(M-1), 2/(M-1), ... , (M-2)/(M-1), 1, 1, 1
//
// control_points:
// 0, 1, 2, ..., N+1
// N+2, ......., 2N+3
// .              .
// .              .
// .              .
// ............, (M+2)(N+2)-1
}
//#####################################################################
// Function Smooth_Fit_Loop
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit_Loop(B_SPLINE_PATCH<TV,3>& bs,const ARRAY<ARRAY<TV> >& X)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Create_Segmented_Curve
//#####################################################################
template<class TV,int d> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT* PhysBAM::
Create_Segmented_Curve(const B_SPLINE_PATCH<TV,d>& spline,bool same_particles)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> B_SPLINE_PATCH<TV,d>* B_SPLINE_PATCH<TV,d>::
Create()
{
    return new B_SPLINE_PATCH;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> B_SPLINE_PATCH<TV,d>* B_SPLINE_PATCH<TV,d>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    return new B_SPLINE_PATCH(particles);
}
//#####################################################################
// Function Evaluate
//#####################################################################
template<class TV,int d> TV B_SPLINE_PATCH<TV,d>::
Evaluate(T s, T t) const
{
    t=clamp(t,knots_t(d-1),knots_t(knots_t.m-d));
    int id_t=std::upper_bound(knots_t.begin(),knots_t.end()-d,t)-knots_t.begin();
    s=clamp(s,knots_s(d-1),knots_s(knots_s.m-d));
    int id_s=std::upper_bound(knots_s.begin(),knots_s.end()-d,s)-knots_s.begin();

    // First solve along t direction...
    TV x[d+1][d+1];
    TV y[d+1][d+1];
    for(int j=0;j<=d;j++){
        for(int i=0;i<=d;i++)
            x[0][i]=particles.X(control_points(j-d+id_s)(i-d+id_t));
        for(int k=0;k<d;k++)
            for(int i=k+1;i<=d;i++){
                T u0=knots_t(i-d+id_t-1),u1=knots_t(i-k+id_t-1),a=u1!=u0?(t-u0)/(u1-u0):0;
                x[k+1][i]=(1-a)*x[k][i-1]+a*x[k][i];}
        y[0][j]=x[d][d];}

    // ...then use those outputs as starting values and solve in s direction.
    for(int k=0;k<d;k++)
        for(int j=k+1;j<=d;j++){
            T u0=knots_s(j-d+id_s-1),u1=knots_s(j-k+id_s-1),a=u1!=u0?(s-u0)/(u1-u0):0;
            y[k+1][j]=(1-a)*y[k][j-1]+a*y[k][j];}

    return y[d][d];
}
//#####################################################################
// Function Fill_Bezier
//#####################################################################
template<class TV> void PhysBAM::
Fill_Bezier(BEZIER_SPLINE_PATCH<TV,3>& bez,const B_SPLINE_PATCH<TV,3>& bs)
{
    typedef typename TV::SCALAR T;

    int segs_t=bs.knots_t.m-5;
    int segs_s=bs.knots_s.m-5;




    bez.particles.Resize((3*segs_s+1)*(3*segs_t+1));
    if(bez.control_points.m!=segs_s*segs_t){
    // Need to reconstruct control_points.
        bez.control_points.Resize(segs_s*segs_t);
        for(int id=0;id<segs_s*segs_t;id++)
            for(int i=0;i<4;i++)
                for(int j=0;j<4;j++)
                    bez.control_points(id)(4*i+j)=(3*(id/segs_t)+i)*(3*segs_t+1)+3*(id%segs_t)+j;}

    ARRAY<ARRAY<TV> > intermediate(3*(bs.knots_s.m-5)+1);
    for(int i=0;i<3*(bs.knots_s.m-5)+1;i++)
        intermediate(i).Resize(bs.knots_t.m-2);

    // Fill intermediate using vertical solves.
    for(int i=0;i<bs.knots_s.m-5;i++){
        T k23=bs.knots_s(i+2)-bs.knots_s(i+1);
        T k14=bs.knots_s(i+3)-bs.knots_s(i);
        T k24=bs.knots_s(i+3)-bs.knots_s(i+1);
        T k34=bs.knots_s(i+3)-bs.knots_s(i+2);
        T k15=bs.knots_s(i+4)-bs.knots_s(i);
        T k25=bs.knots_s(i+4)-bs.knots_s(i+1);
        T k35=bs.knots_s(i+4)-bs.knots_s(i+2);
        T k45=bs.knots_s(i+4)-bs.knots_s(i+3);
        
        T a11=k34*k34/(k14*k24);
        T a13=k23*k23/(k25*k24);
        T a12=1-a11-a13;
        T a22=k35/k25;
        T a23=k23/k25;
        T a32=k45/k25;
        T a33=k24/k25;
        
        for(int j=0;j<bs.knots_t.m-2;j++){
            intermediate(3*i)(j)=a11*bs.particles.X(bs.control_points(i)(j))
                +a12*bs.particles.X(bs.control_points(i+1)(j))+a13*bs.particles.X(bs.control_points(i+2)(j));
            
            intermediate(3*i+1)(j)=a22*bs.particles.X(bs.control_points(i+1)(j))
                +a23*bs.particles.X(bs.control_points(i+2)(j));
            
            intermediate(3*i+2)(j)=a32*bs.particles.X(bs.control_points(i+1)(j))
                +a33*bs.particles.X(bs.control_points(i+2)(j));}
    }
    int m=bs.knots_s.m-4;
    T k25=bs.knots_s(m+2)-bs.knots_s(m-1);
    T k35=bs.knots_s(m+2)-bs.knots_s(m);
    T k45=bs.knots_s(m+2)-bs.knots_s(m+1);
    T k34=bs.knots_s(m+1)-bs.knots_s(m);
    T k36=bs.knots_s(m+3)-bs.knots_s(m);
    T a42=k45*k45/(k35*k25);
    T a44=k34*k34/(k36*k35);
    T a43=1-a42-a44;
    for(int j=0;j<bs.knots_t.m-2;j++)
        intermediate(3*(bs.knots_s.m-5))(j)=a42*bs.particles.X(bs.control_points(m-1)(j))
            +a43*bs.particles.X(bs.control_points(m)(j))+a44*bs.particles.X(bs.control_points(m+1)(j));
    
    // Now cash in intermediate for final Bezier points by working horizontally.
    for(int j=0;j<bs.knots_t.m-5;j++)
    {
        T k23=bs.knots_t(j+2)-bs.knots_t(j+1);
        T k14=bs.knots_t(j+3)-bs.knots_t(j);
        T k24=bs.knots_t(j+3)-bs.knots_t(j+1);
        T k34=bs.knots_t(j+3)-bs.knots_t(j+2);
        T k15=bs.knots_t(j+4)-bs.knots_t(j);
        T k25=bs.knots_t(j+4)-bs.knots_t(j+1);
        T k35=bs.knots_t(j+4)-bs.knots_t(j+2);
        T k45=bs.knots_t(j+4)-bs.knots_t(j+3);
        
        T a11=k34*k34/(k14*k24);
        T a13=k23*k23/(k25*k24);
        T a12=1-a11-a13;
        T a22=k35/k25;
        T a23=k23/k25;
        T a32=k45/k25;
        T a33=k24/k25;
        for(int i=0;i<3*(bs.knots_s.m-5)+1;i++){
            int id=i/3 * (bs.knots_t.m-5)+j;
            if(id>=bez.control_points.m){
                id-=bs.knots_t.m-5;
                bez.particles.X(bez.control_points(id)(4*(3)))=a11*intermediate(i)(j)
                    +a12*intermediate(i)(j+1)+a13*intermediate(i)(j+2);
                bez.particles.X(bez.control_points(id)(4*(3)+1))=a22*intermediate(i)(j+1)
                    +a23*intermediate(i)(j+2);
                bez.particles.X(bez.control_points(id)(4*(3)+2))=a32*intermediate(i)(j+1)
                    +a33*intermediate(i)(j+2);}
            else{
                bez.particles.X(bez.control_points(id)(4*(i%3)))=a11*intermediate(i)(j)
                    +a12*intermediate(i)(j+1)+a13*intermediate(i)(j+2);
                bez.particles.X(bez.control_points(id)(4*(i%3)+1))=a22*intermediate(i)(j+1)
                    +a23*intermediate(i)(j+2);
                bez.particles.X(bez.control_points(id)(4*(i%3)+2))=a32*intermediate(i)(j+1)
                    +a33*intermediate(i)(j+2);}}
    }
    m=bs.knots_t.m-4;
    k25=bs.knots_t(m+2)-bs.knots_t(m-1);
    k35=bs.knots_t(m+2)-bs.knots_t(m);
    k45=bs.knots_t(m+2)-bs.knots_t(m+1);
    k34=bs.knots_t(m+1)-bs.knots_t(m);
    k36=bs.knots_t(m+3)-bs.knots_t(m);
    a42=k45*k45/(k35*k25);
    a44=k34*k34/(k36*k35);
    a43=1-a42-a44;
    for(int i=0;i<3*(bs.knots_s.m-5)+1;i++){
        int id=i/3 * (bs.knots_t.m-5)+bs.knots_t.m-6;
        if(id>=bez.control_points.m)
            bez.particles.X(bez.control_points(id-(bs.knots_t.m-5))(15))=a42*intermediate(i)(m-1);
        else
            bez.particles.X(bez.control_points(id)(4*(i%3)+3))=a42*intermediate(i)(m-1)
                +a43*intermediate(i)(m)+a44*intermediate(i)(m+1);}
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV, int d> std::string B_SPLINE_PATCH<TV,d>::
Name() const
{
    return Static_Name();
}
//#####################################################################
// Function Static_Name
//#####################################################################
template<class TV, int d> std::string B_SPLINE_PATCH<TV,d>::
Static_Name()
{
    return STRING_UTILITIES::string_sprintf("B_SPLINE_PATCH<VECTOR<T,%d> ,%d>",TV::dimension,d);
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV, int d> std::string B_SPLINE_PATCH<TV,d>::
Extension() const
{
    return Static_Extension();
}
//#####################################################################
// Function Static_Extension
//#####################################################################
template<class TV, int d> std::string B_SPLINE_PATCH<TV,d>::
Static_Extension()
{
    return "";
}
namespace PhysBAM{
template class B_SPLINE_PATCH<VECTOR<float,2>,3>;
template class B_SPLINE_PATCH<VECTOR<float,3>,3>;
template class B_SPLINE_PATCH<VECTOR<double,2>,3>;
template class B_SPLINE_PATCH<VECTOR<double,3>,3>;

template void Smooth_Fit<VECTOR<float,2> >(B_SPLINE_PATCH<VECTOR<float,2>,3>&,ARRAY<ARRAY<VECTOR<float,2>,int>,int> const&);
template void Smooth_Fit<VECTOR<double,2> >(B_SPLINE_PATCH<VECTOR<double,2>,3>&,ARRAY<ARRAY<VECTOR<double,2>,int>,int> const&);
template void Fill_Bezier<VECTOR<float,2> >(BEZIER_SPLINE_PATCH<VECTOR<float,2>,3>&,B_SPLINE_PATCH<VECTOR<float,2>,3> const&);
template void Fill_Bezier<VECTOR<double,2> >(BEZIER_SPLINE_PATCH<VECTOR<double,2>,3>&,B_SPLINE_PATCH<VECTOR<double,2>,3> const&);
template void Smooth_Fit<VECTOR<float,3> >(B_SPLINE_PATCH<VECTOR<float,3>,3>&,ARRAY<ARRAY<VECTOR<float,3>,int>,int> const&);
template void Smooth_Fit<VECTOR<double,3> >(B_SPLINE_PATCH<VECTOR<double,3>,3>&,ARRAY<ARRAY<VECTOR<double,3>,int>,int> const&);
template void Fill_Bezier<VECTOR<float,3> >(BEZIER_SPLINE_PATCH<VECTOR<float,3>,3>&,B_SPLINE_PATCH<VECTOR<float,3>,3> const&);
template void Fill_Bezier<VECTOR<double,3> >(BEZIER_SPLINE_PATCH<VECTOR<double,3>,3>&,B_SPLINE_PATCH<VECTOR<double,3>,3> const&);
}
