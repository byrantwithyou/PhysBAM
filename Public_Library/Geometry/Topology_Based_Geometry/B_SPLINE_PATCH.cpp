//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class B_SPLINE_PATCH
//##################################################################### 
#include <Tools/Arrays/IDENTITY_ARRAY.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/BANDED_MATRIX.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE_PATCH.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
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
    for(RANGE_ITERATOR<IV::m> it(control_points.domain);it.Valid();it.Next())
        bs->control_points(it.index)+=offset;
    bs->knots_s=knots_s;
    bs->knots_t=knots_t;
    bs->loop_s=loop_s;
    bs->loop_t=loop_t;
    bs->m=m;
    return bs;
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV,int d> void B_SPLINE_PATCH<TV,d>::
Read(TYPED_ISTREAM& input)
{
    int size;
    Read_Binary(input,m,loop_s,loop_t,control_points,knots_s,knots_t,size);
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
    Write_Binary(output,m,loop_s,loop_t,control_points,knots_s,knots_t,particles.X);
}
//#####################################################################
// Function Smooth_Fit
//#####################################################################
template<class TV> void PhysBAM::
Smooth_Fit(B_SPLINE_PATCH<TV,3>& bs,const ARRAY<TV,VECTOR<int,2>>& X,bool loop_s,bool loop_t)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,2> IV;
    int M=X.domain.max_corner(0); int N=X.domain.max_corner(1);

    bs.control_points.Resize(RANGE<IV>(IV(),IV(M+(loop_s?3:2),N+(loop_t?3:2))));
    const IV max=bs.control_points.domain.max_corner;
    for(int i=0;i<max(0);i++)
        for(int j=0;j<max(1);j++){
            int ii=loop_s?i%M:i;
            int jj=loop_t?j%N:j;
            bs.control_points(i,j)=jj+ii*(N+(loop_t?0:2));}

    for(int i=-2;i<M+2+loop_s;i++){
        T s=(T)i/(M-1+loop_s);
        if(!loop_s) s=clamp(s,(T)0,(T)1);
        bs.knots_s.Append(s);}
    for(int i=-2;i<N+2+loop_t;i++){
        T t=(T)i/(N-1+loop_t);
        if(!loop_t) t=clamp(t,(T)0,(T)1);
        bs.knots_t.Append(t);}

    bs.particles.Add_Elements((M+(loop_s?0:2))*(N+(loop_t?0:2)));

    ARRAY<ARRAY<TV> > rhs1(M);
    for(int i=0;i<M;i++) rhs1(i).Resize(N);

    BANDED_MATRIX<T,3> mat(M);
    mat.diagonal_column=1;
    mat.wrap=loop_s;
    for(int i=0;i<M;i++)
        for(int j=0;j<N;j++)
            rhs1(i)(j)=X(i,j)*12;

    if(!loop_s){
        mat.A(0)=VECTOR<T,3>(0,18,-6);
        if(M<=3){
            PHYSBAM_ASSERT(M>=2);
            if(M==3) mat.A(1)=VECTOR<T,3>(3,6,3);}
        else{
            mat.A(1)=VECTOR<T,3>(3,7,2);
            for(int i=2;i<M-2;i++) mat.A(i)=VECTOR<T,3>(2,8,2);
            mat.A(M-2)=VECTOR<T,3>(2,7,3);}
        mat.A(M-1)=VECTOR<T,3>(-6,18,0);}
    else {
        for(int i=0;i<M;i++) mat.A(i)=VECTOR<T,3>(2,8,2);}

    mat.QR_Solve(rhs1); // Vertical solve to get intermediate points

    ARRAY<ARRAY<TV> > rhs2(N);
    for(int i=0;i<N;i++) rhs2(i).Resize(M+(loop_s?0:2));

    if(!loop_s){
        for(int j=0;j<N;j++){
            rhs2(j)(0)=X(0,j)*12;
            for(int i=1;i<=M;i++) rhs2(j)(i)=rhs1(i-1)(j)*12;
            rhs2(j)(M+1)=X(M-1,j)*12;}}
    else {
        for(int j=0;j<N;j++)
            for(int i=0;i<M;i++)
                rhs2(j)(i)=rhs1(i)(j)*12;}
 

    mat.Resize(N);
    mat.diagonal_column=1;
    mat.wrap=loop_t;

    if(!loop_t){
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
            bs.particles.X(bs.control_points(i,0))=rhs2(0)(i)/12;
            bs.particles.X(bs.control_points(i,N+1))=rhs2(N-1)(i)/12;}}
    else {
        for(int i=0;i<N;i++) mat.A(i)=VECTOR<T,3>(2,8,2);}

    mat.QR_Solve(rhs2); // Horizontal solve to get final control points

    for(int i=0;i<M+(loop_s?0:2);i++)
        for(int j=0;j<N;j++)
            bs.particles.X(bs.control_points(i,j+1-loop_t))=rhs2(j)(i);
    bs.m=(bs.knots_s.m-5)*(bs.knots_t.m-5);
    bs.loop_s=loop_s;
    bs.loop_t=loop_t;
// knots_s: 0, 0, 0, 1/(M-1), 2/(M-1), ... , (M-2)/(M-1), 1, 1, 1
// knots_t: 0, 0, 0, 1/(N-1), 2/(N-1), ... , (N-2)/(N-1), 1, 1, 1
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
// Function Create_Triangulated_Surface
//#####################################################################
template<class TV,int d> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,2>::OBJECT* PhysBAM::
Create_Triangulated_Surface(const B_SPLINE_PATCH<TV,d>& spline,bool same_particles)
{
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,2>::OBJECT T_TRIANGULATED_SURFACE;
    T_TRIANGULATED_SURFACE* tri_patch=0;
    if(same_particles) tri_patch=T_TRIANGULATED_SURFACE::Create(spline.particles);
    else tri_patch=T_TRIANGULATED_SURFACE::Create();
    for(int i=0;i<spline.control_points.domain.max_corner(0)-1;i++){
        for(int j=0;j<spline.control_points.domain.max_corner(1)-1;j++){
            int a00=spline.control_points(i,j),a01=spline.control_points(i,j+1),
                a10=spline.control_points(i+1,j),a11=spline.control_points(i+1,j+1);
            tri_patch->mesh.elements.Append(VECTOR<int,3>(a00,a01,a11));
            tri_patch->mesh.elements.Append(VECTOR<int,3>(a00,a11,a10));}}

    if(!same_particles){
        ARRAY<int> map(spline.particles.X.m,true,-1);
        ARRAY_VIEW<int> av=tri_patch->mesh.elements.Flattened();
        int next=0;
        LOG::printf("av.m = %P\n",av.m);
        for(int i=0;i<av.m;i++)
            if(map(av(i))<0){
                map(av(i))=next++;
                tri_patch->particles.Add_Element();
                tri_patch->particles.X.Last()=spline.particles.X(av(i));}
        av=map.Subset(av);}
    tri_patch->Update_Number_Nodes();
    tri_patch->Print_Statistics(LOG::cout,.1);
    return tri_patch;
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
Evaluate(T s,T t,VECTOR<TV,2>* tangents) const
{
    t=clamp(t,knots_t(d-1),knots_t(knots_t.m-d));
    int id_t=std::upper_bound(knots_t.begin(),knots_t.end()-d,t)-knots_t.begin();
    s=clamp(s,knots_s(d-1),knots_s(knots_s.m-d));
    int id_s=std::upper_bound(knots_s.begin(),knots_s.end()-d,s)-knots_s.begin();

    // First solve along t direction...
    TV x[d+1];
    TV y[d+1];
    TV x_t[d+1];
    TV y_s[d+1];
    TV y_t[d+1];

    for(int j=0;j<=d;j++){
        for(int i=0;i<=d;i++){
            x[i]=particles.X(control_points(IV(j-d+id_s,i-d+id_t)));
            x_t[i]=TV();}
        for(int k=0;k<d;k++)
            for(int i=d;i>=k+1;i--){
                T u0=knots_t(i-d+id_t-1),u1=knots_t(i-k+id_t-1),a=u1!=u0?(t-u0)/(u1-u0):0;
                T a_t=u1!=u0?1/(u1-u0):0;
                x_t[i]=(-a_t)*x[i-1]+a_t*x[i]+(1-a)*x_t[i-1]+a*x_t[i];
                x[i]=(1-a)*x[i-1]+a*x[i];}
        y_t[j]=x_t[d];
        y[j]=x[d];}

    // ...then use those outputs as starting values and solve in s direction.
    for(int k=0;k<d;k++)
        for(int j=d;j>=k+1;j--){
            T u0=knots_s(j-d+id_s-1),u1=knots_s(j-k+id_s-1),a=u1!=u0?(s-u0)/(u1-u0):0;
            T a_s=u1!=u0?1/(u1-u0):0;
            y_s[j]=(-a_s)*y[j-1]+a_s*y[j]+(1-a)*y_s[j-1]+a*y_s[j];
            y_t[j]=(1-a)*y_t[j-1]+a*y_t[j];
            y[j]=(1-a)*y[j-1]+a*y[j];}

    if(tangents) *tangents=VECTOR<TV,2>(y_s[d],y_t[d]);
    return y[d];
}
template<class T,int w1,int w2> void Vec_Outer(VECTOR<T,w1*w2>& u, const VECTOR<T,w1>& v1,const VECTOR<T,w2>& v2)
{
    for(int i=0;i<w1;i++)
        for(int j=0;j<w2;j++)
            u(w2*i+j)=v1(i)*v2(j);

}
//#####################################################################
// Function Calculate_Weights
//#####################################################################
template<class TV, int d> void B_SPLINE_PATCH<TV,d>::
Calculate_Weights(T s,T t,WV& w,WV& dw0,WV& dw1,WV& ddw0,WV& ddw1,WV& ddw2) const
{
    // ddw is d^2w/ds^2 d, d^2w/dsdt, d^2w/dt^2
    t=clamp(t,knots_t(d-1),knots_t(knots_t.m-d));
    int id_t=std::upper_bound(knots_t.begin(),knots_t.end()-d,t)-knots_t.begin();
    s=clamp(s,knots_s(d-1),knots_s(knots_s.m-d));
    int id_s=std::upper_bound(knots_s.begin(),knots_s.end()-d,s)-knots_s.begin();
    typedef VECTOR<T,d+1> WV1;

    WV1 x[d+1];
    WV1 x_s[d+1];
    WV1 x_ss[d+1];

    for(int i=0;i<=d;i++)
        x[i]=WV1::Axis_Vector(i);
    for(int k=0;k<d;k++)
        for(int i=d;i>=k+1;i--){
            T u0=knots_s(i-d+id_s-1),u1=knots_s(i-k+id_s-1),a=u1!=u0?(s-u0)/(u1-u0):0;
            T a_s=u1!=u0?1/(u1-u0):0;
            x_ss[i]=(-2*a_s)*x_s[i-1]+(2*a_s)*x_s[i]+(1-a)*x_ss[i-1]+a*x_ss[i];
            x_s[i]=(-a_s)*x[i-1]+a_s*x[i]+(1-a)*x_s[i-1]+a*x_s[i];
            x[i]=(1-a)*x[i-1]+a*x[i];}

    WV1 y[d+1];
    WV1 y_t[d+1];
    WV1 y_tt[d+1];

    for(int j=0;j<=d;j++)
        y[j]=WV1::Axis_Vector(j);
    for(int k=0;k<d;k++)
        for(int j=d;j>=k+1;j--){
            T u0=knots_t(j-d+id_t-1),u1=knots_t(j-k+id_t-1),a=u1!=u0?(t-u0)/(u1-u0):0;
            T a_t=u1!=u0?1/(u1-u0):0;
            y_tt[j]=(-2*a_t)*y_t[j-1]+(2*a_t)*y_t[j]+(1-a)*y_tt[j-1]+a*y_tt[j];
            y_t[j]=(-a_t)*y[j-1]+a_t*y[j]+(1-a)*y_t[j-1]+a*y_t[j];
            y[j]=(1-a)*y[j-1]+a*y[j];}
    Vec_Outer(w,x[d],y[d]);
    Vec_Outer(dw0,x_s[d],y[d]);
    Vec_Outer(dw1,x[d],y_t[d]);
    Vec_Outer(ddw0,x_ss[d],y[d]);
    Vec_Outer(ddw1,x_s[d],y_t[d]);
    Vec_Outer(ddw2,x[d],y_tt[d]);
}
//#####################################################################
// Function Control_Points_For_Element
//#####################################################################
template<class TV, int d> VECTOR<int,(d+1)*(d+1)> B_SPLINE_PATCH<TV,d>::
Control_Points_For_Element(int element) const
{
    VECTOR<int,(d+1)*(d+1)> a;
    int i=element/(knots_t.m-d-2);
    int j=element%(knots_t.m-d-2);
    const VECTOR<int,2> max=control_points.domain.max_corner;
    for(int ii=0;ii<=d;ii++)
        for(int jj=0;jj<=d;jj++)
            a((d+1)*ii+jj)=control_points((i+ii)%max(0),(j+jj)%max(1));
    return a;
}
//#####################################################################
// Function Range_For_Element
//#####################################################################
template<class TV, int d> RANGE<VECTOR<typename TV::SCALAR,2> > B_SPLINE_PATCH<TV,d>::
Range_For_Element(int element) const
{
    int i=element/(knots_t.m-d-2);
    int j=element%(knots_t.m-d-2);
    return RANGE<VECTOR<T,2>>(VECTOR<T,2>(knots_s(i+d-1), knots_t(j+d-1)), VECTOR<T,2>(knots_s(i+d), knots_t(j+d)));
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
            intermediate(3*i)(j)=a11*bs.particles.X(bs.control_points(i,j))
                +a12*bs.particles.X(bs.control_points(i+1)(j))+a13*bs.particles.X(bs.control_points(i+2)(j));
            
            intermediate(3*i+1)(j)=a22*bs.particles.X(bs.control_points(i+1,j))
                +a23*bs.particles.X(bs.control_points(i+2,j));
            
            intermediate(3*i+2)(j)=a32*bs.particles.X(bs.control_points(i+1,j))
                +a33*bs.particles.X(bs.control_points(i+2,j));}
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
        intermediate(3*(bs.knots_s.m-5))(j)=a42*bs.particles.X(bs.control_points(m-1,j))
            +a43*bs.particles.X(bs.control_points(m,j))+a44*bs.particles.X(bs.control_points(m+1,j));
    
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
    return LOG::sprintf("B_SPLINE_PATCH<VECTOR<T,%d> ,%d>",TV::dimension,d);
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
template class B_SPLINE_PATCH<VECTOR<float,3>,3>;
template class B_SPLINE_PATCH<VECTOR<double,3>,3>;

template TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<float,3>,2>::OBJECT* Create_Triangulated_Surface<VECTOR<float,3>,3>(B_SPLINE_PATCH<VECTOR<float,3>,3> const&,bool);
template TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<double,3>,2>::OBJECT* Create_Triangulated_Surface<VECTOR<double,3>,3>(B_SPLINE_PATCH<VECTOR<double,3>,3> const&,bool);

template void Smooth_Fit<VECTOR<float,3> >(B_SPLINE_PATCH<VECTOR<float,3>,3>&,ARRAY<VECTOR<float,3>,VECTOR<int,2>> const&,bool,bool);
template void Smooth_Fit<VECTOR<double,3> >(B_SPLINE_PATCH<VECTOR<double,3>,3>&,ARRAY<VECTOR<double,3>,VECTOR<int,2>> const&,bool,bool);
}
