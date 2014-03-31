//#####################################################################
// Copyright 2014, Yuting Wang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSISTENT_INTERSECTIONS
//##################################################################### 
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include "CONSISTENT_INTERSECTIONS.h"
using namespace PhysBAM;
using std::abs;
//#####################################################################
// Function Set_Tol
//#####################################################################
template<class T> void CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Set_Tol()
{
    T L_ta=0,L_sc=0;
    for(int i=0;i<ta.mesh.elements.m;i++)
        L_ta=max(L_ta,RANGE<TV>::Bounding_Box(ta.particles.X.Subset(ta.mesh.elements(i))).Edge_Lengths().Max());
    for(int i=0;i<sc.mesh.elements.m;i++)
        L_sc=max(L_sc,RANGE<TV>::Bounding_Box(sc.particles.X.Subset(sc.mesh.elements(i))).Edge_Lengths().Max());
    T L=L_ta+L_sc;
    T eps=std::numeric_limits<T>::epsilon(),sqrt_eps=sqrt(eps);
    L*=(1+5*eps)/(1-7*sqrt_eps);
    L*=1e6;
    T sqrt_eps_L=sqrt_eps*L;

    sigma=(T)6.5*sqrt_eps_L;
    tau=(T)4.5*sqrt_eps_L;
    sigma_hat=(T)5.5*sqrt_eps_L;
    kappa=21*eps*sqr(L);
}
//#####################################################################
// Function Compute_VV
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Compute_VV(int p,int q)
{
    TV A=ta.particles.X(p),B=sc.particles.X(q);
    T d2=(A-B).Magnitude_Squared();
    if(d2>sqr(sigma)) return false;
    hash_vv.Set(I2(p,q));
    return true;
}
//#####################################################################
// Function Compute_VE_Helper
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Compute_VE_Helper(int p,I2 e,ARRAY_VIEW<TV> Xp,ARRAY_VIEW<TV> Xe,T& gamma)
{
    TV P=Xp(p),A=Xe(e.x),B=Xe(e.y),u=B-A;
    T m=u.Normalize();
    if(m<=sigma_hat) return false;
    TV w=P-A;
    T d=u.Cross(w).Magnitude();
    if(d>tau) return false;
    T a_hat=u.Dot(w),a_bar=m-a_hat;
    if(a_hat<0 || a_bar<0) return false;
    gamma=a_hat/m;
    return true;
}
//#####################################################################
// Function Compute_VE
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Compute_VE(int p,I2 e)
{
    if(hash_vv.Contains(I2(p,e.x))) return false;
    if(hash_vv.Contains(I2(p,e.y))) return false;
    T gamma=0;
    e.Sort();
    if(!Compute_VE_Helper(p,e,ta.particles.X,sc.particles.X,gamma)) return false;
    hash_ve.Set(e.Insert(p,0),gamma);
    return true;
}
//#####################################################################
// Function Compute_EV
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Compute_EV(I2 e,int p)
{
    if(hash_vv.Contains(I2(e.x,p))) return false;
    if(hash_vv.Contains(I2(e.y,p))) return false;
    T gamma=0;
    e.Sort();
    if(!Compute_VE_Helper(p,e,sc.particles.X,ta.particles.X,gamma)) return false;
    hash_ev.Set(e.Append(p),gamma);
    return true;
}
//#####################################################################
// Function Compute_EE
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Compute_EE(I2 e,I2 g)
{
    e.Sort();
    g.Sort();
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++)
            if(hash_vv.Contains(I2(e(i),g(j))))
                return false;
        if(hash_ve.Contains(g.Insert(e(i),0))) return false;
        if(hash_ev.Contains(e.Append(g(i)))) return false;}

    TV A=ta.particles.X(e.x),B=ta.particles.X(e.y);
    TV P=sc.particles.X(g.x),Q=sc.particles.X(g.y);
    T a_P=(P-A).Cross(B-A).x;
    T a_Q=(Q-A).Cross(B-A).x;
    if(abs(a_P)<=kappa || abs(a_Q)<=kappa || (a_P<0)==(a_Q<0)) return false;
    T a_A=(A-P).Cross(Q-P).x;
    T a_B=(B-P).Cross(Q-P).x;
    if(abs(a_A)<=kappa || abs(a_B)<=kappa || (a_A<0)==(a_B<0)) return false;

    T gamma_e=a_A/(a_A-a_B);
    T gamma_g=a_P/(a_P-a_Q);
    hash_ee.Set(e.Append_Elements(g),TV(gamma_e,gamma_g));
    return true;
}
//#####################################################################
// Function Compute_FV
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Compute_FV(I3 f,int p)
{
    f.Sort();
    for(int i=0;i<3;i++){
        if(hash_vv.Contains(I2(f(i),p))) return false;
        if(hash_ev.Contains(f.Remove_Index(i).Append(p))) return false;}

    TV A=ta.particles.X(f.x),B=ta.particles.X(f.y),C=ta.particles.X(f.z);
    TV P=sc.particles.X(p);
    T a_A=(B-P).Cross(C-P).x;
    T a_B=(P-A).Cross(C-A).x;
    T a_C=(B-A).Cross(P-A).x;
    if(abs(a_A)<=kappa || abs(a_B)<=kappa || abs(a_C)<=kappa) return false;
    bool a=a_A<0,b=a_B<0,c=a_C<0;
    if(a!=b || b!=c) return false;
    T3 gamma(a_A,a_B,a_C);
    hash_fv.Set(f.Append(p),gamma/gamma.Sum());
    return true;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Compute()
{
    Set_Tol();
    sc.Initialize_Hierarchy();
    ARRAY<int> sc_particles(sc.mesh.elements.Flattened());
    sc_particles.Prune_Duplicates();
    PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > sc_ph(sc.particles.X.Subset(sc_particles),true,0);

    ta.Initialize_Hierarchy();
    SEGMENTED_CURVE<TV> ta_sc(ta.Get_Segment_Mesh(),ta.particles);
    ta_sc.Initialize_Hierarchy();
    ARRAY<int> ta_particles(ta.mesh.elements.Flattened());
    ta_particles.Prune_Duplicates();
    PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > ta_ph(ta.particles.X.Subset(ta_particles),true,0);

    ARRAY<ARRAY<int> > a;
    BOX_VISITOR_TRIVIAL v(a);
    a.Resize(ta_particles.m);
    ta_ph.Intersection_List(sc_ph,v,sigma);
    for(int i=0;i<a.m;i++){
        int p=ta_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VV(p,sc_particles(a(i)(j)));}

    a.Remove_All();
    a.Resize(ta_particles.m);
    ta_ph.Intersection_List(*sc.hierarchy,v,tau);
    for(int i=0;i<a.m;i++){
        int p=ta_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VE(p,sc.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(ta_sc.mesh.elements.m);
    ta_sc.hierarchy->Intersection_List(sc_ph,v,tau);
    for(int i=0;i<a.m;i++){
        I2 e=ta_sc.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_EV(e,sc_particles(a(i)(j)));}

    a.Remove_All();
    a.Resize(ta_sc.mesh.elements.m);
    ta_sc.hierarchy->Intersection_List(*sc.hierarchy,v,0);
    for(int i=0;i<a.m;i++){
        I2 e=ta_sc.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_EE(e,sc.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(ta.mesh.elements.m);
    ta.hierarchy->Intersection_List(sc_ph,v,0);
    for(int i=0;i<a.m;i++){
        I3 f=ta.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_FV(f,sc_particles(a(i)(j)));}
}
//#####################################################################
// Function Set_Tol
//#####################################################################
template<class T> void CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Set_Tol()
{
    T L_tv=0,L_ts=0;
    for(int i=0;i<tv.mesh.elements.m;i++)
        L_tv=max(L_tv,RANGE<TV>::Bounding_Box(tv.particles.X.Subset(tv.mesh.elements(i))).Edge_Lengths().Max());
    for(int i=0;i<ts.mesh.elements.m;i++)
        L_ts=max(L_ts,RANGE<TV>::Bounding_Box(ts.particles.X.Subset(ts.mesh.elements(i))).Edge_Lengths().Max());
    T L=L_tv+L_ts;
    T eps=std::numeric_limits<T>::epsilon(),sqrt_eps=sqrt(eps),sqrt_sqrt_eps=sqrt(sqrt_eps),eps_3_4=sqrt_eps*sqrt_sqrt_eps;
    L*=(1+5*eps)/(1-7*sqrt_sqrt_eps);
    T sqrt_sqrt_eps_L=sqrt_sqrt_eps*L,L2=L*L,L3=L*L2,L4=L2*L2,eps_3_4_L3=eps_3_4*L3,eps_L4=eps*L4;

    sigma=(T)6.5*sqrt_sqrt_eps_L;
    tau=(T)4.5*sqrt_sqrt_eps_L;
    delta=(T)2.25*sqrt_sqrt_eps_L;
    gamma=(T)2.25*sqrt_sqrt_eps_L;
    sigma_hat=(T)5.5*sqrt_sqrt_eps_L;

    mu=24*eps_3_4_L3;
    rho=56*eps_3_4_L3;
    xi=56*eps_3_4_L3;

    lambda=1014*eps_L4;
    nu=(T)6844.5*eps_L4;
    zeta=1317*eps_L4;
}
//#####################################################################
// Function Compute_VV
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_VV(int p,int q)
{
    TV A=tv.particles.X(p),B=ts.particles.X(q);
    T d2=(A-B).Magnitude_Squared();
    if(d2>sqr(sigma)) return false;
    hash_vv.Set(I2(p,q));
    return true;
}
//#####################################################################
// Function Compute_VE_Helper
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_VE_Helper(int p,I2 e,ARRAY_VIEW<TV> Xp,ARRAY_VIEW<TV> Xe,T& gamma)
{
    TV P=Xp(p),A=Xe(e.x),B=Xe(e.y),u=B-A;
    T m=u.Normalize();
    if(m<=sigma_hat) return false;
    TV w=P-A;
    T d2=u.Cross(w).Magnitude_Squared();
    if(d2>sqr(tau)) return false;
    T a_hat=u.Dot(w),a_bar=m-a_hat;
    if(a_hat<0 || a_bar<0) return false;
    gamma=a_hat/m;
    return true;
}
//#####################################################################
// Function Compute_VE
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_VE(int p,I2 e)
{
    if(hash_vv.Contains(I2(p,e.x))) return false;
    if(hash_vv.Contains(I2(p,e.y))) return false;
    T gamma=0;
    e.Sort();
    if(!Compute_VE_Helper(p,e,tv.particles.X,ts.particles.X,gamma)) return false;
    hash_ve.Set(e.Insert(p,0),gamma);
    return true;
}
//#####################################################################
// Function Compute_EV
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_EV(I2 e,int p)
{
    if(hash_vv.Contains(I2(e.x,p))) return false;
    if(hash_vv.Contains(I2(e.y,p))) return false;
    T gamma=0;
    e.Sort();
    if(!Compute_VE_Helper(p,e,ts.particles.X,tv.particles.X,gamma)) return false;
    hash_ev.Set(e.Append(p),gamma);
    return true;
}
//#####################################################################
// Function Compute_EE
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_EE(I2 e,I2 g)
{
    e.Sort();
    g.Sort();
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++)
            if(hash_vv.Contains(I2(e(i),g(j))))
                return false;
        if(hash_ve.Contains(g.Insert(e(i),0))) return false;
        if(hash_ev.Contains(e.Append(g(i)))) return false;}

    TV A=tv.particles.X(e.x),B=tv.particles.X(e.y);
    TV P=ts.particles.X(g.x),Q=ts.particles.X(g.y);
    TV u=B-A,v=Q-P,r=u.Cross(v);
    T m2=r.Magnitude_Squared();
    if(m2<=lambda) return false;
    TV w=P-A;
    T d_hat=r.Dot(w);
    if(sqr(d_hat)>sqr(gamma)*m2) return false;
    TV n=r.Cross(w);
    T a_hat=n.Dot(v),b_hat=n.Dot(u),a_bar=m2-a_hat,b_bar=m2-b_hat;
    if(a_hat<=lambda || a_bar<=lambda || b_hat<=lambda || b_bar<=lambda) return false;
    hash_ee.Set(e.Append_Elements(g),T2(a_hat,b_hat)/m2);
    return true;
}
//#####################################################################
// Function Compute_VF_Helper
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_VF_Helper(int p,I3 f,ARRAY_VIEW<TV> Xp,ARRAY_VIEW<TV> Xf,TV& gamma)
{
    TV A=Xf(f.x),B=Xf(f.y),C=Xf(f.z);
    TV P=Xp(p);
    TV u=B-A,v=C-A,r=u.Cross(v);
    T m2=r.Magnitude_Squared();
    if(m2<=nu) return false;
    TV w=P-A;
    T d_hat=r.Dot(w);
    if(sqr(d_hat)>sqr(delta)*m2) return false;
    TV n=r.Cross(w);
    T b_hat=n.Dot(v),c_hat=-n.Dot(u),a_hat=m2-b_hat-c_hat;
    if(a_hat<=zeta || b_hat<=zeta || c_hat<=zeta) return false;
    gamma=TV(a_hat,b_hat,c_hat)/m2;
    return true;
}
//#####################################################################
// Function Compute_FV
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_FV(I3 f,int p)
{
    f.Sort();
    for(int i=0;i<3;i++){
        if(hash_vv.Contains(I2(f(i),p))) return false;
        if(hash_ev.Contains(f.Remove_Index(i).Append(p))) return false;}

    TV gamma;
    if(!Compute_VF_Helper(p,f,ts.particles.X,tv.particles.X,gamma)) return false;
    hash_fv.Set(f.Append(p),gamma);
    return true;
}
//#####################################################################
// Function Compute_VF
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_VF(int p,I3 f)
{
    f.Sort();
    for(int i=0;i<3;i++){
        if(hash_vv.Contains(I2(p,f(i)))) return false;
        if(hash_ve.Contains(f.Remove_Index(i).Insert(p,0))) return false;}

    TV gamma;
    if(!Compute_VF_Helper(p,f,tv.particles.X,ts.particles.X,gamma)) return false;
    hash_vf.Set(f.Insert(p,0),gamma);
    return true;
}
//#####################################################################
// Function Compute_EF_Helper
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_EF_Helper(I2 e,I3 f,ARRAY_VIEW<TV> Xe,ARRAY_VIEW<TV> Xf,T& gamma,TV& bary)
{
    TV A=Xf(f.x),B=Xf(f.y),C=Xf(f.z);
    TV P=Xe(e.x),Q=Xe(e.y);
    T v_A=(B-Q).Cross(C-Q).Dot(P-Q);
    T v_B=(C-Q).Cross(A-Q).Dot(P-Q);
    T v_C=(A-Q).Cross(B-Q).Dot(P-Q);
    if(abs(v_A)<=mu || abs(v_B)<=mu || abs(v_C)<=mu) return false;
    bool a=v_A<0,b=v_B<0,c=v_C<0;
    if(a!=b || b!=c) return false;

    T v_P=(A-P).Cross(B-P).Dot(C-P);
    T v_Q=(A-Q).Cross(B-Q).Dot(C-Q);
    if(abs(v_P)<=xi || abs(v_Q)<=xi || (v_P<0)==(v_Q<0)) return false;
    TV g(v_A,v_B,v_C);
    bary=g/g.Sum();
    gamma=v_P/(v_P-v_Q);
    return true;
}
//#####################################################################
// Function Compute_FE
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_FE(I3 f,I2 e)
{
    f.Sort();
    e.Sort();
    for(int i=0;i<3;i++){
        for(int j=0;j<2;j++){
            if(hash_vv.Contains(I2(f(i),e(j)))) return false;
            if(hash_ev.Contains(f.Remove_Index(i).Append(e(j)))) return false;}
        if(hash_ve.Contains(e.Insert(f(i),0))) return false;
        if(hash_ee.Contains(f.Remove_Index(i).Append_Elements(e))) return false;}
    for(int j=0;j<2;j++)
        if(hash_fv.Contains(f.Append(e(j)))) return false;

    T gamma=0;
    TV bary;
    if(!Compute_EF_Helper(e,f,ts.particles.X,tv.particles.X,gamma,bary)) return false;
    hash_fe.Set(f.Append_Elements(e),bary.Append(gamma));
    return true;
}
//#####################################################################
// Function Compute_EF
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_EF(I2 e,I3 f)
{
    f.Sort();
    e.Sort();
    for(int i=0;i<3;i++){
        for(int j=0;j<2;j++){
            if(hash_vv.Contains(I2(e(j),f(i)))) return false;
            if(hash_ve.Contains(f.Remove_Index(i).Insert(e(j),0))) return false;}
        if(hash_ev.Contains(e.Append(f(i)))) return false;
        if(hash_ee.Contains(e.Append_Elements(f.Remove_Index(i)))) return false;}
    for(int j=0;j<2;j++)
        if(hash_vf.Contains(f.Insert(e(j),0))) return false;

    T gamma=0;
    TV bary;
    if(!Compute_EF_Helper(e,f,tv.particles.X,ts.particles.X,gamma,bary)) return false;
    hash_ef.Set(e.Append_Elements(f),bary.Insert(gamma,0));
    return true;
}
//#####################################################################
// Function Compute_TV
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_TV(I4 t,int p)
{
    t.Sort();
    for(int i=0;i<4;i++){
        if(hash_vv.Contains(I2(t(i),p))) return false;
        if(hash_fv.Contains(t.Remove_Index(i).Append(p))) return false;}
    for(int i=0;i<4;i++)
        for(int j=i+1;j<4;j++)
            if(hash_ev.Contains(I3(t(i),t(j),p))) return false;

    TV A=tv.particles.X(t(0)),B=tv.particles.X(t(1)),C=tv.particles.X(t(2)),D=tv.particles.X(t(3));
    TV P=ts.particles.X(p);
    T v_A=(B-P).Cross(C-P).Dot(D-P);
    T v_B=(P-A).Cross(C-A).Dot(D-A);
    T v_C=(B-A).Cross(P-A).Dot(D-A);
    T v_D=(B-A).Cross(C-A).Dot(P-A);
    if(abs(v_A)<=rho || abs(v_B)<=rho || abs(v_C)<=rho || abs(v_D)<=rho) return false;
    bool a=v_A<0,b=v_B<0,c=v_C<0,d=v_D<0;
    if(a!=b || b!=c || c!=d) return false;
    T4 gamma(v_A,v_B,v_C,v_D);
    hash_tv.Set(t.Append(p),gamma/gamma.Sum());
    return true;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute()
{
    Set_Tol();
    ts.Initialize_Hierarchy();
    SEGMENTED_CURVE<TV> ts_sc(ts.Get_Segment_Mesh(),ts.particles);
    ts_sc.Initialize_Hierarchy();
    ARRAY<int> ts_particles(ts.mesh.elements.Flattened());
    ts_particles.Prune_Duplicates();
    PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > ts_ph(ts.particles.X.Subset(ts_particles),true,0);

    tv.Initialize_Hierarchy();
    tv.mesh.Initialize_Triangle_Mesh();
    TRIANGULATED_SURFACE<T> tv_ts(*tv.mesh.triangle_mesh,tv.particles);
    tv_ts.Initialize_Hierarchy();
    SEGMENTED_CURVE<TV> tv_sc(tv.Get_Segment_Mesh(),tv.particles);
    tv_sc.Initialize_Hierarchy();
    ARRAY<int> tv_particles(tv.mesh.elements.Flattened());
    tv_particles.Prune_Duplicates();
    PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > tv_ph(tv.particles.X.Subset(tv_particles),true,0);

    ARRAY<ARRAY<int> > a;
    BOX_VISITOR_TRIVIAL v(a);
    a.Resize(tv_particles.m);
    tv_ph.Intersection_List(ts_ph,v,sigma);
    for(int i=0;i<a.m;i++){
        int p=tv_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VV(p,ts_particles(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_particles.m);
    tv_ph.Intersection_List(*ts_sc.hierarchy,v,tau);
    for(int i=0;i<a.m;i++){
        int p=tv_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VE(p,ts_sc.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_particles.m);
    tv_ph.Intersection_List(*ts.hierarchy,v,delta);
    for(int i=0;i<a.m;i++){
        int p=tv_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VF(p,ts.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_sc.mesh.elements.m);
    tv_sc.hierarchy->Intersection_List(ts_ph,v,tau);
    for(int i=0;i<a.m;i++){
        I2 e=tv_sc.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_EV(e,ts_particles(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_sc.mesh.elements.m);
    tv_sc.hierarchy->Intersection_List(*ts_sc.hierarchy,v,gamma);
    for(int i=0;i<a.m;i++){
        I2 e=tv_sc.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_EE(e,ts_sc.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_sc.mesh.elements.m);
    tv_sc.hierarchy->Intersection_List(*ts.hierarchy,v,0);
    for(int i=0;i<a.m;i++){
        I2 e=tv_sc.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_EF(e,ts.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_ts.mesh.elements.m);
    tv_ts.hierarchy->Intersection_List(ts_ph,v,delta);
    for(int i=0;i<a.m;i++){
        I3 f=tv_ts.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_FV(f,ts_particles(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_ts.mesh.elements.m);
    tv_ts.hierarchy->Intersection_List(*ts_sc.hierarchy,v,0);
    for(int i=0;i<a.m;i++){
        I3 f=tv_ts.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_FE(f,ts_sc.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv.mesh.elements.m);
    tv.hierarchy->Intersection_List(ts_ph,v,0);
    for(int i=0;i<a.m;i++){
        I4 t=tv.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_TV(t,ts_particles(a(i)(j)));}
}
//#####################################################################
template class CONSISTENT_INTERSECTIONS<VECTOR<double,2> >;
template class CONSISTENT_INTERSECTIONS<VECTOR<double,3> >;
template class CONSISTENT_INTERSECTIONS<VECTOR<float,2> >;
template class CONSISTENT_INTERSECTIONS<VECTOR<float,3> >;
