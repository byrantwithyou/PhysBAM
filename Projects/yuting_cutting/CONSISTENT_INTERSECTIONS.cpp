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
    L=1e5;
    T eps=std::numeric_limits<T>::epsilon(),sqrt_eps_L=sqrt(eps)*L;
    
    tol_vv[safe]=54*sqrt_eps_L;
    tol_vv[assume]=53*sqrt_eps_L;
    tol_ev[safe]=5*sqrt_eps_L;
    tol_ev[assume]=4*sqrt_eps_L;

    tol_vv[test]=(3*tol_vv[safe]+tol_vv[assume])/4;
    tol_vv[prune]=(tol_vv[safe]+3*tol_vv[assume])/4;
    tol_ev[test]=(3*tol_ev[safe]+tol_ev[assume])/4;
    tol_ev[prune]=(tol_ev[safe]+3*tol_ev[assume])/4;
}
//#####################################################################
// Function Compute_VV
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Compute_VV(int p,int q)
{
    T d2=(ta.particles.X(p)-sc.particles.X(q)).Magnitude_Squared();
    if(d2>sqr(tol_vv[test])) return false;
    hash_vv.Set(I2(p,q));
    return true;
}
//#####################################################################
// Function Compute_VE_Helper
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,2> >::
Compute_VE_Helper(int p,I2 e,ARRAY_VIEW<TV> Xp,ARRAY_VIEW<TV> Xe,T& gamma)
{
    TV P=Xp(p),A=Xe(e.x),B=Xe(e.y),u=B-A,w=P-A;
    T m=u.Normalize();
    if(m<tol_vv[prune]) return false;
    T d=u.Cross(w).Magnitude();
    if(d>tol_ev[test]) return false;
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

    T tol=sqr(tol_ev[prune]);
    TV A=ta.particles.X(e.x),B=ta.particles.X(e.y);
    TV P=sc.particles.X(g.x),Q=sc.particles.X(g.y);
    T area_ABP=TRIANGLE_2D<T>::Signed_Area(A,B,P);
    T area_ABQ=TRIANGLE_2D<T>::Signed_Area(A,B,Q);
    if(abs(area_ABP)<tol || abs(area_ABQ)<tol || (area_ABP<0)==(area_ABQ<0)) return false;
    T area_APQ=TRIANGLE_2D<T>::Signed_Area(A,P,Q);
    T area_BPQ=TRIANGLE_2D<T>::Signed_Area(B,P,Q);
    if(abs(area_APQ)<tol || abs(area_BPQ)<tol || (area_APQ<0)==(area_BPQ<0)) return false;

    T gamma_e=area_APQ/(area_APQ-area_BPQ);
    T gamma_g=area_ABP/(area_ABP-area_ABQ);
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

    T tol=sqr(tol_ev[prune]);
    TV A=ta.particles.X(f.x),B=ta.particles.X(f.y),C=ta.particles.X(f.z);
    TV P=sc.particles.X(p);
    T area_PBC=TRIANGLE_2D<T>::Signed_Area(P,B,C);
    T area_APC=TRIANGLE_2D<T>::Signed_Area(A,P,C);
    T area_ABP=TRIANGLE_2D<T>::Signed_Area(A,B,P);
    if(abs(area_PBC)<tol || abs(area_APC)<tol || abs(area_ABP)<tol) return false;
    bool a=area_PBC<0,b=area_APC<0,c=area_ABP<0;
    if(a!=b || b!=c) return false;
    T3 gamma(area_PBC,area_APC,area_ABP);
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
    ta_ph.Intersection_List(sc_ph,v,tol_vv[test]);
    for(int i=0;i<a.m;i++){
        int p=ta_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VV(p,sc_particles(a(i)(j)));}

    a.Remove_All();
    a.Resize(ta_particles.m);
    ta_ph.Intersection_List(*sc.hierarchy,v,tol_ev[test]);
    for(int i=0;i<a.m;i++){
        int p=ta_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VE(p,sc.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(ta_sc.mesh.elements.m);
    ta_sc.hierarchy->Intersection_List(sc_ph,v,tol_ev[test]);
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
    L*=1e-3;
    T eps=std::numeric_limits<T>::epsilon(),sqrt_sqrt_eps_L=sqrt(sqrt(eps))*L;
    std::cout << "sqrt_eps_L: " << sqrt_sqrt_eps_L << std::endl;
    
    T* tol[4]={tol_vv,tol_ev,tol_ee,tol_fv};
    T safe_scales[4]={7,5,3,3};
    T assume_scales[4]={6,4,2,2};

    for(int i=0;i<4;i++){
        T a=safe_scales[i]*sqrt_sqrt_eps_L,b=assume_scales[i]*sqrt_sqrt_eps_L;
        tol[i][safe]=a;
        tol[i][assume]=b;
        tol[i][test]=(3*a+b)/4;
        tol[i][prune]=(a+3*b)/4;}
}
//#####################################################################
// Function Compute_VV
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_VV(int p,int q)
{
    T d2=(tv.particles.X(p)-ts.particles.X(q)).Magnitude_Squared();
    if(d2>sqr(tol_vv[test])) return false;
    hash_vv.Set(I2(p,q));
    return true;
}
//#####################################################################
// Function Compute_VE_Helper
//#####################################################################
template<class T> bool CONSISTENT_INTERSECTIONS<VECTOR<T,3> >::
Compute_VE_Helper(int p,I2 e,ARRAY_VIEW<TV> Xp,ARRAY_VIEW<TV> Xe,T& gamma)
{
    TV P=Xp(p),A=Xe(e.x),B=Xe(e.y),u=B-A,w=P-A;
    T m=u.Normalize();
    if(m<tol_vv[prune]) return false;
    T d2=u.Cross(w).Magnitude_Squared();
    if(d2>sqr(tol_ev[test])) return false;
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
    TV u=B-A,v=Q-P,w=P-A,r=u.Cross(v);
    T m2=r.Magnitude_Squared();
    if(m2<16*sqr(sqr(tol_ev[prune])-sqr(tol_ee[test]))) return false;
    T d_hat=r.Dot(w);
    if(sqr(d_hat)>sqr(tol_ee[test])*m2) return false;
    TV n=r.Cross(w);
    T a_hat=n.Dot(v),b_hat=n.Dot(u),a_bar=m2-a_hat,b_bar=m2-b_hat;
    if(a_hat<0 || a_bar<0 || b_hat<0 || b_bar<0) return false;
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
    TV u=B-A,v=C-A,w=P-A,r=u.Cross(v),n=r.Cross(w);
    T m2=r.Magnitude_Squared();
    if(m2<108*sqr(sqr(tol_ev[prune])-sqr(tol_fv[test]))) return false;
    T d_hat=r.Dot(w),b_hat=n.Dot(v),c_hat=-n.Dot(u),a_hat=m2-b_hat-c_hat;
    if(sqr(d_hat)>sqr(tol_fv[test])*m2) return false;
    if(a_hat<0 || b_hat<0 || c_hat<0) return false;
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
    T vA=TETRAHEDRON<T>::Signed_Volume(B,C,P,Q);
    T vB=TETRAHEDRON<T>::Signed_Volume(C,A,P,Q);
    T vC=TETRAHEDRON<T>::Signed_Volume(A,B,P,Q);

    T tol1=2*tol_fv[prune]*sqr(tol_ee[prune]);
    if(abs(vA)<tol1 || abs(vB)<tol1 || abs(vC)<tol1) return false;
    bool a=vA<0,b=vB<0,c=vC<0;
    if(a!=b || b!=c) return false;

    T vP=TETRAHEDRON<T>::Signed_Volume(A,B,C,P);
    T vQ=TETRAHEDRON<T>::Signed_Volume(A,B,C,Q);
    T tol2=3*sqrt(3)*tol_fv[prune]*sqr(tol_ee[prune]);
    if(abs(vP)<tol2 || abs(vQ)<tol2) return false;
    if((vP<0)==(vQ<0)) return false;

    TV g(vA,vB,vC);
    bary=g/g.Sum();
    gamma=vP/(vP-vQ);
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

    T tol=sqrt((T)3)*cube(tol_fv[prune]);
    TV A=tv.particles.X(t(0)),B=tv.particles.X(t(1)),C=tv.particles.X(t(2)),D=tv.particles.X(t(3));
    TV P=ts.particles.X(p);
    T vA=TETRAHEDRON<T>::Signed_Volume(P,B,C,D);
    T vB=TETRAHEDRON<T>::Signed_Volume(A,P,C,D);
    T vC=TETRAHEDRON<T>::Signed_Volume(A,B,P,D);
    T vD=TETRAHEDRON<T>::Signed_Volume(A,B,C,P);
    if(abs(vA)<tol || abs(vB)<tol || abs(vC)<tol || abs(vD)<tol) return false;
    bool a=vA<0,b=vB<0,c=vC<0,d=vD<0;
    if(a!=b || b!=c || c!=d) return false;
    T4 gamma(vA,vB,vC,vD);
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
    tv_ph.Intersection_List(ts_ph,v,tol_vv[test]);
    for(int i=0;i<a.m;i++){
        int p=tv_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VV(p,ts_particles(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_particles.m);
    tv_ph.Intersection_List(*ts_sc.hierarchy,v,tol_ev[test]);
    for(int i=0;i<a.m;i++){
        int p=tv_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VE(p,ts_sc.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_particles.m);
    tv_ph.Intersection_List(*ts.hierarchy,v,tol_fv[test]);
    for(int i=0;i<a.m;i++){
        int p=tv_particles(i);
        for(int j=0;j<a(i).m;j++)
            Compute_VF(p,ts.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_sc.mesh.elements.m);
    tv_sc.hierarchy->Intersection_List(ts_ph,v,tol_ev[test]);
    for(int i=0;i<a.m;i++){
        I2 e=tv_sc.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_EV(e,ts_particles(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_sc.mesh.elements.m);
    tv_sc.hierarchy->Intersection_List(*ts_sc.hierarchy,v,tol_ee[test]);
    for(int i=0;i<a.m;i++){
        I2 e=tv_sc.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_EE(e,ts_sc.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_sc.mesh.elements.m);
    tv_sc.hierarchy->Intersection_List(*ts.hierarchy,v,tol_fv[test]);
    for(int i=0;i<a.m;i++){
        I2 e=tv_sc.mesh.elements(i);
        for(int j=0;j<a(i).m;j++)
            Compute_EF(e,ts.mesh.elements(a(i)(j)));}

    a.Remove_All();
    a.Resize(tv_ts.mesh.elements.m);
    tv_ts.hierarchy->Intersection_List(ts_ph,v,tol_fv[test]);
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
