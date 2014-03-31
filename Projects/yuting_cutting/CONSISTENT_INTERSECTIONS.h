//#####################################################################
// Copyright 2014, Yuting Wang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSISTENT_INTERSECTIONS
//##################################################################### 
#ifndef __CONSISTENT_INTERSECTIONS__
#define __CONSISTENT_INTERSECTIONS__

#include <Tools/Data_Structures/HASHTABLE.h>
namespace PhysBAM{

template<class TV> class CONSISTENT_INTERSECTIONS;
template<class T> class TRIANGULATED_AREA;
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class TV> class SEGMENTED_CURVE;
template<class T> class TRIANGULATED_SURFACE;

template<class T>
class CONSISTENT_INTERSECTIONS<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> I2;
    typedef VECTOR<int,3> I3;
    typedef VECTOR<int,4> I4;
    typedef VECTOR<T,3> T3;
    enum WORKAROUND {d=TV::m,assume=0,prune=1,test=2,safe=3};
public:

    TRIANGULATED_AREA<T>& ta;
    SEGMENTED_CURVE<TV>& sc;

    // ta is always before sc
    HASHTABLE<I2> hash_vv;
    HASHTABLE<I3,T> hash_ve,hash_ev;
    HASHTABLE<I4,TV> hash_ee;
    HASHTABLE<I4,T3> hash_fv;

    T sigma,tau,sigma_hat,kappa;

    CONSISTENT_INTERSECTIONS(TRIANGULATED_AREA<T>& ta_,SEGMENTED_CURVE<TV>& sc_)
        :ta(ta_),sc(sc_)
    {}

    ~CONSISTENT_INTERSECTIONS(){}

    void Set_Tol();
    bool Compute_VV(int p,int q);
    bool Compute_VE_Helper(int p,I2 e,ARRAY_VIEW<TV> Xp,ARRAY_VIEW<TV> Xe,T& gamma);
    bool Compute_VE(int p,I2 e);
    bool Compute_EV(I2 e,int p);
    bool Compute_EE(I2 e,I2 g);
    bool Compute_FV(I3 f,int p);
    void Compute();

//#####################################################################
};   

template<class T>
class CONSISTENT_INTERSECTIONS<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,2> I2;
    typedef VECTOR<int,3> I3;
    typedef VECTOR<int,4> I4;
    typedef VECTOR<int,5> I5;
    typedef VECTOR<T,2> T2;
    typedef VECTOR<T,4> T4;
    enum WORKAROUND {d=TV::m,assume=0,prune=1,test=2,safe=3};
public:

    TETRAHEDRALIZED_VOLUME<T>& tv;
    TRIANGULATED_SURFACE<T>& ts;

    // ta is always before sc
    HASHTABLE<I2> hash_vv;
    HASHTABLE<I3,T> hash_ve,hash_ev;
    HASHTABLE<I4,T2> hash_ee;
    HASHTABLE<I4,TV> hash_fv,hash_vf;
    HASHTABLE<I5,T4> hash_fe,hash_ef;
    HASHTABLE<I5,T4> hash_tv;

    T sigma,tau,delta,gamma;
    T sigma_hat,lambda,mu,nu,rho,xi,zeta;

    CONSISTENT_INTERSECTIONS(TETRAHEDRALIZED_VOLUME<T>& tv_,TRIANGULATED_SURFACE<T>& ts_)
        :tv(tv_),ts(ts_)
    {}

    ~CONSISTENT_INTERSECTIONS(){}

    void Set_Tol();
    bool Compute_VV(int p,int q);
    bool Compute_VE_Helper(int p,I2 e,ARRAY_VIEW<TV> Xp,ARRAY_VIEW<TV> Xe,T& gamma);
    bool Compute_VE(int p,I2 e);
    bool Compute_EV(I2 e,int p);
    bool Compute_EE(I2 e,I2 g);
    bool Compute_VF_Helper(int p,I3 f,ARRAY_VIEW<TV> Xp,ARRAY_VIEW<TV> Xf,TV& gamma);
    bool Compute_FV(I3 f,int p);
    bool Compute_VF(int p,I3 f);
    bool Compute_EF_Helper(I2 e,I3 f,ARRAY_VIEW<TV> Xe,ARRAY_VIEW<TV> Xf,T& gamma,TV& bary);
    bool Compute_FE(I3 f,I2 e);
    bool Compute_EF(I2 e,I3 f);
    bool Compute_TV(I4 t,int p);
    void Compute();

//#####################################################################
};   
}
#endif
