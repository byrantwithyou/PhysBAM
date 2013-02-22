//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include "SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV>::
SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV>::
~SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL()
{}
//#####################################################################
// Function Compute_Kernal_Centers_And_Transformation
//#####################################################################
template<class TV> void SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV>::
Compute_Kernal_Centers_And_Transformation(const GEOMETRY_PARTICLES<TV>& particles,const T h,const T r,const T lambda,const int N_eps,const T kr,const T ks,const T kn,ARRAY<TV>& Xbar,ARRAY<MATRIX<T,TV::m> >& G) const
{
    LOG::cout<<"1 "<<r<<std::endl;
    ARRAY<T> total_weight_on_particle(particles.number);
    ARRAY<TV> xw(particles.X);
    Xbar.Resize(particles.number);

    HASHTABLE<TV_INT,ARRAY<int> > buckets;
    for(int i=0;i<particles.number;i++)
        buckets.Get_Or_Insert(TV_INT(floor(particles.X(i)/r))).Append(i);

    LOG::cout<<"a"<<std::endl;
    T r2=r*r;
    ARRAY<TRIPLE<int,int,T> > pairs;
    for(typename HASHTABLE<TV_INT,ARRAY<int> >::ITERATOR it(buckets);it.Valid();it.Next()){
        const ARRAY<int>& l0=it.Data();
        for(RANGE_ITERATOR<TV::m> it2(RANGE<TV_INT>(it.Key()-1,it.Key()+2));it2.Valid();it2.Next()){
            if(const ARRAY<int>* l1=buckets.Get_Pointer(it2.index)){
                for(int i=0;i<l0.m;i++)
                    for(int j=0;j<l1->m;j++)
                        if(l0(i)<(*l1)(j)){
                            TRIPLE<int,int,T> tr(l0(i),(*l1)(j),(particles.X(l0(i))-particles.X((*l1)(j))).Magnitude_Squared());
                            if(tr.z<r2){
                                tr.z=1-cube(sqrt(tr.z)/r);
                                pairs.Append(tr);}}}}}
    LOG::cout<<"b"<<std::endl;

    for(int i=0;i<pairs.m;i++){
        total_weight_on_particle(pairs(i).x)+=pairs(i).z;
        total_weight_on_particle(pairs(i).y)+=pairs(i).z;
        xw(pairs(i).y)+=pairs(i).z*particles.X(pairs(i).x);
        xw(pairs(i).x)+=pairs(i).z*particles.X(pairs(i).y);}
    for(int i=0;i<particles.number;i++){
        xw(i)/=total_weight_on_particle(i)+1;
        Xbar(i)=((T)1-lambda)*particles.X(i)+lambda*xw(i);}
    LOG::cout<<"2"<<std::endl;
    ARRAY<int> neighbor_count(particles.number);
    ARRAY<MATRIX<T,TV::m> > C(particles.number);
    for(int i=0;i<particles.number;i++) C(i)=MATRIX<T,TV::m>::Outer_Product(particles.X(i)-xw(i),particles.X(i)-xw(i));
    for(int i=0;i<pairs.m;i++){
        neighbor_count(pairs(i).x)++;
        neighbor_count(pairs(i).y)++;
        TV d1=particles.X(pairs(i).x)-xw(pairs(i).y),d2=particles.X(pairs(i).y)-xw(pairs(i).x);
        C(pairs(i).x)+=MATRIX<T,TV::m>::Outer_Product(d2,pairs(i).z*d2);
        C(pairs(i).y)+=MATRIX<T,TV::m>::Outer_Product(d1,pairs(i).z*d1);}
    for(int i=0;i<particles.number;i++) C(i)/=total_weight_on_particle(i)+1;
    LOG::cout<<"3"<<std::endl;
    G.Resize(particles.number);
    T one_over_h=(T)1/h;
    for(int i=0;i<particles.number;i++){
        DIAGONAL_MATRIX<T,3> eigenvalues,Sigma_wave;
        MATRIX<T,3> eigenvectors;
        SYMMETRIC_MATRIX<T,3>(C(i).Symmetric_Part()).Solve_Eigenproblem(eigenvalues,eigenvectors);
        int index_max=eigenvalues.To_Vector().Arg_Max();
        T max_ev=eigenvalues(index_max,index_max);
        if(neighbor_count(i)<=N_eps) Sigma_wave+=kn;
        else{
            if(kr<1){Sigma_wave+=ks*max_ev/kr;Sigma_wave(index_max,index_max)=ks*max_ev;}
            else for(int d=0;d<TV::m;d++) Sigma_wave(d,d)=ks*max(eigenvalues(d,d),max_ev/kr);}

        G(i)=one_over_h*(eigenvectors*Sigma_wave.Inverse()).Times_Transpose(eigenvectors);}
}
//#####################################################################
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<float,3> >;
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<double,3> >;
}
