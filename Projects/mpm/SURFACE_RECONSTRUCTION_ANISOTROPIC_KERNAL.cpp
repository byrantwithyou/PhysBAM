//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
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
// Function - kernel stuff
//#####################################################################
template<class T> T
P(const T x)
{
    if(x<(T)1) return (T)0.25*cube((T)2-x)-cube((T)1-x);
    else if(x<(T)2) return (T)0.25*cube((T)2-x);
    else return (T)0;
}
static double kernel_sigma[3]={0.666666666666667,0.4547284088339866859,0.3183098861837906912};
//#####################################################################
// Function Compute_Kernal_Centers_And_Transformation_And_Density
//#####################################################################
template<class TV> void SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV>::
Compute_Kernal_Centers_And_Transformation_And_Density(const ARRAY_VIEW<TV>& X,const ARRAY_VIEW<T>& m,const T h,const T r,const T lambda,const int N_eps,const T kr,const T ks,const T kn,ARRAY<TV>& Xbar,ARRAY<MATRIX<T,TV::m> >& G,ARRAY<T>& density) const
{
    // spatial hash
    HASHTABLE<TV_INT,ARRAY<int> > buckets;
    for(int i=0;i<X.m;i++)
        buckets.Get_Or_Insert(TV_INT(floor(X(i)/r))).Append(i);
    // close pairs
    T r2=r*r;
    RANGE<TV_INT> range(TV_INT()-1,TV_INT()+2);
    ARRAY<TRIPLE<int,int,T> > pairs;
    for(typename HASHTABLE<TV_INT,ARRAY<int> >::ITERATOR it(buckets);it.Valid();it.Next()){
        const ARRAY<int>& l0=it.Data();
        TV_INT key=it.Key();
        for(RANGE_ITERATOR<TV::m> it2(range);it2.Valid();it2.Next()){
            if(const ARRAY<int>* l1=buckets.Get_Pointer(it2.index+key)){
                for(int i=0;i<l0.m;i++)
                    for(int j=0;j<l1->m;j++)
                        if(l0(i)<(*l1)(j)){
                            TRIPLE<int,int,T> tr(l0(i),(*l1)(j),(X(l0(i))-X((*l1)(j))).Magnitude_Squared());
                            if(tr.z<r2){
                                tr.z=1-cube(sqrt(tr.z)/r);
                                pairs.Append(tr);}}}}}
    // Xbar and density
    ARRAY<T> total_weight_on_particle(X.m);
    ARRAY<TV> xw(X);
    Xbar.Resize(X.m);
    density.Resize(X.m);
    T one_over_h=Inverse(h);
    T one_over_h_d=one_over_h;for(int d=1;d<TV::m;d++) one_over_h_d*=one_over_h;
    for(int i=0;i<X.m;i++) density(i)=m(i)*(T)kernel_sigma[TV::m-1]*one_over_h_d;
    for(int i=0;i<pairs.m;i++){
        total_weight_on_particle(pairs(i).x)+=pairs(i).z;
        total_weight_on_particle(pairs(i).y)+=pairs(i).z;
        xw(pairs(i).y)+=pairs(i).z*X(pairs(i).x);
        xw(pairs(i).x)+=pairs(i).z*X(pairs(i).y);
        T r_over_h=(X(pairs(i).x)-X(pairs(i).y)).Magnitude()*one_over_h;
        T kernel_evaluation=(T)kernel_sigma[TV::m-1]*one_over_h_d*P(r_over_h);
        density(pairs(i).x)+=m(pairs(i).y)*kernel_evaluation;
        density(pairs(i).y)+=m(pairs(i).x)*kernel_evaluation;}
    for(int i=0;i<X.m;i++){
        xw(i)/=total_weight_on_particle(i)+1;
        Xbar(i)=((T)1-lambda)*X(i)+lambda*xw(i);}
    // C
    ARRAY<int> neighbor_count(X.m);
    ARRAY<MATRIX<T,TV::m> > C(X.m);
    for(int i=0;i<X.m;i++){
        TV dd=X(i)-xw(i);
        C(i)=MATRIX<T,TV::m>::Outer_Product(dd,dd);}
    for(int i=0;i<pairs.m;i++){
        neighbor_count(pairs(i).x)++;
        neighbor_count(pairs(i).y)++;
        TV d1=X(pairs(i).x)-xw(pairs(i).y),d2=X(pairs(i).y)-xw(pairs(i).x);
        C(pairs(i).x)+=MATRIX<T,TV::m>::Outer_Product(d2,pairs(i).z*d2);
        C(pairs(i).y)+=MATRIX<T,TV::m>::Outer_Product(d1,pairs(i).z*d1);}
    for(int i=0;i<X.m;i++) C(i)/=total_weight_on_particle(i)+1;
    // G
    G.Resize(X.m);
    for(int i=0;i<X.m;i++){
        DIAGONAL_MATRIX<T,TV::m> eigenvalues,Sigma_wave;
        MATRIX<T,TV::m> eigenvectors;
        SYMMETRIC_MATRIX<T,TV::m>(C(i).Symmetric_Part()).Solve_Eigenproblem(eigenvalues,eigenvectors);
        int index_max=eigenvalues.To_Vector().Arg_Max();
        T max_ev=eigenvalues(index_max,index_max);
        if(neighbor_count(i)<=N_eps) Sigma_wave+=kn;
        else for(int d=0;d<TV::m;d++) Sigma_wave(d,d)=ks*max(eigenvalues(d,d),max_ev/kr);
        G(i)=one_over_h*(eigenvectors*Sigma_wave.Inverse()).Times_Transpose(eigenvectors);
        // G(i)=one_over_h*MATRIX<T,TV::m>::Identity_Matrix();
    }
}
//#####################################################################
// Function Build_Scalar_Field
//#####################################################################
template<class TV> void SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV>::
Build_Scalar_Field(const ARRAY<TV>& Xbar,const ARRAY_VIEW<T>& m,const ARRAY_VIEW<T>& density,const ARRAY<MATRIX<T,TV::m> >& G,const GRID<TV>& grid,ARRAY<T,TV_INT>& phi) const
{
    phi.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));phi.Fill(T(0));
    
    for(int i=0;i<Xbar.m;i++){
        TV_INT particle_cell=grid.Cell(Xbar(i),0);
        int layer=1;
        while(1){
            bool should_go_to_next_layer=false;
            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(particle_cell-layer,particle_cell+layer+1));it.Valid();it.Next()){
                bool boundary=false;
                for(int d=0;d<TV::m;d++){
                    if(it.index(d)==(particle_cell-layer)(d) || it.index(d)==(particle_cell+layer+1)(d)){
                        boundary=true;
                        break;}}
                if(boundary && (G(i)*(grid.Node(it.index)-Xbar(i))).Magnitude_Squared()<(T)4)
                    should_go_to_next_layer=true;}
            if(should_go_to_next_layer) layer++;
            else break;}
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(particle_cell-layer,particle_cell+layer+1));it.Valid();it.Next())
            phi(it.index)+=m(i)/density(i)*kernel_sigma[TV::m-1]*G(i).Determinant()*P((G(i)*(grid.Node(it.index)-Xbar(i))).Magnitude());}
}
//#####################################################################
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<float,2> >;
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<double,2> >;
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<float,3> >;
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<double,3> >;
}
