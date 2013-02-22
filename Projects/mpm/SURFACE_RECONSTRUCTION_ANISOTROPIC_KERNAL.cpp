//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_3D_SPATIAL_PARTITION.h>
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
template<class TV> static void SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV>::
Compute_Kernal_Centers_And_Transformation(const GEOMETRY_PARTICLES<TV>& particles,const T h,const T r,const T lambda,const int N_eps,const T kr,const T ks,const T kn,ARRAY<TV>& Xbar,ARRAY<MATRIX<T,TV::m> >& G)
{
    ARRAY<T> total_weight_on_particle(particles.number);total_weight_on_particle.Fill(T(1));
    ARRAY<TV> xw(particles.number);for(int i=0;i<particles.number;i++) xw(i)=particles.X(i);
    Xbar.Resize(particles.number);
    PARTICLE_3D_SPATIAL_PARTITION<T> partition_spatial_partition(particles,r);
    particle_spatial_partition.Reinitialize();particle_spatial_partition.Reset_Pair_Finder();
    int index1;ARRAY<int> nearby_particle_indices;nearby_particle_indices.Preallocate(100);
    while(particle_spatial_partition.Get_Next_Particles_Potentially_Within_Interaction_Radius(index1,nearby_particle_indices)){
        TV position1=particles.X(index1);
        for(int k=0;k<nearby_particle_indices.m;k++){
            int index2=nearby_particle_indices(k);
            if(index1==index2) PHYSBAM_FATAL_ERROR;
            TV position2=particles.X(index2);
            T distance=(position2-position1).Magnitude();
            T w=(distance<r)?((T)1-cube(distance/r)):(T)0;
            total_weight_on_particle(index1)+=w;
            total_weight_on_particle(index2)+=w;
            xw(index1)+=w*position2;
            xw(index2)+=w*position1;}}
    for(int i=0;i<particles.number;i++){
        xw(i)/=total_weight_on_particle(i);
        Xbar(i)=((T)1-lambda)*partcles.X(i)+lambda*xw(i);}

    ARRAY<int> neighbor_count(particles.number);neighbor_count.Fill(0);
    ARRAY<MATRIX<T,TV::m> > C(particles.number);
    for(int i=0;i<particles.number;i++) C(i)=MATRIX<T,TV::m>::Outer_Product(particles.X(i)-xw(i),particles.X(i)-xw(i));
    while(particle_spatial_partition.Get_Next_Particles_Potentially_Within_Interaction_Radius(index1,nearby_particle_indices)){
        TV position1=particles.X(index1);
        for(int k=0;k<nearby_particle_indices.m;k++){
            int index2=nearby_particle_indices(k);
            if(index1==index2) PHYSBAM_FATAL_ERROR;
            TV position2=particles.X(index2);
            T distance=(position2-position1).Magnitude();
            T w=(distance<r)?((T)1-cube(distance/r)):(T)0;
            neighbor_count(index1)++;
            neighbor_count(index2)++;
            C(index1)+=w*MATRIX<T,TV::m>::Outer_Product(position2-xw(index1),position2-xw(index1));
            C(index2)+=w*MATRIX<T,TV::m>::Outer_Product(position1-xw(index2),position1-xw(index2));}}
    for(int i=0;i<particles.number;i++) C(i)/=total_weight_on_particle(i);

    G.resize(particles.number);
    T one_over_h=(T)1/h;
    for(int i=0;i<particles.number;i++){
        MATRIX<T,TV::m> U,V;
        DIAGONAL_MATRIX<T,TV::m> Sigma,Sigma_wave;
        C(i).Fast_Singular_Value_Decompositon(U,Sigma,V);
        for(int d=0;d<TV::m-1;d++) if(Sigma(d,d)<Sigma(d+1,d+1)) PHYSBAM_FATAL_ERROR;
        if(Sigma(TV::m-1,TV::m-1)<0) PHYSBAM_FATAL_ERROR;
        if((U-V).Frobenius_Norm_Squared()>1e-10) PHYSBAM_FATEL_ERROR;
        if(neighbor_count(i)<=N_eps) Sigma_wave=DIAGONAL_MATRIX<T,TV::m>::Identity_Matrix()*kn;
        else{
            Sigma_wave(0,0)=ks*Sigma(0,0);
            for(int d=1;d<TV::m;d++) Sigma_wave(d,d)=ks*max(Sigma(d,d),Sigma(0,0)/kr);}
        G(i)=one_over_h*U*Sigma_wave.Inverse()*V.Transposed();}
}
//#####################################################################
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<float,3> >;
template class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<VECTOR<double,3> >;
}
