//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/MATRIX.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include <Hybrid_Methods/Seeding/MPM_PARTICLE_SOURCE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PARTICLE_SOURCE<TV>::
MPM_PARTICLE_SOURCE(POISSON_DISK<TV>& poisson_disk,RANDOM_NUMBERS<T>& random,
    const TV& X0,const TV& n,IMPLICIT_OBJECT<TV>* io,
    std::function<void(TV X,T ts,T t,SOURCE_PATH<TV>& p)> path)
    :X0(X0),n(n.Normalized()),io(io),poisson_disk(poisson_disk),random(random),
    path(path)
{
    frame.t=X0;
    frame.r=ROTATION<TV>::From_Rotated_Vector(TV::Axis_Vector(0),n);
    io->Update_Box();
    sample_box=ORIENTED_BOX<TV>(io->box,frame.Inverse()).Axis_Aligned_Bounding_Box();
    sample_box.min_corner.x=0;
    sample_box.max_corner.x=0;
}
//#####################################################################
// Function Seed
//#####################################################################
template<class TV> void MPM_PARTICLE_SOURCE<TV>::
Seed(T t0,T dt,ARRAY<TV>& X,ARRAY<TV>& V,ARRAY<MATRIX<T,TV::m> >* dV)
{
    SOURCE_PATH<TV> init_p;
    path(X0,t0,t0,init_p);
    T seed_velocity=init_p.vp.Dot(n);
    if(seed_velocity<=0) return;
    sample_box.max_corner.x=seed_velocity*dt;
    int seed_X_m=seed_X.m;
    poisson_disk.Sample(random,sample_box,seed_X);

    for(int i=seed_X_m;i<seed_X.m;i++){
        TV Z=seed_X(i);
        T ts=Z.x/seed_velocity;
        Z.x=0;
        Z=frame*Z;
        if(!io->Lazy_Inside(Z)) continue;
        SOURCE_PATH<TV> p;
        path(Z,t0+ts,t0+dt,p);
        X.Append(p.xp);
        V.Append(p.vp);
        if(dV){
            MATRIX<T,TV::m> dvp_dz=p.vp_x+MATRIX<T,TV::m>::Outer_Product(p.vp_ts/seed_velocity-p.vp_x*n,n);
            MATRIX<T,TV::m> dxp_dz=p.xp_x+MATRIX<T,TV::m>::Outer_Product(p.xp_ts/seed_velocity-p.xp_x*n,n);
            dV->Append(dvp_dz*dxp_dz.Inverse());}}
    for(int i=seed_X.m-1;i>=0;i--){
        seed_X(i).x-=sample_box.max_corner.x;
        if(seed_X(i).x<-2*poisson_disk.h)
            seed_X.Remove_Index_Lazy(i);}
}
//#####################################################################
// Function Seed_Points
//#####################################################################
template<class TV> void MPM_PARTICLE_SOURCE<TV>::
Seed_Points(ARRAY_VIEW<const TV> X)
{
    for(int i=0;i<X.m;i++){
        TV Y=frame.Inverse_Times(X(i));
        Y.x=-Y.x;
        if(Y.x>=-2*poisson_disk.h)
            seed_X.Append(Y);}
}
namespace PhysBAM{
template class MPM_PARTICLE_SOURCE<VECTOR<float,3> >;
template class MPM_PARTICLE_SOURCE<VECTOR<float,2> >;
template class MPM_PARTICLE_SOURCE<VECTOR<double,3> >;
template class MPM_PARTICLE_SOURCE<VECTOR<double,2> >;
}
