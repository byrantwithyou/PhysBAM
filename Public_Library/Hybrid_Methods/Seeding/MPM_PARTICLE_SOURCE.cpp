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
    const TV& X0,const TV& n,const TV& V0,const TV& g,IMPLICIT_OBJECT<TV>* io)
    :X0(X0),n(n.Normalized()),V0(V0),g(g),io(io),poisson_disk(poisson_disk),random(random)
{
    frame.t=X0;
    frame.r=ROTATION<TV>::From_Rotated_Vector(TV::Axis_Vector(0),n);
    io->Update_Box();
    sample_box=ORIENTED_BOX<TV>(io->box,frame.Inverse()).Axis_Aligned_Bounding_Box();
    sample_box.min_corner.x=0;
    sample_box.max_corner.x=0;
    PHYSBAM_ASSERT(V0.Dot(n)>0);
}
//#####################################################################
// Function Seed
//#####################################################################
template<class TV> void MPM_PARTICLE_SOURCE<TV>::
Seed(T dt,ARRAY<TV>& X,ARRAY<TV>& V,ARRAY<MATRIX<T,TV::m> >* dV)
{
    T seed_velocity=V0.Dot(n);
    sample_box.max_corner.x=seed_velocity*dt;
    int seed_X_m=seed_X.m;
    poisson_disk.Sample(random,sample_box,seed_X);

    for(int i=seed_X_m;i<seed_X.m;i++){
        TV Z=seed_X(i);
        T ts=Z.x/seed_velocity;
        Z.x=0;
        Z=frame*Z;
        if(!io->Lazy_Inside(Z)) continue;
        T t=dt-ts;
        X.Append(Z+t*V0+t*t/2*g);
        V.Append(V0+t*g);
        if(dV) dV->Append(Outer_Product(g,n/(V0.Dot(n)+t*g.Dot(n))));}
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
        bool inside=true;
        for(int j=1;j<TV::m;j++){
            if(Y(j)<sample_box.min_corner(j) || Y(j)>sample_box.max_corner(j)){
                inside=false;
                break;}}
        if(Y.x>=-2*poisson_disk.h && inside)
            seed_X.Append(Y);}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void MPM_PARTICLE_SOURCE<TV>::
Write(TYPED_OSTREAM output) const
{
    Write_Binary(output,seed_X);
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void MPM_PARTICLE_SOURCE<TV>::
Read(TYPED_ISTREAM input)
{
    Read_Binary(input,seed_X);
}
namespace PhysBAM{
template class MPM_PARTICLE_SOURCE<VECTOR<float,3> >;
template class MPM_PARTICLE_SOURCE<VECTOR<float,2> >;
template class MPM_PARTICLE_SOURCE<VECTOR<double,3> >;
template class MPM_PARTICLE_SOURCE<VECTOR<double,2> >;
}
