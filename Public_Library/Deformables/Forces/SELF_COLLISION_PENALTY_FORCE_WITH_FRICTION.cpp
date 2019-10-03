//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Math_Tools/cyclic_shift.h>
#include <Core/Matrices/IDENTITY_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/ZERO_MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/RELAX_ATTACHMENT_MESH.h>
#include <Deformables/Forces/SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input)
    :BASE(particles_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
~SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION()
{
}
namespace{

template<class TV,class CP,class TV_INT>
void Evalute_Derivatives(TV& dY,TV& dw,const CP& c,const TV& dZ,ARRAY_VIEW<const TV> dX,const ARRAY<TV_INT>& elements)
{
    dY=dX.Subset(elements(c.e0)).Weighted_Sum(c.w0);
    TV dYp;
    for(int i=0;i<c.ram.diff_entry.m;i++){
        dYp=dY;
        const auto& de=c.ram.diff_entry(i);
        dY=de.dYdI(0)*dY+de.dYdI(1)*dZ;
        for(int j=0;j<TV::m;j++)
            dY+=de.dYdI(2+j)*dX(elements(de.e)(j));}

    if(c.ram.diff_entry.m>0){
        dw=c.ram.dwdI(0)*dYp+c.ram.dwdI(1)*dZ;
        for(int i=0;i<TV::m;i++)
            dw+=c.ram.dwdI(2+i)*dX(elements(c.ram.diff_entry.Last().e)(i));}
    else dw=TV();
}

template<class TV,class CP,class TV_INT>
void Apply_Derivatives_Transpose(TV dY,const TV& dw,const CP& c,TV& dZ,ARRAY_VIEW<TV> F,const ARRAY<TV_INT>& elements)
{
    TV dYp;
    if(c.ram.diff_entry.m>0){
        for(int i=0;i<TV::m;i++){
            F(elements(c.ram.diff_entry.Last().e)(i))+=c.ram.dwdI(2+i).Transpose_Times(dw);}
        dYp=c.ram.dwdI(0).Transpose_Times(dw);
        dZ=c.ram.dwdI(1).Transpose_Times(dw);}
    else dZ=TV();

    for(int i=c.ram.diff_entry.m-1;i>=0;i--){
        const auto& de=c.ram.diff_entry(i);
        for(int j=0;j<TV::m;j++)
            F(elements(de.e)(j))+=de.dYdI(2+j).Transpose_Times(dY);
        dZ+=de.dYdI(1).Transpose_Times(dY);
        dY=de.dYdI(0).Transpose_Times(dY)+dYp;
        dYp=TV();}

    for(int i=0;i<TV::m;i++)
        F(elements(c.e0)(i))+=c.w0(i)*dY;
}

template<class T,class TV,class CP>
void Visualize(const CP& c,const TV& Z,const TRIANGULATED_SURFACE<T>& ts)
{
    Dump_Surface(ts,TV(.5,.5,.5));
    TV Y=c.Y0;
    Add_Debug_Particle(Y,TV(1,0,0));
    Add_Debug_Particle(Z,TV(1,0,0));
    for(int i=0;i<c.ram.diff_entry.m;i++){
        const auto& de=c.ram.diff_entry(i);
        Add_Debug_Object(VECTOR<TV,2>(Y,de.Y),TV(1,0,1));
        Add_Debug_Particle(de.Y,TV(1,1,0));
        Y=de.Y;}
}

template<class TV,class CP,class SURF>
void Test_Weights(const CP& c,const SURF& ts)
{
    TV Y0=c.Y0,Y1=c.Y0;
    if(c.ram.diff_entry.m>0){
        const auto& le=c.ram.diff_entry.Last();
        Y1=ts.particles.X.Subset(ts.mesh.elements(le.e)).Weighted_Sum(c.ram.w);
        Y0=le.Y;}
    LOG::printf("Y vs w: %P %P   %P\n",Y0,Y1,(Y0-Y1).Magnitude());
}

template<class TV,class TV_INT,class CP,class T>
void Test_Diff(const CP& c0,const CP& c1,
    const TV& dZ,const ARRAY<TV>& dX,const ARRAY<TV_INT>& elements,T eps)
{
    if(c0.ram.diff_entry.m!=c1.ram.diff_entry.m){
        LOG::printf("SKIP Test_Diff; paths differ: len %i != %i\n",
            c0.ram.diff_entry.m,c1.ram.diff_entry.m);
        return;}
    for(int i=0;i<c0.ram.diff_entry.m;i++)
        if(c0.ram.diff_entry(i).e!=c1.ram.diff_entry(i).e){
            LOG::printf("SKIP Test_Diff; paths differ: tri(%i) %i != %i\n",
                i,c0.ram.diff_entry(i).e,c1.ram.diff_entry(i).e);
            return;}

    TV dY0,dY1,dw0,dw1;
    Evalute_Derivatives<TV>(dY0,dw0,c0,dZ,dX,elements);
    Evalute_Derivatives<TV>(dY1,dw1,c1,dZ,dX,elements);
    TV da=(c1.ram.Y-c0.ram.Y)/eps;
    TV db=(dY0+dY1)/(2*eps);
    LOG::printf("diff test: %.16P %.16P %P\n",da.Magnitude(),db.Magnitude(),(da-db).Magnitude());

    TV dc=(c1.ram.w-c0.ram.w)/eps;
    TV dd=(dw0+dw1)/(2*eps);
    LOG::printf("diff weight: %.16P %.16P %P\n",dc.Magnitude(),dd.Magnitude(),(dc-dd).Magnitude());
}
}
//#####################################################################
// Function Test_Relax
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Test_Relax(int cp)
{
    T eps=1e-6;
    
    COLLISION_PAIR c1=collision_pairs(cp),c2=collision_pairs(cp);
    auto& ts=*surfaces(c1.s);

    RANDOM_NUMBERS<T> rand;
    ARRAY<TV> dX(particles.X.m),X0(particles.X);
    rand.Fill_Uniform(dX,-eps,eps);

    TV Z=particles.X(c1.p);
    c1.ram.Relax(c1.e0,c1.w0,Z,ts,c1.p,friction);

    Test_Weights<TV>(c1,ts);
    LOG::printf("ret %i\n",c1.ram.active);
    particles.X+=dX;
    TV Z2=particles.X(c2.p),dZ=Z2-Z;
    c2.ram.Relax(c2.e0,c2.w0,Z2,ts,c2.p,friction);
    particles.X=X0;
    LOG::printf("ret %i\n",c2.ram.active);
    Test_Diff(c1,c2,dZ,dX,ts.mesh.elements,eps);
    TV Y=c1.ram.Y;
    // if(c1.active!=(torus.Normal(Y).Dot(Z-Y)<=0)){
    //     LOG::printf("INSIDE TEST DOES NOT MATCH %i %g\n",c1.active,torus.Signed_Distance(Z));
    //     Visualize(c1,Z,ts);
    // }
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.ram.active){
            const ARRAY<TV_INT>& elements=surfaces(c.s)->mesh.elements;
            TV_INT e=elements(c.ram.e);
            TV j=stiffness_coefficient*(particles.X(c.p)-c.ram.Y);
            F(c.p)-=j;
            for(int k=0;k<TV::m;k++)
                F(e(k))+=c.ram.w(k)*j;}}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    get_candidates();

    for(int i=0;i<collision_pairs.m;i++)
        Relax_Attachment(i);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.ram.active){
            const ARRAY<TV_INT>& elements=surfaces(c.s)->mesh.elements;
            TV_INT e=elements(c.ram.e);
            TV j=stiffness_coefficient*(particles.X(c.p)-c.ram.Y);
            if(transpose){
                TV dj=-V(c.p),dw,dZ;
                for(int k=0;k<TV::m;k++){
                    dj+=V(e(k))*c.ram.w(k);
                    dw(k)+=j.Dot(V(e(k)));}
                TV dY=-stiffness_coefficient*dj;
                F(c.p)-=dY;
                Apply_Derivatives_Transpose(dY,dw,c,dZ,F,elements);
                F(c.p)+=dZ;}
            else{
                TV dY,dw,dZ=V(c.p);
                Evalute_Derivatives(dY,dw,c,dZ,V,elements);
                TV dj=stiffness_coefficient*(V(c.p)-dY);
                F(c.p)-=dj;
                for(int k=0;k<TV::m;k++)
                    F(e(k))+=c.ram.w(k)*dj+dw(k)*j;}}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Potential_Energy(const T time) const
{
    T pe=0;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.ram.active)
            pe+=(T).5*stiffness_coefficient*(particles.X(c.p)-c.ram.Y).Magnitude_Squared();}
    return pe;
}
//#####################################################################
// Function Relax_Attachment
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Relax_Attachment(int cp)
{
    COLLISION_PAIR& c=collision_pairs(cp);
    TV Z=particles.X(c.p);
    c.ram.Relax(c.e0,c.w0,Z,*surfaces(c.s),c.p,friction);
    if(surfaces(c.s)->mesh.elements(c.ram.e).Contains(c.p))
        c.ram.active=false;
}
//#####################################################################
// Function Update_Attachments_And_Prune_Pairs
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Attachments_And_Prune_Pairs()
{
    int k=0;
    for(int i=0;i<collision_pairs.m;i++){
        COLLISION_PAIR c=collision_pairs(i);
        if(c.ram.active){
            if(c.ram.diff_entry.m>0){
                c.w0=c.ram.w;
                c.e0=c.ram.diff_entry.Last().e;
                c.Y0=c.ram.diff_entry.Last().Y;}
            collision_pairs(k++)=c;}
        else hash.Delete({c.p,c.s});}
    collision_pairs.Resize(k);
}
//#####################################################################
// Function Add_Pair
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Pair(int p,int s,const TV& w0,int e0)
{
    // TODO: Interpolate X^n and X^(n+1) to choose surface point.
    if(hash.Contains({p,s})) return;
    const auto& ts=*surfaces(s);
    if(ts.mesh.elements(e0).Contains(p)) return;
    TV Y=ts.particles.X.Subset(ts.mesh.elements(e0)).Weighted_Sum(w0);
    COLLISION_PAIR c={p,s,w0,e0,Y};
    collision_pairs.Append(c);
    hash.Insert({p,s});
}
//#####################################################################
// Function Add_Surface
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Surface(T_SURFACE& surface)
{
    surface.mesh.Initialize_Incident_Elements();
    surface.mesh.Initialize_Neighbor_Nodes();
    surface.mesh.Initialize_Adjacent_Elements();
    int si=surfaces.Append(&surface);
    for(int i=0;i<surface.mesh.elements.m;i++)
        object_from_element.Set(surface.mesh.elements(i).Sorted(),{si,i});
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
CFL_Strain_Rate() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Read(TYPED_ISTREAM input)
{
    ARRAY<PAIR<int,int> > keys;
    ARRAY<TV_INT> faces;
    Read_Binary(input,collision_pairs,keys,faces);
    hash.Set_All(keys);
    for(int k=0;k<faces.m;k++)
        Read_Binary(input,object_from_element.Get_Or_Insert(faces(k)));
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Write(TYPED_OSTREAM output) const
{
    ARRAY<PAIR<int,int> > keys;
    hash.Get_Keys(keys);
    keys.Sort();
    ARRAY<TV_INT> faces;
    object_from_element.Get_Keys(faces);
    faces.Sort(LEXICOGRAPHIC_COMPARE());
    Write_Binary(output,collision_pairs,keys,faces);
    for(int k=0;k<faces.m;k++)
        Write_Binary(output,object_from_element.Get(faces(k)));
}
namespace PhysBAM{
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,2> >;
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,3> >;
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,2> >;
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,3> >;
}
