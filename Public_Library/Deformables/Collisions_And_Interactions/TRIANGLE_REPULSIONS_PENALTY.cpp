//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_PENALTY.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int bp=0;bp<bad_pairs.m;bp++){
        int i=bad_pairs(bp);
        for(int j=0;j<TV::m+1;j++)
            F(interaction_pairs(i).nodes(j))-=grad_pe(bp)(j);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Lagged_Update_Position_Based_State(const T time)
{
    volume.Remove_All();
    for(int i=0;i<interaction_pairs.m;i++){
        TV u=particles.X(interaction_pairs(i).nodes(1))-particles.X(interaction_pairs(i).nodes(0));
        TV v=particles.X(interaction_pairs(i).nodes(2))-particles.X(interaction_pairs(i).nodes(0));
        TV w=particles.X(interaction_pairs(i).nodes(3))-particles.X(interaction_pairs(i).nodes(0));
        T vol = u.Dot(v.Cross(w));
        volume.Append((vol>=0?1:-1)*maxabs(vol,(T)1e-10));}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const
{
    for(int bp=0;bp<bad_pairs.m;bp++){
        int i=bad_pairs(bp);
        const VECTOR<int,TV::m+1>& n=interaction_pairs(i).nodes;
        for(int j=0;j<n.m;j++){
            int p=n(j);
            for(int k=0;k<n.m;k++)
                F(p)-=H_pe(bp)(j)(k)*V(n(k))*scale;}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR TRIANGLE_REPULSIONS_PENALTY<TV>::
Potential_Energy(const T time) const
{
    return pe;
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    bad_pairs.Remove_All();
    pe=0;
    grad_pe.Remove_All();
    H_pe.Remove_All();
    for(int i=0;i<interaction_pairs.m;i++){
        const VECTOR<int,TV::m+1>& n=interaction_pairs(i).nodes;
        TV u=particles.X(n(1))-particles.X(n(0));
        TV v=particles.X(n(2))-particles.X(n(0));
        TV w=particles.X(n(3))-particles.X(n(0));
        if(1<u.Dot(v.Cross(w)/volume(i)))
            continue;
        bad_pairs.Append(i);
        T e;
        VECTOR<TV,TV::m+1> de;
        VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1> he;
        Penalty(volume(i),n,particles.X,e,de,he);
        pe+=e;
        grad_pe.Append(de);
        H_pe.Append(he);}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR TRIANGLE_REPULSIONS_PENALTY<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
}
//#####################################################################
//#####################################################################
// Function Penalty
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS_PENALTY<TV>::
Penalty(T original_volume,const VECTOR<int,4>& nodes,const ARRAY_VIEW<TV>&X,T& e,VECTOR<TV,4>& de,VECTOR<VECTOR<MATRIX<T,TV::m>,4>,4>& he)
{
    auto u=From_Var<3,0>(X(nodes(1))-X(nodes(0)));
    auto v=From_Var<3,1>(X(nodes(2))-X(nodes(0)));
    auto w=From_Var<3,2>(X(nodes(3))-X(nodes(0)));
    auto a=u.Dot(v.Cross(w));
    auto d=1-a/From_Const<TV,3>(original_volume);
    auto ee=stiffness*sqr(d)*d;
    e=ee.x;
    for(int i=1;i<4;i++){
        TV t=ee.dx(i-1);
        de(i)=t;
        de(0)-=t;}
    for(int i=1;i<4;i++){
        for(int j=1;j<4;j++){
            MATRIX<T,3> t=ee.ddx(i-1,j-1);
            he(i)(j)=t;
            he(i)(0)-=t;
            he(0)(0)+=t;}
        he(0)(i)=he(i)(0).Transposed();}
}
namespace PhysBAM{
template class TRIANGLE_REPULSIONS_PENALTY<VECTOR<float,3> >;
template class TRIANGLE_REPULSIONS_PENALTY<VECTOR<double,3> >;
}
