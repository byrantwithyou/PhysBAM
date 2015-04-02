//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION
//##################################################################### 
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION<TV>::
FLUID_TO_SOLID_INTERPOLATION(const COLLISION_AWARE_INDEX_MAP<TV>& map,const DEFORMABLE_PARTICLES<TV>& particles_input)
    :FLUID_TO_SOLID_INTERPOLATION_BASE<TV>(map),max_dist(2),particles(particles_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_TO_SOLID_INTERPOLATION<TV>::
~FLUID_TO_SOLID_INTERPOLATION()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Compute_Weights(const TV& X,int axis,ARRAY<ENTRY>& array)
{
    T limit=max_dist*index_map.grid.dX.Max(),total_weight=0;
    for(FACE_ITERATOR<TV> it(index_map.grid,RANGE<TV_INT>(index_map.grid.Index(X)).Thickened(3),axis);it.Valid();it.Next()){
        T dist=(it.Location()-X).Magnitude();
        if(dist>limit) continue;
        int index=index_map.face_indices(it.Full_Index());
        if(index<0) continue;
        ENTRY e={LEVELSET_UTILITIES<T>::Delta(dist,limit),index};
        total_weight+=e.w;
        array.Append(e);}
    PHYSBAM_ASSERT(total_weight>0);
    for(int i=0;i<array.m;i++) array(i).w/=total_weight;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Compute(const int ghost_cells)
{
    entries.Resize(coupled_particles.m);
    for(int i=0;i<entries.m;i++)
        for(int a=0;a<TV::m;a++)
            Compute_Weights(particles.X(coupled_particles(i)),a,entries(i)(a));
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Times_Add(const ARRAY<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity) const
{
    for(int j=0;j<coupled_particles.m;j++){int p=coupled_particles(j);
        for(int a=0;a<TV::m;a++){
            T& v=solid_velocity.V.array(p)(a);
            const ARRAY<ENTRY>& array=entries(p)(a);
            for(int i=0;i<array.m;i++)
                v+=array(i).w*fluid_velocity(array(i).i);}}
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Transpose_Times_Add(const GENERALIZED_VELOCITY<TV>& solid_force,ARRAY<T>& fluid_force) const
{
    for(int j=0;j<coupled_particles.m;j++){int p=coupled_particles(j);
        for(int a=0;a<TV::m;a++){
            T v=solid_force.V.array(p)(a);
            const ARRAY<ENTRY>& array=entries(p)(a);
            for(int i=0;i<array.m;i++)
                fluid_force(array(i).i)+=array(i).w*v;}}
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void FLUID_TO_SOLID_INTERPOLATION<TV>::
Print_Each_Matrix(int n,int fluid_faces,GENERALIZED_VELOCITY<TV>& G) const
{
    OCTAVE_OUTPUT<T> oo(LOG::sprintf("H-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("H",G.Raw_Size(),fluid_faces);
    ARRAY<int> reverse_map_deformable(G.V.array.Size());
    reverse_map_deformable.Subset(G.V.indices)=IDENTITY_ARRAY<>(G.V.Size());

    for(int j=0;j<coupled_particles.m;j++){int p=coupled_particles(j);
        for(int a=0;a<TV::m;a++){
            const ARRAY<ENTRY>& array=entries(p)(a);
            for(int i=0;i<array.m;i++)
                oo.Add_Sparse_Entry((reverse_map_deformable(p)-1)*TV::m+a,array(i).i,array(i).w);}}

    oo.End_Sparse_Matrix();
}
//#####################################################################
namespace PhysBAM{
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<float,1> >;
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<float,2> >;
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<float,3> >;
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<double,1> >;
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<double,2> >;
template class FLUID_TO_SOLID_INTERPOLATION<VECTOR<double,3> >;
}
