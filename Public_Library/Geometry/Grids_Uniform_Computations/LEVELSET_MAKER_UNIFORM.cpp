//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Math_Tools/RANGE.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
//#####################################################################
// Function Compute_Level_Set
//#####################################################################
template<class TV> void LEVELSET_MAKER_UNIFORM<TV>::
Compute_Level_Set(T_SURFACE& surface,GRID<TV>& grid,int ghost_cells,ARRAY<T,TV_INT>& phi)
{
    PHYSBAM_ASSERT(phi.domain.min_corner.Max()<=-ghost_cells);
    if(!surface.mesh.adjacent_elements) surface.mesh.Initialize_Adjacent_Elements();
    phi.Fill(FLT_MAX);
    T dx=grid.dX.Max();

    ARRAY<TV_INT> seed_indices;
    for(int i=0;i<surface.mesh.elements.m;i++){
        typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX simplex(surface.particles.X.Subset(surface.mesh.elements(i)));
        RANGE<TV_INT>  box(grid.Clamp_To_Cell(simplex.Bounding_Box(),ghost_cells).Thickened(1).Intersect(phi.domain));
        for(RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next()){
            TV X=grid.X(it.index);
            T dist=simplex.Distance_To_Element(X),&p=phi(it.index);
            if(dist<abs(p)+dx*(T)1e-4 && dist<dx){
                if(p==FLT_MAX) seed_indices.Append(it.index);
                bool new_sign=simplex.Signed_Distance(X)<0;
                if(abs(dist-abs(p))<dx*(T)1e-4 && new_sign!=(p<0))
                    new_sign=surface.Inside(X);
                if(abs(dist)<abs(p)) p=dist;
                p=abs(p)*(new_sign?-1:1);}}}

    ARRAY<TV_INT> todo;
    for(int i=0;i<seed_indices.m;i++)
        if(phi(seed_indices(i))<0)
            todo.Append(seed_indices(i));
    while(todo.m){
        ARRAY<TV_INT> next;
        for(int i=0;i<todo.m;i++)
            for(int a=0;a<TV::m;a++)
                for(int s=-1;s<2;s+=2){
                    TV_INT index=todo(i);
                    index(a)+=s;
                    if(phi.domain.Lazy_Inside_Half_Open(index) && phi(index)==FLT_MAX){
                        phi(index)=-FLT_MAX;
                        next.Append(index);}}
        next.Exchange(todo);}

    LEVELSET<TV> levelset(grid,phi);
    FAST_MARCHING_METHOD_UNIFORM<TV> fmm(levelset,ghost_cells);
    fmm.Fast_Marching_Method(phi,0,&seed_indices);
}
namespace PhysBAM{
template class LEVELSET_MAKER_UNIFORM<VECTOR<float,2> >;
template class LEVELSET_MAKER_UNIFORM<VECTOR<float,3> >;
template class LEVELSET_MAKER_UNIFORM<VECTOR<double,2> >;
template class LEVELSET_MAKER_UNIFORM<VECTOR<double,3> >;
}
