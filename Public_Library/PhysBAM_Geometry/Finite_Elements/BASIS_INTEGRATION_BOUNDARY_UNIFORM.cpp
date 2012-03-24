//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
#include <PhysBAM_Tools/Symbolics/MULTIVARIATE_POLYNOMIAL.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BASIS_INTEGRATION_BOUNDARY_UNIFORM<TV>::
BASIS_INTEGRATION_BOUNDARY_UNIFORM(const GRID<TV>& grid_input,const BASIS_STENCIL_UNIFORM<TV>& stencil_input,const BASIS_STENCIL_BOUNDARY_UNIFORM<TV>& boundary_input,CELL_MAPPING<TV>& cm_input)
    :grid(grid_input),stencil(stencil_input),boundary(boundary_input),cm(cm_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BASIS_INTEGRATION_BOUNDARY_UNIFORM<TV>::
~BASIS_INTEGRATION_BOUNDARY_UNIFORM()
{
}
//#####################################################################
// Function Compute_Matrix
//#####################################################################
template<class TV> void BASIS_INTEGRATION_BOUNDARY_UNIFORM<TV>::
Compute_Matrix(SYSTEM_MATRIX_HELPER<T>& helper)
{
    helper.New_Block();
    ARRAY<T_FACE> clipped_simplices;
    for(int e=0;e<boundary.object.mesh.elements.m;e++){
        T_FACE face=boundary.object.Get_Element(e);
        RANGE<TV> box=face.Bounding_Box();
        for(int i=0;i<stencil.diced.m;i++){
            const typename BASIS_STENCIL_UNIFORM<TV>::DICED& entry=stencil.diced(i);
            RANGE<TV> T_cell_range((box-grid.domain.min_corner)*grid.one_over_dX - (T).5*RANGE<TV>(entry.range+(TV_INT()+1)));
            RANGE<TV_INT> cell_range(TV_INT(ceil(T_cell_range.min_corner+1e-14)),TV_INT(floor(T_cell_range.max_corner-1e-14))+1);
            for(UNIFORM_ARRAY_ITERATOR<TV::m> it(cell_range);it.Valid();it.Next()){
                TV_INT index=it.index+entry.index_offset;
                RANGE<TV> domain=RANGE<TV>(entry.range+(1+2*it.index))*((T).5*grid.dX)+grid.domain.min_corner;
                clipped_simplices.Remove_All();
                face.Clip_To_Box(domain,clipped_simplices);
                T integral=0;
                for(int i=0;i<clipped_simplices.m;i++){
                    for(int j=0;j<TV::m;j++)
                        clipped_simplices(i).X(j)-=grid.Center(it.index);
                    integral+=entry.polynomial.Integrate_Over_Primitive(reinterpret_cast<const VECTOR<TV,TV::m>&>(clipped_simplices(i).X(0)));}
                if(fabs(integral)<1e-14) continue;
                helper.data.Append(TRIPLE<int,int,T>(cm.Get_Index(index,0),e,integral));
                helper.data.Append(TRIPLE<int,int,T>(cm.Get_Index(index,1),e,-integral));}}}
}
template class BASIS_INTEGRATION_BOUNDARY_UNIFORM<VECTOR<float,2> >;
template class BASIS_INTEGRATION_BOUNDARY_UNIFORM<VECTOR<float,3> >;
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class BASIS_INTEGRATION_BOUNDARY_UNIFORM<VECTOR<double,2> >;
template class BASIS_INTEGRATION_BOUNDARY_UNIFORM<VECTOR<double,3> >;
#endif
