//#####################################################################
// Copyright 2010, Craig Schroeder, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED
//##################################################################### 
#include <Core/Arrays/CONSTANT_ARRAY.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <Dynamics/Coupled_Evolution/GENERALIZED_FLUID_MASS.h>
#include <Dynamics/Coupled_Evolution/MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input)
    :MATRIX_FLUID_INTERPOLATION_BASE<TV>(index_map_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
~MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED()
{
}
//#####################################################################
// Function Number_Of_Constraints
//#####################################################################
template<class TV> COUPLING_CONSTRAINT_ID MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
Number_Of_Constraints() const
{
    return COUPLING_CONSTRAINT_ID(stencils.m*d);
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
Compute(int ghost_cells)
{
    stencils.Resize(entries.m);
    for(int i=0;i<entries.m;i++){
        int sign=(entries(i).side==0)?-1:1;
        FACE_INDEX<d> face=entries(i).face_index;
        TV_INT cell=face.Cell_Index(entries(i).side);
        T face_X=index_map.grid.Face(face)(face.axis),X=entries(i).X(face.axis),cell_X=index_map.grid.X(cell)(face.axis);
        T alpha=(X-face_X)*sign/index_map.grid.dX(face.axis);
        PHYSBAM_ASSERT(index_map.face_indices(face));
        stencils(i).s(face.axis)(0)=PAIR<int,T>(index_map.face_indices(face),1-alpha);
        face.index(face.axis)+=sign;
        int face_index=index_map.face_indices(face);
        PHYSBAM_ASSERT(face_index);
        stencils(i).s(face.axis)(0)=PAIR<int,T>(face_index,alpha);
        stencils(i).s(face.axis)(1)=PAIR<int,T>(face_index,0);
        stencils(i).s(face.axis)(2)=PAIR<int,T>(face_index,0);        

        alpha=(X-cell_X)*sign/index_map.grid.dX(face.axis);
        T ia[2]={(T)(0.5*(1-alpha)),(T)(0.5*alpha)};
        for(int a=0;a<d;a++) if(a!=face.axis) for(int s=0;s<2;s++) for(int j=0;j<2;j++){
            FACE_INDEX<d> current_face(a,cell);
            current_face.index(face.axis)+=sign*j;
            current_face.index(a)+=s;
            PHYSBAM_ASSERT(index_map.face_indices(current_face));
            stencils(i).s(a)(s*2+j)=PAIR<int,T>(index_map.face_indices(current_face),ia[j]);}}
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
Times_Add(const ARRAY<T>& faces,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const
{
    for(int i=0;i<stencils.m;i++) for(int a=0;a<d;a++){
        T& constraint=constraints(COUPLING_CONSTRAINT_ID(d*i+a));
        const VECTOR<PAIR<int,T>,4>& e=stencils(i).s(a);
        for(int j=0;j<4;j++) constraint+=faces(e(j).x)*e(j).y;}
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
Transpose_Times_Add(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,ARRAY<T>& faces) const
{
    for(int i=0;i<stencils.m;i++) for(int a=0;a<d;a++){
        T constraint=constraints(COUPLING_CONSTRAINT_ID(d*i+a));
        const VECTOR<PAIR<int,T>,4>& e=stencils(i).s(a);
        for(int j=0;j<4;j++) faces(e(j).x)+=constraint*e(j).y;}
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
Print() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
Print_Each_Matrix(int n) const
{
    OCTAVE_OUTPUT<T> oo(LOG::sprintf("W-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("W",Value(Number_Of_Constraints()),index_map.Number_Faces());

    for(int i=0;i<stencils.m;i++) for(int a=0;a<d;a++){
        const VECTOR<PAIR<int,T>,4>& e=stencils(i).s(a);
        for(int j=0;j<4;j++) oo.Add_Sparse_Entry(d*i+a,e(j).x,e(j).y);}

    oo.End_Sparse_Matrix();
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    for(int i=0;i<stencils.m;i++) for(int a=0;a<d;a++){
        const VECTOR<PAIR<int,T>,4>& e=stencils(i).s(a);
        for(int j=0;j<4;j++) data.Append(TRIPLE<int,int,T>(d*i+a,e(j).x,e(j).y));}
}
//#####################################################################
// Function Add_Diagonal
//#####################################################################
template<class TV> void MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::
Add_Diagonal(ARRAY<T,COUPLING_CONSTRAINT_ID>& diagonal,const GENERALIZED_FLUID_MASS<TV>& fluid_mass) const
{
    for(int i=0;i<stencils.m;i++) for(int a=0;a<d;a++){
        T& diag=diagonal(COUPLING_CONSTRAINT_ID(d*i+a));
        const VECTOR<PAIR<int,T>,4>& e=stencils(i).s(a);
        for(int j=0;j<4;j++) diag+=fluid_mass.one_over_fluid_mass_at_faces(e(j).x)*sqr(e(j).y);}
}
//#####################################################################
namespace PhysBAM{
template class MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<VECTOR<float,1> >;
template class MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<VECTOR<float,2> >;
template class MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<VECTOR<float,3> >;
template class MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<VECTOR<double,1> >;
template class MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<VECTOR<double,2> >;
template class MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<VECTOR<double,3> >;
}
