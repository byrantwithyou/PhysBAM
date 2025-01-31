//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Incompressible/Incompressible_Flows/LEVELSET_INDEX_MAP_UNIFORM.h>
#include <Incompressible/Solids_And_Fluids/BOUNDARY_CONDITIONS_CALLBACKS.h>
using namespace PhysBAM;

template<class TV> LEVELSET_INDEX_MAP_UNIFORM<TV>::
//#####################################################################
// Constructor
//#####################################################################
LEVELSET_INDEX_MAP_UNIFORM(const GRID<TV>& grid_input,BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input)
    :grid(grid_input),callback(callback_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_INDEX_MAP_UNIFORM<TV>::
~LEVELSET_INDEX_MAP_UNIFORM()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void LEVELSET_INDEX_MAP_UNIFORM<TV>::
Compute(int axis,VECTOR<bool,d> periodic_boundary_input)
{
    inside.Resize(grid,1);
    face_to_index.Resize(grid);
    face_to_index.Fill(0);
    index_to_face.Remove_All();
    inside.Fill(false);
    callback->Mark_Outside(inside);
    periodic_boundary=periodic_boundary_input;
    for(FACE_ITERATOR<TV> it(grid,0,GRID<TV>::WHOLE_REGION,-1,axis);it.Valid();it.Next()){FACE_INDEX<d> index=it.Full_Index();
        bool& b=inside(index);
        b=!b;
        if(b && !(periodic_boundary[index.axis] && index.index[index.axis]==0))
            face_to_index(index)=index_to_face.Append(index);}

    for(int a=0;a<d;a++) if(periodic_boundary[a])
        for(FACE_ITERATOR<TV> it(grid,0,GRID<TV>::BOUNDARY_REGION,2*a,a);it.Valid();it.Next()){
            FACE_INDEX<d> index=it.Full_Index();
            FACE_INDEX<d> matching_index(index);matching_index.index(a)=1;
            assert(face_to_index(index));
            face_to_index(matching_index)=face_to_index(index);}
}
//#####################################################################
// Function Gather
//#####################################################################
template<class TV> void LEVELSET_INDEX_MAP_UNIFORM<TV>::
Gather(const ARRAY<T,FACE_INDEX<d> >& u,ARRAY<T>& v) const
{
    for(int i=0;i<index_to_face.m;i++) v(i)=u(index_to_face(i));
}
//#####################################################################
// Function Scatter
//#####################################################################
template<class TV> void LEVELSET_INDEX_MAP_UNIFORM<TV>::
Scatter(const ARRAY<T>& u,ARRAY<T,FACE_INDEX<d> >& v) const
{
    for(int i=0;i<index_to_face.m;i++) v(index_to_face(i))=u(i);
    for(int a=0;a<d;a++) if(periodic_boundary[a])
        for(FACE_ITERATOR<TV> it(grid,0,GRID<TV>::BOUNDARY_REGION,2*a,a);it.Valid();it.Next()){
            FACE_INDEX<d> index=it.Full_Index();
            FACE_INDEX<d> matching_index(index);matching_index.index(a)=1;
            v(matching_index)=v(index);}
}
namespace PhysBAM{
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<float,1> >;
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<float,2> >;
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<float,3> >;
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<double,1> >;
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<double,2> >;
template class LEVELSET_INDEX_MAP_UNIFORM<VECTOR<double,3> >;
}
