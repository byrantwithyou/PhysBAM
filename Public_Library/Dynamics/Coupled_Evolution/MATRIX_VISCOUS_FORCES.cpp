//#####################################################################
// Copyright 2010, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_VISCOUS_FORCES
//##################################################################### 
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Log/LOG.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <Dynamics/Coupled_Evolution/MATRIX_VISCOUS_FORCES.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_VISCOUS_FORCES<TV>::
MATRIX_VISCOUS_FORCES(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input)
    :grid(index_map_input.grid),index_map(index_map_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_VISCOUS_FORCES<TV>::
~MATRIX_VISCOUS_FORCES()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Compute(const T dt,const ARRAY<bool,FACE_INDEX<d> >& psi_N,T mu)
{
    // Setting up the correct gradient
    // TODO: Dirchlet cells & Neumann faces
    entries.Remove_All();

    T coefficient=sqrt(dt*mu*grid.One_Over_Cell_Size());
    TV face_areas=grid.Face_Sizes();
    last_id=VISCOUS_FORCE_ID();
    if(!mu) return;

    for(FACE_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index();
        FACE_INDEX<d> this_face(iterator.Full_Index());
        for(int other_axis=0;other_axis<TV::dimension;other_axis++){
            last_id++;
            FACE_INDEX<d> other_face(axis,face_index+TV_INT::Axis_Vector(other_axis));
            T weight=coefficient*face_areas(axis);
            int this_index=-1,other_index=-1;
            if(index_map.face_indices.Valid_Index(this_face)) this_index=index_map.face_indices(this_face);
            else if(index_map.constraint_indices.Contains(SIDED_FACE_INDEX<TV::dimension>(1,this_face))) this_index=index_map.indexed_faces.m+index_map.constraint_indices.Get(SIDED_FACE_INDEX<TV::dimension>(1,this_face));
            if(index_map.face_indices.Valid_Index(other_face)) other_index=index_map.face_indices(other_face);
            else if(index_map.constraint_indices.Contains(SIDED_FACE_INDEX<TV::dimension>(0,other_face))) other_index=index_map.indexed_faces.m+index_map.constraint_indices.Get(SIDED_FACE_INDEX<TV::dimension>(0,other_face));
            if((this_index>index_map.indexed_faces.m && other_index>index_map.indexed_faces.m) || !this_index || !other_index){last_id--;continue;}
            if(this_index>=0) entries.Append(ENTRY(-weight,this_index,last_id));
            if(other_index>=0) entries.Append(ENTRY(weight,other_index,last_id));}}
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Times_Add(const ARRAY<T>& velocities,ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients) const
{
    for(int i=0;i<entries.m;i++){
        const ENTRY& entry=entries(i);
        viscous_force_coefficients(entry.viscous_id)+=entry.weight*velocities(entry.face_index);}
}
//#####################################################################
// Function Times
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Times(const ARRAY<T>& velocities,ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients) const
{
    viscous_force_coefficients.Fill(T());
    Times_Add(velocities,viscous_force_coefficients);
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Transpose_Times_Add(const ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients,ARRAY<T>& velocities) const
{
    for(int i=0;i<entries.m;i++){
        const ENTRY& entry=entries(i);
        velocities(entry.face_index)+=entry.weight*viscous_force_coefficients(entry.viscous_id);}
}
//#####################################################################
// Function Transpose_Times
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Transpose_Times(const ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients,ARRAY<T>& velocities) const
{
    velocities.Fill(T());
    Transpose_Times_Add(viscous_force_coefficients,velocities);
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
// can depend on position too
template<class TV> VISCOUS_FORCE_ID MATRIX_VISCOUS_FORCES<TV>::
Viscous_Forces_Size() const
{
    return last_id;
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Test_Matrix() const
{
    RANDOM_NUMBERS<T> random;
    ARRAY<T,VISCOUS_FORCE_ID> K(last_id),K1(K);
    ARRAY<T> V(index_map.Number_Faces()),V1(V);

    random.Fill_Uniform(K,-1,1);
    random.Fill_Uniform(V,-1,1);
    Times(V,K1);
    Transpose_Times(K,V1);

    T inner_V=V.Dot(V1);
    T inner_K=K.Dot(K1);
    LOG::cout<<"MATRIX_VISCOUS_FORCES Symmetry Test: "<<inner_V<<"  vs  "<<inner_K<<"  relative  "<<abs(inner_V-inner_K)/maxabs((T)1e-30,inner_V,inner_K)<<std::endl;
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Print_Each_Matrix(int n) const
{
    OCTAVE_OUTPUT<T> oo(STRING_UTILITIES::string_sprintf("N-%i.txt",n).c_str());
    oo.Begin_Sparse_Matrix("N",Value(last_id),index_map.Number_Faces());
    for(int i=0;i<entries.m;i++){
        const ENTRY& entry=entries(i);
        oo.Add_Sparse_Entry(Value(entry.viscous_id),entry.face_index,entry.weight);}
    oo.End_Sparse_Matrix();
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void MATRIX_VISCOUS_FORCES<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    for(int i=0;i<entries.m;i++){
        const ENTRY& entry=entries(i);
        data.Append(TRIPLE<int,int,T>(Value(entry.viscous_id),entry.face_index,entry.weight));}
}
//#####################################################################
namespace PhysBAM{
template class MATRIX_VISCOUS_FORCES<VECTOR<float,1> >;
template class MATRIX_VISCOUS_FORCES<VECTOR<float,2> >;
template class MATRIX_VISCOUS_FORCES<VECTOR<float,3> >;
template class MATRIX_VISCOUS_FORCES<VECTOR<double,1> >;
template class MATRIX_VISCOUS_FORCES<VECTOR<double,2> >;
template class MATRIX_VISCOUS_FORCES<VECTOR<double,3> >;
}
