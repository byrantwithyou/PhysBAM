//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_GRADIENT_BASE
//##################################################################### 
#include <Tools/Log/LOG.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT_BASE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_FLUID_GRADIENT_BASE<TV>::
MATRIX_FLUID_GRADIENT_BASE(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input)
    :index_map(index_map_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_FLUID_GRADIENT_BASE<TV>::
~MATRIX_FLUID_GRADIENT_BASE()
{
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT_BASE<TV>::
Times_Add(const ARRAY<T>& cells,ARRAY<T>& faces) const
{
    gradient.Times_Add(cells,faces);
}
//#####################################################################
// Function Times
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT_BASE<TV>::
Times(const ARRAY<T>& cells,ARRAY<T>& faces) const
{
    gradient.Times(cells,faces);
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT_BASE<TV>::
Transpose_Times_Add(const ARRAY<T>& faces,ARRAY<T>& cells) const
{
    gradient.Transpose_Times_Add(faces,cells);
}
//#####################################################################
// Function Times
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT_BASE<TV>::
Transpose_Times(const ARRAY<T>& faces,ARRAY<T>& cells) const
{
    gradient.Transpose_Times(faces,cells);
}
//#####################################################################
// Function Collect_Maxabs_Velocity
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT_BASE<TV>::
Collect_Maxabs_Velocity(const ARRAY<T>& faces,ARRAY<T>& cells) const
{
    cells.Fill(0);
    int index=gradient.offsets(0);
    for(int i=0;i<gradient.m;i++){
        int end=gradient.offsets(i+1);T y=faces(i);
        for(;index<end;index++) cells(gradient.A(index).j)=max(cells(gradient.A(index).j),abs(y*gradient.A(index).a));}
}
//#####################################################################
// Function Test_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT_BASE<TV>::
Test_Matrix() const
{
    RANDOM_NUMBERS<T> random;
    ARRAY<T> a(gradient.m),a2(gradient.m),b(gradient.n),b2(gradient.n);
    random.Fill_Uniform(a,-1,1);
    random.Fill_Uniform(b,-1,1);

    Times(b,a2);
    Transpose_Times(a,b2);

    T inner1=ARRAY<T>::Dot_Product(a,a2);
    T inner2=ARRAY<T>::Dot_Product(b,b2);
    LOG::cout<<"MATRIX_FLUID_GRADIENT Test: "<<inner1<<"  vs  "<<inner2<<"  relative  "<<abs(inner1-inner2)/maxabs((T)1e-30,inner1,inner2)<<std::endl;
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT_BASE<TV>::
Print_Each_Matrix(int n) const
{
    OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("G-%i.txt",n).c_str()).Write("G",gradient);
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT_BASE<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    for(int i=0;i<gradient.m;i++){
        int s=gradient.offsets(i),e=gradient.offsets(i+1);
        for(int j=s;j<e;j++)
            data.Append(TRIPLE<int,int,T>(i,gradient.A(j).j,gradient.A(j).a));}
}
namespace PhysBAM{
template class MATRIX_FLUID_GRADIENT_BASE<VECTOR<float,1> >;
template class MATRIX_FLUID_GRADIENT_BASE<VECTOR<float,2> >;
template class MATRIX_FLUID_GRADIENT_BASE<VECTOR<float,3> >;
template class MATRIX_FLUID_GRADIENT_BASE<VECTOR<double,1> >;
template class MATRIX_FLUID_GRADIENT_BASE<VECTOR<double,2> >;
template class MATRIX_FLUID_GRADIENT_BASE<VECTOR<double,3> >;
}
