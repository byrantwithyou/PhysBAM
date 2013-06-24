//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Read_Write/GNUPLOT_OUTPUT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
GNUPLOT_OUTPUT::
GNUPLOT_OUTPUT()
{}
//#####################################################################
// Destructor
//#####################################################################
GNUPLOT_OUTPUT::
~GNUPLOT_OUTPUT()
{}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void GNUPLOT_OUTPUT::
Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str());
    for(int i=0;i<m;i++){VECTOR<T,1> position=grid.X(VECTOR<int,1>(i));Matlab_Output<<position.x<<"\t"<<output(i)<<std::endl;}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void GNUPLOT_OUTPUT::
Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x,n=output.domain.max_corner.y;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str());
    for(int i=0;i<m;i++){for(int j=0;j<n;j++){
            VECTOR<T,2> position=grid.X(VECTOR<int,2>(i,j));Matlab_Output<<position.x<<"\t"<<position.y<<"\t"<<output(i,j)<<std::endl;}
        Matlab_Output<<std::endl;}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void GNUPLOT_OUTPUT::
Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x,n=output.domain.max_corner.y,mn=output.domain.max_corner.z;
    std::ofstream Matlab_Output;Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str());
    for(int i=0;i<m;i++){for(int j=0;j<n;j++){for(int ij=0;ij<mn;ij++){
            VECTOR<T,3> position=grid.X(VECTOR<int,3>(i,j,ij));Matlab_Output<<position.x<<"\t"<<position.y<<"\t"<<position.z<<"\t"<<output(i,j,ij)<<std::endl;}
            Matlab_Output<<std::endl;}
        Matlab_Output<<std::endl;}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T,int d> void GNUPLOT_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,d> >& X,const int stepnumber)
{
    std::ofstream Matlab_Output;
    Matlab_Output.open(STRING_UTILITIES::string_sprintf("%s.%d",file_name.c_str(),stepnumber).c_str());
    for(int i=0;i<X.Size();i++){for(int a=0;a<d;a++){
        Matlab_Output<<X(i)[a];if(a!=d) Matlab_Output<<"\t";}
        Matlab_Output<<std::endl;}
    Matlab_Output.close();
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,1> >& X,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,2> >& X,const int stepnumber); \
    template void GNUPLOT_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,3> >& X,const int stepnumber);

INSTANTIATION_HELPER(float);
INSTANTIATION_HELPER(double);

