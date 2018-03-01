//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Grid_Tools/Debug_Output/MATLAB_OUTPUT.h>
#include <Grid_Tools/Grids/GRID.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
MATLAB_OUTPUT::
MATLAB_OUTPUT()
{
    Check_Endian();
}
//#####################################################################
// Destructor
//#####################################################################
MATLAB_OUTPUT::
~MATLAB_OUTPUT()
{}
//#####################################################################
// Function Check_Endian
//#####################################################################
void MATLAB_OUTPUT::
Check_Endian()
{
    unsigned short one=0x0001;
    if(*((char*)&one)) little_endian=1; // little endian (e.g. PC), 01 comes before 00
    else little_endian=0;                    // big endian (e.g. SGI), 00 comes before 01
}
//#####################################################################
// Function Convert_Bytes
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Convert_Bytes(T& data) const
{
    union {T data;char raw[sizeof(T)];} input,output;
    input.data=data;
    for(int k=0;k<(int)sizeof(T);k++) output.raw[k]=input.raw[sizeof(T)-k-1];
    data=output.data;
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const int stepnumber)
{
    int m=grid.counts.x;
    ARRAY<T,VECTOR<int,1> > x(grid.counts);
    for(int i=0;i<m;i++) x(i)=grid.X(VECTOR<int,1>(i)).x;
    Write_Header_File(file_name,x,stepnumber);
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& x,const int stepnumber)
{
    int m=x.domain.max_corner.x,data_int;double data_double;
    std::ofstream Matlab_Output;
    Matlab_Output.open(LOG::sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    data_int=m;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    for(int i=0;i<m;i++){data_double=x(i);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const int stepnumber)
{
    int m=grid.counts.x,n=grid.counts.y;
    ARRAY<T,VECTOR<int,2> > x(grid.counts),y(grid.counts);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){x(i,j)=grid.X(VECTOR<int,2>(i,j)).x;y(i,j)=grid.X(VECTOR<int,2>(i,j)).y;}
    Write_Header_File(file_name,x,y,stepnumber);
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& x,const ARRAY<T,VECTOR<int,2> >& y,const int stepnumber)
{
    int m=x.domain.max_corner.x,n=x.domain.max_corner.y,data_int;double data_double;
    std::ofstream Matlab_Output;
    Matlab_Output.open(LOG::sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    data_int=m;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    data_int=n;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    int i;for(i=0;i<m;i++) for(int j=0;j<n;j++){data_double=x(i,j);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    for(i=0;i<m;i++) for(int j=0;j<n;j++){data_double=y(i,j);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const int stepnumber)
{
    ARRAY<T,VECTOR<int,3> > x(grid.counts),y(grid.counts),z(grid.counts);
    for(RANGE_ITERATOR<3> it(grid.Domain_Indices());it.Valid();it.Next()){x(it.index)=grid.X(it.index).x;y(it.index)=grid.X(it.index).y;z(it.index)=grid.X(it.index).z;}
    for(RANGE_ITERATOR<3> it(grid.Domain_Indices());it.Valid();it.Next()){x(it.index)=grid.X(it.index).x;y(it.index)=grid.X(it.index).y;z(it.index)=grid.X(it.index).z;}
    Write_Header_File(file_name,x,y,z,stepnumber);
}
//#####################################################################
// Function Write_Header_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& x,const ARRAY<T,VECTOR<int,3> >& y,const ARRAY<T,VECTOR<int,3> >& z,const int stepnumber)
{
    int m=x.domain.max_corner.x,n=x.domain.max_corner.y,mn=x.domain.max_corner.z,data_int;double data_double;
    std::ofstream Matlab_Output;Matlab_Output.open(LOG::sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    data_int=m;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    data_int=n;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    data_int=mn;if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    int i;for(i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++){
        data_double=x(i,j,ij);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    for(i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++){
        data_double=y(i,j,ij);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    for(i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++){
        data_double=z(i,j,ij);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x;double data_double;
    std::ofstream Matlab_Output;Matlab_Output.open(LOG::sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    for(int i=0;i<m;i++){data_double=output(i);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber)
{
    int m=output.domain.max_corner.x,n=output.domain.max_corner.y;double data_double;
    std::ofstream Matlab_Output;Matlab_Output.open(LOG::sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){data_double=output(i,j);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T> void MATLAB_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber)
{
    double data_double;
    std::ofstream Matlab_Output;Matlab_Output.open(LOG::sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    for(RANGE_ITERATOR<3> it(output.domain);it.Valid();it.Next()){
        data_double=output(it.index);if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
// Function Write_Output_File
//#####################################################################
template<class T,int d> void MATLAB_OUTPUT::
Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,d> >& X,const int stepnumber)
{
    std::ofstream Matlab_Output;
    Matlab_Output.open(LOG::sprintf("%s.%d",file_name.c_str(),stepnumber).c_str(),std::ios::out|std::ios::binary);
    int data_int=X.Size();if(!little_endian) Convert_Bytes(data_int);Matlab_Output.write((const char*)&data_int,4);
    for(int a=0;a<d;a++) for(int k=0;k<X.Size();k++){
        double data_double=X(k)[a];if(!little_endian) Convert_Bytes(data_double);Matlab_Output.write((const char*)&data_double,8);}
    Matlab_Output.close();
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template void MATLAB_OUTPUT::Convert_Bytes(T& data) const; \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& x,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& x,const ARRAY<T,VECTOR<int,2> >& y,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Header_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& x,const ARRAY<T,VECTOR<int,3> >& y,const ARRAY<T,VECTOR<int,3> >& z,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,1> >& X,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,2> >& X,const int stepnumber); \
    template void MATLAB_OUTPUT::Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,3> >& X,const int stepnumber);

INSTANTIATION_HELPER(float);
INSTANTIATION_HELPER(double);
