//#####################################################################
// Copyright 2002-2006, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Tools/Images/BMP_FILE.h>
#include <Tools/Images/BMP_HEADER.h>
#include <Tools/Images/IMAGE.h>
using namespace PhysBAM;
//#####################################################################
// Function Read
//#####################################################################
template<class T> void BMP_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("BMP_FILE: Cannot read alpha channel");
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void BMP_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    std::istream* input=Safe_Open_Input_Raw(filename);
    BMP_HEADER header;Read_Binary<T>(*input,header);
    image.Resize(VECTOR<int,2>(header.w,header.h),no_init);
    int line_width=header.w*3,line_padding=((line_width+3)&~3)-line_width;
    input->seekg(header.offset,std::ios::beg);
    VECTOR<unsigned char,3> color_byte;
    for(int j=0;j<header.h;j++){
        for(int i=0;i<header.w;i++){Read_Binary<T>(*input,color_byte);image(i,j)=IMAGE<T>::Byte_Color_To_Scalar_Color(VECTOR<unsigned char,3>(color_byte.z,color_byte.y,color_byte.x));}
        input->seekg(line_padding,std::ios::cur);}
    delete input;
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void BMP_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{  
    PHYSBAM_ASSERT(image.domain.min_corner.x==0 && image.domain.min_corner.y==0);
    std::ostream* output=Safe_Open_Output_Raw(filename,true,false); // no compression
    VECTOR<int,2> counts=image.domain.Edge_Lengths();
    BMP_HEADER header;
    header.Initialize(counts.x,counts.y);
    Write_Binary<T>(*output,header);
    int line_width=header.w*3,line_padding=((line_width+3)&~3)-line_width;
    for(int j=0;j<header.h;j++){
        for(int i=0;i<header.w;i++){VECTOR<unsigned char,d> byte=IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j));Write_Binary<T>(*output,VECTOR<unsigned char,3>(byte[2],byte[1],byte[0]));}
        for(int i=0;i<line_padding;i++) Write_Binary<T>(*output,(unsigned char)0);}
    delete output;
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool BMP_FILE<T>::
Is_Supported()
{
    return true;
}
//#####################################################################
namespace PhysBAM{
template class BMP_FILE<float>;
template void BMP_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&);
template void BMP_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&);
template class BMP_FILE<double>;
template void BMP_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&);
template void BMP_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&);
}
