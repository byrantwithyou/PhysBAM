//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PNG_FILE  
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Images/PNG_FILE.h>
#ifdef USE_LIBPNG
#include <png.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Read/Write stubs for case of no libpng
//#####################################################################
#ifndef USE_LIBPNG

template<class T> void PNG_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("Not compiled with USE_LIBPNG.  Cannot read png image.");
}
template<class T> void PNG_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("Not compiled with USE_LIBPNG.  Cannot read png image.");
}

template<class T> template<int d> void PNG_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("Not compiled with USE_LIBPNG.  Cannot write png image.");
}

template<class T> bool PNG_FILE<T>::
Is_Supported()
{return false;}

#else
//#####################################################################
// Function Read
//#####################################################################
template<class T> void PNG_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    FILE* file=fopen(filename.c_str(),"rb");
    if(!file) throw READ_ERROR(LOG::sprintf("Failed to open %s for reading",filename.c_str()));

    png_structp png_ptr=png_create_read_struct(PNG_LIBPNG_VER_STRING,0,0,0);
    if(!png_ptr){
        fclose(file);
        png_destroy_read_struct(&png_ptr,0,0);
        throw READ_ERROR(LOG::sprintf("Error reading png file %s",filename.c_str()));}
    png_infop info_ptr=png_create_info_struct(png_ptr);
    if(!info_ptr){
        fclose(file);
        png_destroy_read_struct(&png_ptr,&info_ptr,0);
        throw READ_ERROR(LOG::sprintf("Error reading png file %s",filename.c_str()));}
    if(setjmp(png_jmpbuf(png_ptr))){
        fclose(file);
        png_destroy_read_struct(&png_ptr,&info_ptr,0);
        throw READ_ERROR(LOG::sprintf("Error reading png file %s",filename.c_str()));}
    png_init_io(png_ptr,file);
    png_read_png(png_ptr,info_ptr,PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_STRIP_ALPHA | PNG_TRANSFORM_PACKING,0);
    int width=png_get_image_width(png_ptr,info_ptr),height=png_get_image_height(png_ptr,info_ptr);
    int color_type=png_get_color_type(png_ptr,info_ptr);
    if(color_type!=PNG_COLOR_TYPE_RGB && color_type!=PNG_COLOR_TYPE_RGBA) PHYSBAM_FATAL_ERROR("PNG read only supports RGB and RGBA");
    image.Resize(VECTOR<int,2>(width,height));
    VECTOR<unsigned char,3>** row_pointers=(VECTOR<unsigned char,3>**)png_get_rows(png_ptr,info_ptr);
    for(int i=0;i<width;i++)for(int j=0;j<height;j++)image(i,j)=IMAGE<T>::Byte_Color_To_Scalar_Color(row_pointers[height-j-1][i]);
        
    png_destroy_read_struct(&png_ptr,&info_ptr,0);
    fclose(file);
    return;
}
template<class T> void PNG_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void PNG_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{  
    FILE* file=fopen(filename.c_str(),"wb");
    if(!file) PHYSBAM_FATAL_ERROR(LOG::sprintf("Failed to open %s for writing",filename.c_str()));

    png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
    if(!png_ptr) PHYSBAM_FATAL_ERROR(LOG::sprintf("Error writing png file %s",filename.c_str()));
    png_infop info_ptr=png_create_info_struct(png_ptr);
    if(!info_ptr) PHYSBAM_FATAL_ERROR(LOG::sprintf("Error writing png file %s",filename.c_str()));
    if(setjmp(png_jmpbuf(png_ptr))) PHYSBAM_FATAL_ERROR(LOG::sprintf("Error writing png file %s",filename.c_str()));
    png_init_io(png_ptr,file);
    int color_type;
    switch(d){case 3:color_type=PNG_COLOR_TYPE_RGB;break;case 4:color_type=PNG_COLOR_TYPE_RGBA;break;
        default:PHYSBAM_FATAL_ERROR("Invalid number of channels for png write");break;}
    VECTOR<int,2> counts=image.domain.Edge_Lengths();
    png_set_IHDR(png_ptr,info_ptr,counts.x,counts.y,8,color_type,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);

    VECTOR<unsigned char,d>* byte_data=new VECTOR<unsigned char,d>[counts.y*counts.x];
    VECTOR<unsigned char,d>** row_pointers=new VECTOR<unsigned char,d>*[counts.y];
    for(int j=0;j<counts.y;j++){
        row_pointers[counts.y-j-1]=byte_data+counts.x*(counts.y-j-1);
        for(int i=0;i<counts.x;i++) row_pointers[counts.y-j-1][i]=IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j));}
    png_set_rows(png_ptr,info_ptr,(png_byte**)row_pointers);
    png_write_png(png_ptr,info_ptr,PNG_TRANSFORM_IDENTITY,0);
    delete[] row_pointers;delete[] byte_data;

    png_destroy_write_struct(&png_ptr,&info_ptr);fclose(file);return;
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool PNG_FILE<T>::
Is_Supported()
{
    return true;
}
//#####################################################################
#endif
namespace PhysBAM{
template class PNG_FILE<float>;
template void PNG_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&);
template void PNG_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&);
template class PNG_FILE<double>;
template void PNG_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&);
template void PNG_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&);
}
