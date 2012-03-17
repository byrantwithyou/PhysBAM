//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Struct BMP_HEADER
//#####################################################################
// A group of fields with the layout of a Windows bitmap file header. Used by the class BMP_FILE.
//#####################################################################
#ifndef __BMP_HEADER__
#define __BMP_HEADER__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
namespace PhysBAM{

struct BMP_HEADER{
    typedef int HAS_UNTYPED_READ_WRITE;
    char file_type[2]; // must be "BM"
    int file_size;
    short reserved1;
    short reserved2;
    int offset;//#bytes from the beginning of the file to the bitmap
    int info_header_size;
    int w;
    int h;
    short number_of_bitplanes;//must be 1
    short bits_per_pixel;//e.g. 1,4,8, or 24
    int type_of_compression; //0 <=> no compression
    int bitmap_size;//w*h*3 if no compression
    int x_pixels_per_meter;
    int y_pixels_per_meter;
    int number_of_colors;
    int number_of_important_colors;

    void Initialize(const int w_input,const int h_input)
    {w=w_input;h=h_input;bitmap_size=w*h*3; //w*h*3 if no compression
    file_type[0]='B';file_type[1]='M';
    file_size=54+w*h*3;
    reserved1=0;reserved2=0;
    offset=54;//#bytes from the beginning of the file to the bitmap
    info_header_size=40;//probably unimportant
    number_of_bitplanes=1;//must be 1
    bits_per_pixel=24;//e.g. 1,4,8, or 24
    type_of_compression=0; //0 <=> no compression
    x_pixels_per_meter=0;//?
    y_pixels_per_meter=0;//?
    number_of_colors=0;
    number_of_important_colors=0;} 

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,file_type[0],file_type[1],file_size,reserved1,reserved2,offset,info_header_size);
    Read_Binary<RW>(input,w,h,number_of_bitplanes,bits_per_pixel,type_of_compression);
    Read_Binary<RW>(input,bitmap_size,x_pixels_per_meter,y_pixels_per_meter,number_of_colors,number_of_important_colors);

    //check validity
    if(file_type[0]!='B' || file_type[1]!='M') PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Illegal file type: %c%c",file_type[0],file_type[1]));
    if(info_header_size != 40) LOG::cerr<<"Warning: weird info_header_size: "<<info_header_size<<" (expected 40)"<<std::endl;
    if(number_of_bitplanes != 1) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Illegal number of bitplanes: %d (must be 1)",number_of_bitplanes));
    if(bits_per_pixel != 24) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Number of bits per pixel: %d (the only supported number is 24)",bits_per_pixel));
    if(type_of_compression != 0) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Type of compression %d (the only supported type is 0)",type_of_compression));}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,file_type[0],file_type[1],file_size,reserved1,reserved2,offset,info_header_size);
    Write_Binary<RW>(output,w,h,number_of_bitplanes,bits_per_pixel,type_of_compression);
    Write_Binary<RW>(output,bitmap_size,x_pixels_per_meter,y_pixels_per_meter,number_of_colors,number_of_important_colors);}
//#####################################################################
};
template<class T>
inline std::ostream& operator<<(std::ostream& output,const BMP_HEADER& header)
{output<<"BMP_HEADER:"<<std::endl
    <<"file_type: "<<header.file_type[0]<<" "<<header.file_type[1]<<std::endl
    <<"file_size  "<<int(header.file_size)<<std::endl
    <<"reserved1  "<<int(header.reserved1)<<std::endl
    <<"reserved2  "<<int(header.reserved2)<<std::endl
    <<"offset     "<<int(header.offset)<<std::endl
    <<"info_header_size    "<<header.info_header_size<<std::endl
    <<"w "<<header.w<<", h "<<header.h<<std::endl
    <<"number_of_bitplanes "<<int(header.number_of_bitplanes)<<std::endl
    <<"bits_per_pixel      "<<int(header.bits_per_pixel)<<std::endl
    <<"type_of_compression "<<header.type_of_compression<<std::endl
    <<"bitmap_size         "<<header.bitmap_size<<std::endl
    <<"x_pixels_per_meter  "<<header.x_pixels_per_meter<<std::endl
    <<"y_pixels_per_meter  "<<header.y_pixels_per_meter<<std::endl
    <<"number_of_colors    "<<header.number_of_colors<<std::endl
    <<"number_of_important_colors "<<header.number_of_important_colors<<std::endl;
return output;}
}
#endif
