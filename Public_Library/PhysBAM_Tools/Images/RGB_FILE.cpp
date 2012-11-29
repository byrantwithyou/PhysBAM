//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Images/RGB_FILE.h>
#include <PhysBAM_Tools/Images/RGB_HEADER.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <cstring>
using namespace PhysBAM;

//#####################################################################
// Function Read
//#####################################################################
template<class T> void RGB_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    std::istream* input(FILE_UTILITIES::Safe_Open_Input(filename,true));
    RGB_HEADER header;Read_Binary<T>(*input,header);
    image.Resize(0,header.width,0,header.height);unsigned char byte;
    VECTOR<int,2> counts=image.domain.Edge_Lengths();
    if(!header.compression){
        VECTOR<unsigned char,3> color_byte;
        for(int j=0;j<counts.y;j++) for(int i=0;i<counts.x;i++){Read_Binary<T>(*input,byte);image(i,j).x=IMAGE<T>::Byte_Color_To_Scalar_Color(byte);}
        for(int j=0;j<counts.y;j++) for(int i=0;i<counts.x;i++){Read_Binary<T>(*input,byte);image(i,j).y=IMAGE<T>::Byte_Color_To_Scalar_Color(byte);}
        for(int j=0;j<counts.y;j++) for(int i=0;i<counts.x;i++){Read_Binary<T>(*input,byte);image(i,j).z=IMAGE<T>::Byte_Color_To_Scalar_Color(byte);}}
    else{
        unsigned char pixel,count;ARRAY<unsigned int> offset(header.height*header.channels),length(header.height*header.channels);
        unsigned int max_offset=-1,max_offset_index=-1;
        for(int k=0;k<counts.y*header.channels;k++){Read_Binary<T>(*input,offset(k));Swap_Endianity(offset(k));if(offset(k)>max_offset){max_offset=offset(k);max_offset_index=k;}}
        for(int k=0;k<counts.y*header.channels;k++){Read_Binary<T>(*input,length(k));Swap_Endianity(length(k));}
        // read the rest of the file into memory...
        ARRAY<unsigned char> data(offset(max_offset_index)+counts.y*2);
        int current_byte=512+counts.y*header.channels*2*4;
        while(current_byte<data.m) Read_Binary<T>(*input,data(current_byte++));
        // unpack runs
        for(int band=0;band<3;band++)for(int j=0;j<counts.y;j++){
            int index=offset(j+counts.y*band),end_index=index+length(j+counts.y*band);int column=1;
            for(int i=0;index<end_index;i++){
                pixel=data(index++);count=pixel & 0x7f;if(count==0) break;
                if(pixel & 0x80){for(int k=0;k<count;k++){Read_Binary<T>(*input,byte);image(column,j)[band]=IMAGE<T>::Byte_Color_To_Scalar_Color(data(index++));column++;}}
                else{T float_pixel=IMAGE<T>::Byte_Color_To_Scalar_Color(data(index++));for(int k=0;k<count;k++){image(column,j)[band]=float_pixel;column++;}}}}}
    delete input;
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void RGB_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("RGB_FILE: cannot read alpha");
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void RGB_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{  
    PHYSBAM_ASSERT(image.domain.min_corner.x==0 && image.domain.min_corner.y==0);
    VECTOR<int,2> counts=image.domain.Edge_Lengths();
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename,true,false); // no compression
    RGB_HEADER header;header.Initialize(counts.x,counts.y);header.channels=d;
    Write_Binary<T>(*output,header);
    for(int j=0;j<counts.y;j++) for(int i=0;i<counts.x;i++) Write_Binary<T>(*output,IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j)[0]));
    for(int j=0;j<counts.y;j++) for(int i=0;i<counts.x;i++) Write_Binary<T>(*output,IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j)[1]));
    for(int j=0;j<counts.y;j++) for(int i=0;i<counts.x;i++) Write_Binary<T>(*output,IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j)[2]));
    if(d==4) for(int j=0;j<counts.y;j++) for(int i=0;i<counts.x;i++) Write_Binary<T>(*output,IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j)[3]));
    delete output;
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool RGB_FILE<T>::
Is_Supported()
{
    return true;
}
//#####################################################################
template class RGB_FILE<float>;
template void RGB_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&);
template void RGB_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&);
template class RGB_FILE<double>;
template void RGB_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&);
template void RGB_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&);
