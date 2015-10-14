//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMAGE
//#####################################################################
#include <Tools/Images/BMP_FILE.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Images/JPG_FILE.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Images/PPM_FILE.h>
#include <Tools/Images/RGB_FILE.h>
#include <Tools/Log/LOG.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Read
//#####################################################################
template<class T> template<int d> void IMAGE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{
    std::string extension=FILE_UTILITIES::Get_File_Extension(filename);
    if(extension=="bmp") BMP_FILE<T>::Read(filename,image);
    else if(extension=="jpg") JPG_FILE<T>::Read(filename,image);
    else if(extension=="ppm") PPM_FILE<T>::Read(filename,image);
    else if(extension=="rgb") RGB_FILE<T>::Read(filename,image);
    else if(extension=="png") PNG_FILE<T>::Read(filename,image);
    else if(extension=="pbi") FILE_UTILITIES::Read_From_File<float>(filename,image);
    else PHYSBAM_FATAL_ERROR(LOG::sprintf("Unknown image file extension  from filename '%s' extension '%s'",filename.c_str(),extension.c_str()));
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void IMAGE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image,const T gamma,const T dither_amplitude)
{
    RANDOM_NUMBERS<T> random;random.Set_Seed(324032); // want uniform seed so noise perturbation pattern is temporally coherent
    ARRAY<VECTOR<T,d> ,VECTOR<int,2> > *corrected_image=0; // may need to shift image to be (0,0) based, or to gamma correct it
    VECTOR<int,2> counts=image.domain.Edge_Lengths();
    if(gamma!=1 || image.domain.min_corner.x!=0 || image.domain.min_corner.y!=0 || dither_amplitude>0){
        corrected_image=new ARRAY<VECTOR<T,d> ,VECTOR<int,2> >(0,counts.x,0,counts.y,false);
        ARRAY<VECTOR<T,d> ,VECTOR<int,2> >::Shifted_Get(*corrected_image,image,VECTOR<int,2>(image.domain.min_corner.x,image.domain.min_corner.y));
        T one_over_gamma=1/gamma;
        for(int t=0;t<corrected_image->array.Size();t++){
            VECTOR<T,d> color=corrected_image->array(t);
            for(int channel=0;channel<d;channel++) corrected_image->array(t)[channel]=pow(color[channel],one_over_gamma);
            if(dither_amplitude){
                VECTOR<T,d> pixel_values((T)256*corrected_image->array(t));
                VECTOR<int,d> floored_values;
                for(int channel=0;channel<d;channel++) floored_values[channel]=(int)pixel_values[channel];
                VECTOR<T,d> random_stuff=random.Get_Uniform_Vector(VECTOR<T,d>(),VECTOR<T,d>::All_Ones_Vector());
                VECTOR<T,d> normalized_values=pixel_values-VECTOR<T,d>(floored_values);
                for(int k=0;k<d;k++)
                    if(random_stuff(k)>normalized_values(k)) corrected_image->array(t)[k]=(floored_values[k]+(T).5001)/(T)256; // use normal quantized floor
                    else corrected_image->array(t)[k]=(floored_values[k]+(T)1.5001)/(T)256;}}} // jump to next value
    const ARRAY<VECTOR<T,d> ,VECTOR<int,2> > &image_to_write=corrected_image?*corrected_image:image;
    std::string extension=FILE_UTILITIES::Get_File_Extension(filename);
    if(extension=="bmp") BMP_FILE<T>::Write(filename,image_to_write);
    else if(extension=="jpg") JPG_FILE<T>::Write(filename,image_to_write);
    else if(extension=="ppm") PPM_FILE<T>::Write(filename,image_to_write);
    else if(extension=="png") PNG_FILE<T>::Write(filename,image_to_write);
    else if(extension=="rgb") RGB_FILE<T>::Write(filename,image_to_write);
    else if(extension=="pbi") FILE_UTILITIES::Write_To_File<float>(filename,image_to_write);
    else PHYSBAM_FATAL_ERROR(LOG::sprintf("Unknown image file extension from filename '%s' extension '%s'",filename.c_str(),extension.c_str()));
    delete corrected_image;
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool IMAGE<T>::
Is_Supported(const std::string& filename)
{
    std::string extension=FILE_UTILITIES::Get_File_Extension(filename);
    if(extension=="bmp") return BMP_FILE<T>::Is_Supported();
    else if(extension=="jpg") return JPG_FILE<T>::Is_Supported();
    else if(extension=="ppm") return PPM_FILE<T>::Is_Supported();
    else if(extension=="png") return PNG_FILE<T>::Is_Supported();
    else if(extension=="rgb") return RGB_FILE<T>::Is_Supported();
    else if(extension=="pbi") return true;
    else return false;
}
//#####################################################################
// Function Byte_Color_To_Scalar_Color
//#####################################################################
template<class T> template<int d> VECTOR<T,d> IMAGE<T>::
Byte_Color_To_Scalar_Color(const VECTOR<unsigned char,d> color_in)
{
    return (VECTOR<T,d>(color_in)+VECTOR<T,d>::All_Ones_Vector())/(T)256;
}
//#####################################################################
// Function Byte_Color_To_Scalar_Color
//#####################################################################
template<class T> T IMAGE<T>::
Byte_Color_To_Scalar_Color(const unsigned char color_in)
{
    return ((T)color_in+(T).5)/256;
}
//#####################################################################
// Function Scalar_Color_To_Byte_Color
//#####################################################################
template<class T> template<int d> VECTOR<unsigned char,d> IMAGE<T>::
Scalar_Color_To_Byte_Color(const VECTOR<T,d> color_in)
{
    return VECTOR<unsigned char,d>(clamp(VECTOR<int,d>(VECTOR<T,d>((T)256*color_in)),VECTOR<int,d>(),VECTOR<int,d>(255*VECTOR<int,d>::All_Ones_Vector())));
}
//#####################################################################
// Function Scalar_Color_To_Byte_Color
//#####################################################################
template<class T> unsigned char IMAGE<T>::
Scalar_Color_To_Byte_Color(const T color_in)
{
    return clamp((int)((T)256*color_in),0,255);
}
//#####################################################################
// Function Flip_X
//#####################################################################
template<class T> template<int d> void IMAGE<T>::
Flip_X(ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{
    for(int i=0;i<image.m/2;i++)for(int j=0;j<image.n;j++)exchange(image(i,j),image(image.m-i,j));
}
//#####################################################################
// Function Flip_Y
//#####################################################################
template<class T> template<int d> void IMAGE<T>::
Flip_Y(ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{
    for(int j=0;j<image.n/2;j++)for(int i=0;i<image.m;i++)exchange(image(i,j),image(i,image.n-j));
}
//#####################################################################
// Function Invert
//#####################################################################
template<class T> template<int d> void IMAGE<T>::
Invert(ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{
    for(int i=image.domain.min_corner.x;i<image.domain.max_corner.x;i++) for(int j=image.domain.min_corner.y;j<image.domain.max_corner.y;j++) image(i,j)=VECTOR<T,d>::All_Ones_Vector()-image(i,j);
}
//#####################################################################
// Function Threshold
//#####################################################################
template<class T> void IMAGE<T>::
Threshold(ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image,const T threshold,const VECTOR<T,3>& low_color,const VECTOR<T,3>& high_color)
{
    for(int i=image.domain.min_corner.x;i<image.domain.max_corner.x;i++) for(int j=image.domain.min_corner.y;j<image.domain.max_corner.y;j++) if(image(i,j).Magnitude()<threshold) image(i,j)=low_color;else image(i,j)=high_color;
}
//#####################################################################
namespace PhysBAM{
template class IMAGE<float>;
template void IMAGE<float>::Read(const std::string&,ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&);
template void IMAGE<float>::Read(const std::string&,ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&);
template void IMAGE<float>::Write(const std::string&,const ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&,const float,const float);
template void IMAGE<float>::Write(const std::string&,const ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&,const float,const float);
template class IMAGE<double>;
template void IMAGE<double>::Read(const std::string&,ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&);
template void IMAGE<double>::Read(const std::string&,ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&);
template void IMAGE<double>::Write(const std::string&,const ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&,const double,const double);
template void IMAGE<double>::Write(const std::string&,const ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&,const double,const double);
template VECTOR<double,3> IMAGE<double>::Byte_Color_To_Scalar_Color<3>(VECTOR<unsigned char,3>);
template VECTOR<float,3> IMAGE<float>::Byte_Color_To_Scalar_Color<3>(VECTOR<unsigned char,3>);
template VECTOR<unsigned char,3> IMAGE<double>::Scalar_Color_To_Byte_Color<3>(VECTOR<double,3>);
template VECTOR<unsigned char,3> IMAGE<float>::Scalar_Color_To_Byte_Color<3>(VECTOR<float,3>);
template VECTOR<unsigned char,4> IMAGE<double>::Scalar_Color_To_Byte_Color<4>(VECTOR<double,4>);
template VECTOR<unsigned char,4> IMAGE<float>::Scalar_Color_To_Byte_Color<4>(VECTOR<float,4>);
template void IMAGE<float>::Invert<3>(ARRAY<VECTOR<float,3>,VECTOR<int,2> >&);
template void IMAGE<double>::Invert<3>(ARRAY<VECTOR<double,3>,VECTOR<int,2> >&);
}
