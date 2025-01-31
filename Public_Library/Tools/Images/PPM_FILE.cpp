//#####################################################################
// Copyright 2002-2005, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Tools/Images/IMAGE.h>
#include <Tools/Images/PPM_FILE.h>
using namespace PhysBAM;
//#####################################################################
// Function Read
//#####################################################################
template<class T> void PPM_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("PPM_FILE::Read not implemented");
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void PPM_FILE<T>::
Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image)
{
    PHYSBAM_FATAL_ERROR("PPM_FILE::Read not implemented");
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<int d> void PPM_FILE<T>::
Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
{  
    VECTOR<int,2> counts=image.domain.Edge_Lengths();
    assert(image.domain.min_corner.x==0 && image.domain.min_corner.y==0);
    std::ostream* output=Safe_Open_Output_Raw(filename,true,false); // no compression
    *output<<"P6"<<std::endl;
    *output<<"# Generated using PPM_FILE::Write"<<std::endl;
    *output<<counts.x<<" "<<counts.y<<std::endl;
    *output<<255<<std::endl;
    for(int j=counts.y-1;j>=0;j--)for(int i=0;i<counts.x;i++){VECTOR<unsigned char,d> pixel=IMAGE<T>::Scalar_Color_To_Byte_Color(image(i,j));Write_Binary<T>(*output,pixel[0],pixel[1],pixel[2]);}
    delete output;
}
//#####################################################################
// Function Is_Supported
//#####################################################################
template<class T> bool PPM_FILE<T>::
Is_Supported()
{
    return true;
}
//#####################################################################
namespace PhysBAM{
template class PPM_FILE<float>;
template void PPM_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,3> ,VECTOR<int,2> >&);
template void PPM_FILE<float>::Write(const std::string&,const ARRAY<VECTOR<float,4> ,VECTOR<int,2> >&);
template class PPM_FILE<double>;
template void PPM_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,3> ,VECTOR<int,2> >&);
template void PPM_FILE<double>::Write(const std::string&,const ARRAY<VECTOR<double,4> ,VECTOR<int,2> >&);
}
