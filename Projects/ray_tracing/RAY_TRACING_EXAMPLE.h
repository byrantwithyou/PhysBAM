//#####################################################################
// Copyright 2002-2005. Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RAY_TRACING_EXAMPLE  
//##################################################################### 
#ifndef __RAY_TRACING_EXAMPLE__
#define __RAY_TRACING_EXAMPLE__

#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T> class RENDER_WORLD;

template<class T>
class RAY_TRACING_EXAMPLE
{
public:
    int frame;
    std::string output_filename,alpha_filename;
    bool keep_old_files;
    int pixels_between_output_directory_existence_checks;//if 0, don't check at all
    bool use_spatial_partition;
    T gamma_correction;
    RANGE<VECTOR<int,2> > clipping_region;

    RAY_TRACING_EXAMPLE()
        :frame(1),output_filename("rgb"),keep_old_files(false),
        pixels_between_output_directory_existence_checks(0),
        use_spatial_partition(true),
        gamma_correction((T)2.2),clipping_region(RANGE<VECTOR<int,2> >::Centered_Box()*10000)
    {}

    virtual ~RAY_TRACING_EXAMPLE()
    {}

    std::string Get_Output_Filename(const int frame)
    {std::string filename;
    filename=LOG::sprintf(output_filename.c_str(),frame);
    if(keep_old_files && File_Exists(filename)){
        std::string test_filename;
        for(int i=1;;i++){
            test_filename=Get_Basename(filename)+"-"+Number_To_String(i)+"."+Get_File_Extension(filename);
            if(!File_Exists(test_filename)) break;}
        filename=test_filename;}
    return filename;}

    std::string Get_Alpha_Filename(const int frame)
    {std::string filename;
    filename=LOG::sprintf(alpha_filename.c_str(),frame);
    if(keep_old_files && File_Exists(filename)){
        std::string test_filename;
        for(int i=1;;i++){
            test_filename=Get_Basename(filename)+"-"+Number_To_String(i)+"."+Get_File_Extension(filename);
            if(!File_Exists(test_filename)) break;}
        filename=test_filename;}
    return filename;}

//#####################################################################
    virtual void Initialize_Scene(RENDER_WORLD<T>& world,const int frame){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
