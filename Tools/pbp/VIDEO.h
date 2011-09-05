//#####################################################################
// Copyright 2006, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VIDEO
//#####################################################################
#ifndef __VIDEO__
#define __VIDEO__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{

template<class T>
class VIDEO_READER
{
public:
    virtual ~VIDEO_READER()
    {}

    virtual bool Valid() const=0;
    virtual int Frame_Rate() const=0;
    virtual int Frame_Count() const=0;
    virtual int Width() const=0;
    virtual int Height() const=0;
    virtual void Frame(const int frame,ARRAY<VECTOR<T,3>,VECTOR<int,2> >& image)=0;
    virtual void Frame(const int frame,unsigned char* image)=0;

    static VIDEO_READER<T>* Read(const std::string& filename);
    static VIDEO_READER<T>* Read(const std::string& format_string,const int min_frame,const int max_frame);
};

template<class T>
class VIDEO_WRITER
{
public:
    virtual ~VIDEO_WRITER()
    {}

    virtual int Frame(ARRAY<VECTOR<T,3>,VECTOR<int,2> >& image)=0;
    virtual bool Valid() const=0;


    static VIDEO_WRITER<T>* Write(const std::string& filename,const std::string& codec,const int bitrate,const int frame_rate,const int width,const int height);
    static bool Enabled();
};

}
#endif
