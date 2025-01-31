//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PNG_FILE (requires libpng)
//#####################################################################
#ifndef __PNG_FILE__
#define __PNG_FILE__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
#include <string>
namespace PhysBAM{

template<class T>
class PNG_FILE
{
public:
    PNG_FILE()
    {}
    
//#####################################################################
    static void Read(const std::string& filename,ARRAY<VECTOR<T,3>,VECTOR<int,2> >& image);
    static void Read(const std::string& filename,ARRAY<VECTOR<T,4>,VECTOR<int,2> >& image);
    template<int d> static void Write(const std::string& filename,const ARRAY<VECTOR<T,d>,VECTOR<int,2> >& image);
    static bool Is_Supported();
//#####################################################################
};
}
#endif
