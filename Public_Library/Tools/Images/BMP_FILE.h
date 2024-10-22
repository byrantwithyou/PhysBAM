//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BMP_FILE 
//#####################################################################
#ifndef __BMP_FILE__
#define __BMP_FILE__

#include <Core/Vectors/VECTOR_FORWARD.h>
#include <Grid_Tools/Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <string>
namespace PhysBAM{

template<class T>
class BMP_FILE
{
public:
    BMP_FILE()
    {}
    
//#####################################################################
    static void Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image);
    static void Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image);
    template<int d> static void Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image);
    static bool Is_Supported();
//#####################################################################
};
}
#endif
