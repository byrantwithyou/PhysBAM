//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMAGE
//#####################################################################
#ifndef __IMAGE__
#define __IMAGE__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/exchange.h>
#include <Core/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class IMAGE
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;

    IMAGE()
    {}

//#####################################################################
    template<int d> static VECTOR<T,d> Byte_Color_To_Scalar_Color(const VECTOR<unsigned char,d> color_in);
    static T Byte_Color_To_Scalar_Color(const unsigned char color_in);
    template<int d> static VECTOR<unsigned char,d> Scalar_Color_To_Byte_Color(const VECTOR<T,d> color_in);
    static unsigned char Scalar_Color_To_Byte_Color(const T color_in);
    template<int d> static void Flip_X(ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image);
    template<int d> static void Flip_Y(ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image);
    template<int d> static void Invert(ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image);
    static void Threshold(ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image,const T threshold,const VECTOR<T,3>& low_color,
        const VECTOR<T,3>& high_color);
    template<int d> static void Read(const std::string& filename,ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image);
    template<int d> static void Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image,const T gamma=1,const T dither_amplitude=0);
    static bool Is_Supported(const std::string& filename);
//#####################################################################
};
}
#endif
