//####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __IRRADIANCE_SAMPLE__
#define __IRRADIANCE_SAMPLE__

#include <Core/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class IRRADIANCE_SAMPLE
{
    typedef VECTOR<T,3> TV;
public:
    TV location;
    TV direction;
    TV irradiance;
    T harmonic_mean_distance;
    
    IRRADIANCE_SAMPLE()
    {}

    IRRADIANCE_SAMPLE(const TV& location_input,const TV& direction_input,const TV& irradiance_input,const T harmonic_mean_distance_input)
        :location(location_input),direction(direction_input),irradiance(irradiance_input),harmonic_mean_distance(harmonic_mean_distance_input)
    {}
    
//#####################################################################
};
}
#endif
