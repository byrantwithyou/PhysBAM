//#####################################################################
// Copyright 2004-2006, Jiayi Chong, Igor Neverov, Andrew Selle, Mike Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE
//#####################################################################
#ifndef __SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE__
#define __SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE__

#include <Core/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE
{
public:
    VECTOR<T,3> position;
    VECTOR<T,3> normal;
    VECTOR<T,3> transmitted_irradiance;
    VECTOR<T,3> transmitted_irradiance_product;
    T area;

    SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE()
        :area(0)
    {}

    SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE(const VECTOR<T,3>& position_input,const VECTOR<T,3>& normal_input,const T area_input)
        :position(position_input),normal(normal_input),area(area_input)
    {}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,position,normal,transmitted_irradiance,transmitted_irradiance_product,area);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,position,normal,transmitted_irradiance,transmitted_irradiance_product,area);}

//#####################################################################
};
}
#endif
