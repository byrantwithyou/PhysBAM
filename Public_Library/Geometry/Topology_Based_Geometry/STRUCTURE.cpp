//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
template<class TV> STRUCTURE<TV>& Representative(const std::string& name);
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STRUCTURE<TV>::
STRUCTURE()
    :update_every_frame(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STRUCTURE<TV>::
~STRUCTURE()
{}
//#####################################################################
// Function Create_From_Name
//#####################################################################
template<class TV> STRUCTURE<TV>* STRUCTURE<TV>::
Create_From_Name(const std::string& name)
{
    STRUCTURE* structure=STRUCTURE_REGISTRY<TV>::Name_To_Factory(name)->Create();
    if(!structure) throw VALUE_ERROR(STRING_UTILITIES::string_sprintf("%s has no Create() function.",name.c_str()));
    return structure;
}
//#####################################################################
// Function Create_From_Name
//#####################################################################
template<class TV> STRUCTURE<TV>* STRUCTURE<TV>::
Create_From_Name(const std::string& name,GEOMETRY_PARTICLES<TV>& particles)
{
    STRUCTURE* structure=STRUCTURE_REGISTRY<TV>::Name_To_Factory(name)->Create(particles);
    if(!structure) throw VALUE_ERROR(STRING_UTILITIES::string_sprintf("%s has no Create(GEOMETRY_PARTICLES<TV>& particles) function.",name.c_str()));
    return structure;
}
//#####################################################################
// Function Create_From_Fxtension
//#####################################################################
template<class TV> STRUCTURE<TV>* STRUCTURE<TV>::
Create_From_Extension(const std::string& extension)
{
    STRUCTURE* structure=STRUCTURE_REGISTRY<TV>::Extension_To_Factory(extension)->Create();
    if(!structure) throw VALUE_ERROR(STRING_UTILITIES::string_sprintf("No Create() function matching extension %s",extension.c_str()));
    return structure;
}
//#####################################################################
// Function Rescale
//#####################################################################
template<class TV> void STRUCTURE<TV>::
Rescale(const T scaling_factor)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV> STRUCTURE<TV>* STRUCTURE<TV>::
Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& particles,ARRAY<int>* particle_indices) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
namespace PhysBAM{
template class STRUCTURE<VECTOR<float,1> >;
template class STRUCTURE<VECTOR<float,2> >;
template class STRUCTURE<VECTOR<float,3> >;
template class STRUCTURE<VECTOR<double,1> >;
template class STRUCTURE<VECTOR<double,2> >;
template class STRUCTURE<VECTOR<double,3> >;
}
