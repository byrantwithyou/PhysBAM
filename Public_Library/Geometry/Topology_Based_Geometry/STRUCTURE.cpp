//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
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
    if(!structure) throw VALUE_ERROR(LOG::sprintf("%s has no Create() function.",name.c_str()));
    return structure;
}
//#####################################################################
// Function Create_From_Name
//#####################################################################
template<class TV> STRUCTURE<TV>* STRUCTURE<TV>::
Create_From_Name(const std::string& name,GEOMETRY_PARTICLES<TV>& particles)
{
    STRUCTURE* structure=STRUCTURE_REGISTRY<TV>::Name_To_Factory(name)->Create(particles);
    if(!structure) throw VALUE_ERROR(LOG::sprintf("%s has no Create(GEOMETRY_PARTICLES<TV>& particles) function.",name.c_str()));
    return structure;
}
//#####################################################################
// Function Create_From_Fxtension
//#####################################################################
template<class TV> STRUCTURE<TV>* STRUCTURE<TV>::
Create_From_Extension(const std::string& extension)
{
    STRUCTURE* structure=STRUCTURE_REGISTRY<TV>::Extension_To_Factory(extension)->Create();
    if(!structure) throw VALUE_ERROR(LOG::sprintf("No Create() function matching extension %s",extension.c_str()));
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
// Function Update_Number_Nodes
//#####################################################################
template<class TV> void STRUCTURE<TV>::
Update_Number_Nodes()
{
}
//#####################################################################
// Function Mark_Nodes_Referenced
//#####################################################################
template<class TV> void STRUCTURE<TV>::
Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const
{
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> std::string STRUCTURE<TV>::
Name() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return Static_Name();
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV> std::string STRUCTURE<TV>::
Extension() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return Static_Extension();
}
//#####################################################################
// Function Read_Structure
//#####################################################################
template<class TV> void STRUCTURE<TV>::
Read_Structure(TYPED_ISTREAM& input)
{
    std::string name;
    Read_Binary(input,name);
    if(name!=Name()){
        LOG::cerr<<"Trying to read in a "<<name<<" as a "<<Name()<<std::endl;
        PHYSBAM_FATAL_ERROR();}
    Read(input);
}
//#####################################################################
// Function Write_Structure
//#####################################################################
template<class TV> void STRUCTURE<TV>::
Write_Structure(TYPED_OSTREAM& output)
{
    Write_Binary(output,Name());
    Write(output);
}
//#####################################################################
// Function Create_Structure
//#####################################################################
template<class TV> STRUCTURE<TV>*  STRUCTURE<TV>::
Create_Structure(TYPED_ISTREAM& input,GEOMETRY_PARTICLES<TV>& particles)
{
    std::string name;
    Read_Binary(input,name);
    STRUCTURE<TV>* structure=STRUCTURE<TV>::Create_From_Name(name,particles);
    structure->Read(input);
    return structure;
}
//#####################################################################
// Function Create_Structure
//#####################################################################
template<class TV> STRUCTURE<TV>*  STRUCTURE<TV>::
Create_Structure(TYPED_ISTREAM& input)
{
    std::string name;
    Read_Binary(input,name);
    STRUCTURE<TV>* structure=STRUCTURE<TV>::Create_From_Name(name);
    structure->Read(input);
    return structure;
}
//#####################################################################
// Function Create_From_File
//#####################################################################
template<class TV> template<class RW> STRUCTURE<TV>*  STRUCTURE<TV>::
Create_From_File(const std::string& filename)
{
    STRUCTURE<TV>* structure=STRUCTURE<TV>::Create_From_Extension(FILE_UTILITIES::Get_File_Extension(filename));
    FILE_UTILITIES::template Read_From_File<RW>(filename,*structure);
    return structure;
}
//#####################################################################
namespace PhysBAM{
template class STRUCTURE<VECTOR<float,1> >;
template class STRUCTURE<VECTOR<float,2> >;
template class STRUCTURE<VECTOR<float,3> >;
template class STRUCTURE<VECTOR<double,1> >;
template class STRUCTURE<VECTOR<double,2> >;
template class STRUCTURE<VECTOR<double,3> >;
template STRUCTURE<VECTOR<double,1> >* STRUCTURE<VECTOR<double,1> >::Create_From_File<double>(std::string const&);
template STRUCTURE<VECTOR<double,1> >* STRUCTURE<VECTOR<double,1> >::Create_From_File<float>(std::string const&);
template STRUCTURE<VECTOR<double,2> >* STRUCTURE<VECTOR<double,2> >::Create_From_File<double>(std::string const&);
template STRUCTURE<VECTOR<double,2> >* STRUCTURE<VECTOR<double,2> >::Create_From_File<float>(std::string const&);
template STRUCTURE<VECTOR<double,3> >* STRUCTURE<VECTOR<double,3> >::Create_From_File<double>(std::string const&);
template STRUCTURE<VECTOR<double,3> >* STRUCTURE<VECTOR<double,3> >::Create_From_File<float>(std::string const&);
template STRUCTURE<VECTOR<float,1> >* STRUCTURE<VECTOR<float,1> >::Create_From_File<double>(std::string const&);
template STRUCTURE<VECTOR<float,1> >* STRUCTURE<VECTOR<float,1> >::Create_From_File<float>(std::string const&);
template STRUCTURE<VECTOR<float,2> >* STRUCTURE<VECTOR<float,2> >::Create_From_File<double>(std::string const&);
template STRUCTURE<VECTOR<float,2> >* STRUCTURE<VECTOR<float,2> >::Create_From_File<float>(std::string const&);
template STRUCTURE<VECTOR<float,3> >* STRUCTURE<VECTOR<float,3> >::Create_From_File<double>(std::string const&);
template STRUCTURE<VECTOR<float,3> >* STRUCTURE<VECTOR<float,3> >::Create_From_File<float>(std::string const&);
}
