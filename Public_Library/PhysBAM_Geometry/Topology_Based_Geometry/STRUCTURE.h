//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRUCTURE
//##################################################################### 
#ifndef __STRUCTURE__
#define __STRUCTURE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/TYPED_STREAM.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <string>
#include <typeinfo>
namespace PhysBAM{

template<class TV> class GEOMETRY_PARTICLES;


template<class TV>
class STRUCTURE:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    typedef int HAS_TYPED_READ_WRITE;
    typedef T SCALAR;
    typedef TV VECTOR_T;

    bool update_every_frame;

protected:
    STRUCTURE();
public:
    virtual ~STRUCTURE();

    virtual void Read(TYPED_ISTREAM& input)=0;

    void Read_Structure(TYPED_ISTREAM& input)
    {std::string name;Read_Binary(input,name);
    if(name!=Name()){LOG::cerr<<"Trying to read in a "<<name<<" as a "<<Name()<<std::endl;PHYSBAM_FATAL_ERROR();}
    Read(input);}

    virtual void Write(TYPED_OSTREAM& output) const=0;

    void Write_Structure(TYPED_OSTREAM& output)
    {Write_Binary(output,Name());Write(output);}

    static STRUCTURE<TV>* Create_Structure(TYPED_ISTREAM& input,GEOMETRY_PARTICLES<TV>& particles)
    {std::string name;Read_Binary(input,name);
    STRUCTURE<TV>* structure=STRUCTURE<TV>::Create_From_Name(name,particles);
    structure->Read(input);return structure;}

    template<class RW>
    static STRUCTURE<TV>* Create_From_File(const std::string& filename)
    {STRUCTURE<TV>* structure=STRUCTURE<TV>::Create_From_Extension(FILE_UTILITIES::Get_File_Extension(filename));
    FILE_UTILITIES::template Read_From_File<RW>(filename,*structure);return structure;}

//#####################################################################
    virtual std::string Name() const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return Static_Name();}
    virtual std::string Extension() const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return Static_Extension();}
    static std::string Static_Name() {return "";}
    static std::string Static_Extension() {return "";}
    static STRUCTURE* Create_From_Name(const std::string& name);
    static STRUCTURE* Create_From_Name(const std::string& name,GEOMETRY_PARTICLES<TV>& particles);
    static STRUCTURE* Create_From_Extension(const std::string& extension);
    virtual void Rescale(const T scaling_factor);
    virtual STRUCTURE* Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& particles,ARRAY<int>* particle_indices=0) const;
    virtual void Update_Number_Nodes(){}
    virtual void Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const {}
//#####################################################################
};
}
#endif
