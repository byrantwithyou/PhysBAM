//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRUCTURE
//##################################################################### 
#ifndef __STRUCTURE__
#define __STRUCTURE__

#include <Core/Arrays/ARRAY.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Read_Write/TYPED_STREAM.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
#include <string>
#include <typeinfo>
namespace PhysBAM{

template<class TV> class GEOMETRY_PARTICLES;


template<class TV>
class STRUCTURE
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
    STRUCTURE(const STRUCTURE&) = delete;
    void operator=(const STRUCTURE&) = delete;
    virtual ~STRUCTURE();

    virtual void Read(TYPED_ISTREAM input)=0;

    virtual void Write(TYPED_OSTREAM output) const=0;

//#####################################################################
    void Read_Structure(TYPED_ISTREAM input);
    void Write_Structure(TYPED_OSTREAM output);
    static STRUCTURE<TV>* Create_Structure(TYPED_ISTREAM input,GEOMETRY_PARTICLES<TV>& particles);
    static STRUCTURE<TV>* Create_Structure(TYPED_ISTREAM input);
    static STRUCTURE<TV>* Create_From_File(const std::string& filename);
    virtual std::string Name() const;
    virtual std::string Extension() const;
    static std::string Static_Name() {return "";}
    static std::string Static_Extension() {return "";}
    static STRUCTURE* Create_From_Name(const std::string& name);
    static STRUCTURE* Create_From_Name(const std::string& name,GEOMETRY_PARTICLES<TV>& particles);
    static STRUCTURE* Create_From_Extension(const std::string& extension);
    virtual void Rescale(const T scaling_factor);
    virtual STRUCTURE* Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& particles,ARRAY<int>* particle_indices=0) const;
    virtual void Update_Number_Nodes();
    virtual void Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const;
//#####################################################################
};
}
#endif
