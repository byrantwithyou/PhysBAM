//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLES
//#####################################################################
#ifndef __PARTICLES__
#define __PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_BASE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Point_Clouds/PARTICLES_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <memory>
namespace PhysBAM{

template<class TV>
class PARTICLES:public CLONEABLE<PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
public:
    typedef T SCALAR;
    typedef TV VECTOR_T;
    typedef int HAS_UNTYPED_READ_WRITE;

    ARRAY_COLLECTION* array_collection;

    PARTICLES(ARRAY_COLLECTION* array_collection_input)
        :array_collection(array_collection_input)
    {}

    PARTICLES()
        :array_collection(new ARRAY_COLLECTION())
    {}

    virtual ~PARTICLES()
    {delete array_collection;}

    void Clone_Helper(const PARTICLES<TV>& particles)
    {array_collection->Initialize(*particles.array_collection);}

    bool operator==(const PARTICLES<TV>& particles) const
    {return (this==&particles || *array_collection==*particles.array_collection);}

    bool operator!=(const PARTICLES<TV>& particles) const
    {return !(*this==particles);}

    template<class RW> void Read(std::istream& input)
    {int version;
    Read_Binary<RW>(input,version);
    if(version!=1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized particle version %d",(int)version));
    Read_Binary<RW>(input,*array_collection);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,1,*array_collection);}

    void Print(std::ostream& output,const int p)
    {Read_Write<ARRAY_COLLECTION,float>::Print(output,*array_collection,p);}
};
}
#endif
