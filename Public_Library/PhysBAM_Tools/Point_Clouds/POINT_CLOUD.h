//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_CLOUD
//#####################################################################
#ifndef __POINT_CLOUD__
#define __POINT_CLOUD__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_BASE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <memory>
namespace PhysBAM{

template<class TV>
class POINT_CLOUD:public CLONEABLE<POINT_CLOUD<TV> >
{
    typedef typename TV::SCALAR T;
public:
    typedef T SCALAR;
    typedef TV VECTOR_T;

    ARRAY_COLLECTION* array_collection;

    POINT_CLOUD(ARRAY_COLLECTION* array_collection_input)
        :array_collection(array_collection_input)
    {}

    POINT_CLOUD()
        :array_collection(new ARRAY_COLLECTION())
    {}

    virtual ~POINT_CLOUD()
    {delete array_collection;}

    void Clone_Helper(const POINT_CLOUD<TV>& particles)
    {array_collection->Initialize(*particles.array_collection);}

    bool operator==(const POINT_CLOUD<TV>& particles) const
    {return (this==&particles || *array_collection==*particles.array_collection);}

    bool operator!=(const POINT_CLOUD<TV>& particles) const
    {return !(*this==particles);}
};
}
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#endif
