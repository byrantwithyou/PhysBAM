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
    ARRAY_VIEW<TV> X;

    POINT_CLOUD(ARRAY_COLLECTION* array_collection_input)
        :array_collection(array_collection_input),X(0,0)
    {Initialize_Array_Collection();}

    POINT_CLOUD()
        :array_collection(new ARRAY_COLLECTION()),X(0,0)
    {Initialize_Array_Collection();}

    virtual ~POINT_CLOUD()
    {delete array_collection;}

    void Initialize_Array_Collection()
    {array_collection->Add_Array(ATTRIBUTE_ID_X,&X);}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& V,const T dt)
    {X.Subset(indices)+=dt*V.Subset(indices);}

    TV Weighted_Center(ARRAY<T> weights)
    {typename TV::SCALAR total=weights.Sum();return total?X.Weighted_Sum(weights)/total:TV();}

    RANGE<TV> Compute_Bounding_Box()
    {return RANGE<TV>::Bounding_Box(X.array);}

    TV Centroid()
    {return X.Average();};

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
