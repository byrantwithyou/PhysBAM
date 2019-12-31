//#####################################################################
// Copyright 2019, Craig Schroeder, Yunxin Sun.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace OBJECT_PLACEMENT
//#####################################################################
#ifndef __OBJECT_PLACEMENT__
#define __OBJECT_PLACEMENT__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Vectors/TWIST.h>

namespace PhysBAM{

template<class T> class RANDOM_NUMBERS;
template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class OBJECT_PLACEMENT
{
public:
    typedef typename TV::SCALAR T;

    RANDOM_NUMBERS<T>& random;
    T min_sep=0;
    int max_intersect_tries=10; // intersects something
    int max_seed_tries=100; // find point in seed volume
    IMPLICIT_OBJECT<TV>* io; // seeding volume
    RANGE<TV> seed_box;

    // User updates these; new objects must not intersect these
    ARRAY<TV> pts;
    ARRAY<RANGE<TV> > boxes;

    struct OBJECT
    {
        RANGE<TV> box;
        INTERVAL<T> scale;
        T max_angulare_vel;
        RANGE<TV> velocity_range;
    };
    ARRAY<OBJECT> objects;

    struct SEEDING
    {
        int index;
        FRAME<TV> frame;
        TWIST<TV> twist;
        T scale;
    };

    OBJECT_PLACEMENT(RANDOM_NUMBERS<T>& random,IMPLICIT_OBJECT<TV>* io);
    ~OBJECT_PLACEMENT()=default;

    bool Seed_Object(SEEDING& obj);
    void Seed_Objects(ARRAY<SEEDING>& objs,int max_objects);
    bool Find_Seed_Location(FRAME<TV>& frame,T& scale,int& index);
    bool Test_Location(const FRAME<TV>& frame,T scale,int index);
};
//#####################################################################
}
#endif
