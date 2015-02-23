//#####################################################################
// Copyright 2015, Andre Pradhana, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_DISK
//#####################################################################
#ifndef __POISSON_DISK__
#define __POISSON_DISK__
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV> class GRID;
template<class TV> class IMPLICIT_OBJECT;;

template<class TV>
class POISSON_DISK
{
public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const T min_distance;
    const int max_attemps;
    const int ghost;
    T h;
    RANDOM_NUMBERS<T> ran;
    POISSON_DISK(const T min_distance,const int max_attemps=30);
    ~POISSON_DISK();
    void Sample(IMPLICIT_OBJECT<TV>* object,ARRAY<TV>& XXX);
private:
    TV Generate_Random_Point_Around_Annulus(RANDOM_NUMBERS<T>& random,TV& center) const;
    bool Check_Distance(const GRID<TV>& grid,ARRAY<int,TV_INT>& grid_array,const TV& point,ARRAY<TV>& XX) const;
//#####################################################################
};

}
#endif
