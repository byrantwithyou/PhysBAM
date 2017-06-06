//#####################################################################
// Copyright 2015, Andre Pradhana, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_DISK
//#####################################################################
#ifndef __POISSON_DISK__
#define __POISSON_DISK__
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV> class GRID;
template<class TV> class IMPLICIT_OBJECT;;

template<class TV>
class POISSON_DISK
{
public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    T min_distance;
    int max_attemps;
    int ghost;
    T h;
    VECTOR<bool,TV::m> is_periodic;

    POISSON_DISK(T min_distance,int max_attemps=30);
    ~POISSON_DISK();
    void Sample(RANDOM_NUMBERS<T>& random,IMPLICIT_OBJECT<TV>& object,ARRAY<TV>& X);
    void Sample(RANDOM_NUMBERS<T>& random,const RANGE<TV>& box,ARRAY<TV>& X);
    void Set_Distance_By_Volume(T volume_per_sample);
private:
    TV Generate_Random_Point_Around_Annulus(RANDOM_NUMBERS<T>& random,TV& center) const;
    bool Check_Distance(const GRID<TV>& grid,ARRAY<int,TV_INT>& grid_array,const TV& point,ARRAY<TV>& X) const;
//#####################################################################
};

}
#endif
