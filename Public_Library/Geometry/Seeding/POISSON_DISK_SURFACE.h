//#####################################################################
// Copyright 2020, Craig Schroeder, Yunxin Sun.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_DISK_SURFACE
//#####################################################################
#ifndef __POISSON_DISK_SURFACE__
#define __POISSON_DISK_SURFACE__
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV> class GRID;
template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class POISSON_DISK_SURFACE
{
public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    T min_distance=0; // Must set before using (Set_Distance)
    int max_attemps=30;
    int max_proj_attemps=5;
    int ghost=3;
    T h=0;

    POISSON_DISK_SURFACE()=default;
    ~POISSON_DISK_SURFACE()=default;
    void Sample(RANDOM_NUMBERS<T>& random,IMPLICIT_OBJECT<TV>* io,ARRAY<TV>& X);
    void Set_Distance(T distance);
    void Set_Distance_By_Volume(T volume_per_sample);
private:
    TV Generate_Random_Point_Around_Annulus(RANDOM_NUMBERS<T>& random,const TV& center) const;
    bool Check_Distance(const GRID<TV>& grid,ARRAY<int,TV_INT>& grid_array,const TV& point,ARRAY<TV>& X) const;
    bool Adjust_To_Surface(IMPLICIT_OBJECT<TV>* io,const TV& seed,TV& pt) const;
//#####################################################################
};

}
#endif
