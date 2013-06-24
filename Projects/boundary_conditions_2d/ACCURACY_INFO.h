#ifndef __ACCURACY_INFO__
#define __ACCURACY_INFO__
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Vectors/VECTOR.h>
using namespace PhysBAM;

template<int d>
struct ACCURACY_INFO
{
    typedef VECTOR<int,d> TV_INT;
    int base_resolution,resolution;
    RANGE<TV_INT> sample_domain;
    ARRAY<TV_INT> cell_samples;
    ARRAY<FACE_INDEX<d> > face_samples;

    void Compute();
    template<class T> void Print(const char* label,ARRAY<T,TV_INT>& a) const;
    template<class T> void Print(const char* label,ARRAY<T,FACE_INDEX<d> >& a) const;
    template<class TV> void Print_Locations(const GRID<TV>& grid) const;
};

#endif
