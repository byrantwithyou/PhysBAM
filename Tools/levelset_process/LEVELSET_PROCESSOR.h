//#####################################################################
// Copyright 2003-2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_PROCESSOR
//#####################################################################
#ifndef __LEVELSET_PROCESSOR__
#define __LEVELSET_PROCESSOR__

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
namespace PhysBAM{

template<class T>
class LEVELSET_PROCESSOR
{
public:
    LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset;
    ARRAYS<VECTOR<T,3> >& phi;
    GRID<VECTOR<T,3> >& grid;
    int verbose;

    LEVELSET_PROCESSOR(LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset_input,const int verbose_input=1)
        :levelset(levelset_input),phi(levelset_input.phi),grid(levelset_input.grid),
        verbose(verbose_input)
    {}

//#####################################################################
    void Smooth_By_Curvature(const int iters=1,const T sigma=1);
    template<class T2> void Window_Sum_3D(ARRAYS<VECTOR<T2,3> >& phi,const int xradius,const int yradius,const int zradius);
    void Create_Erosion_Kernel(ARRAYS<VECTOR<int,3> >& kernel,const int xradius,const int yradius,const int zradius);
    void Dilation(const int xradius,const int yradius,const int zradius);
    void Erosion(const int xradius,const int yradius,const int zradius);
    void Closing(const int xradius,const int yradius,const int zradius);
    void Opening(const int xrdius, const int yradius, const int zradius);
    void Smooth(const int xradius,const int yradius,const int zradius);
    void Subtract_Set_Same_Size(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset_input);
    void Subtract_Set(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset_input);
    void Add_Set_Same_Size(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset_input);
    void Add_Set(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset_input);
    void Output_Subset(LEVELSET_3D<GRID<VECTOR<T,3> > >& subset,const int m0,const int m1,const int n0,const int n1,const int mn0,const int mn1) const;
    int Num_Points(const int start,const int end,const int stride,const bool include_end) const;
    void Output_Subsample(LEVELSET_3D<GRID<VECTOR<T,3> > >& subsample,const int x_stride,const int y_stride,const int z_stride,bool include_ends) const;
    void Output_Resample(LEVELSET_3D<GRID<VECTOR<T,3> > >& supersample,const int m,const int n,const int mn) const;
    void Output_Resample(LEVELSET_3D<GRID<VECTOR<T,3> > >& resample) const;
//#####################################################################
};
}
#endif
