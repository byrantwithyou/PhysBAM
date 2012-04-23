//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_MANAGER
//#####################################################################
#ifndef __CELL_MANAGER__
#define __CELL_MANAGER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
class CELL_DOMAIN_INTERFACE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    TV_INT a;
    int b;
    ARRAY<int> flat_base;

public:

    const GRID<TV>& grid;
    const TV_INT size;
    const TV_INT coarse_range;
    const int padding;
    const int flat_size;
    const int coarse_factor;
    const int interface_elements;
    
    CELL_DOMAIN_INTERFACE(const GRID<TV>& grid_input,int padding_input,int coarse_factor_input,int interface_elements_input):
        grid(grid_input),padding(padding_input),size(grid.counts+2*padding),flat_size(size.Product()),coarse_factor(coarse_factor_input),
        coarse_range(TV_INT()+coarse_factor),interface_elements(interface_elements_input)
    {
        a(TV::m-1)=1;
        for(int i=TV::m-2;i>=0;i--) a(i)=a(i+1)*size(i+1);
        b=(TV_INT()+padding).Dot(a);
        flat_base.Resize(interface_elements_input);
    }

    inline int Flatten(const TV_INT& index){return index.Dot(a)+b;}
    inline int Flatten_Diff(const TV_INT& index){return index.Dot(a);}
    inline int Get_Flat_Base(int e){return flat_base(e);}
    
    void Set_Flat_Base(int start,int end,const TV_INT& index)
    {int flat=Flatten(index);for(int i=start;i<end;i++) flat_base(i)=flat;}
};

template<class TV>
class CELL_MANAGER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:

    const CELL_DOMAIN_INTERFACE<TV>& cdi;
    ARRAY<bool> active_cells[2];  // inside and outside

    CELL_MANAGER(const CELL_DOMAIN_INTERFACE<TV>& cdi_input):cdi(cdi_input)
    {for(int s=0;s<2;s++) active_cells[s].Resize(cdi.flat_size);}
};
}
#endif
