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
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
class CELL_DOMAIN_INTERFACE:NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    TV_INT a;
    int b;
    ARRAY<int> flat_base; // flat index interface element reference cell (min corner cell of "coarse_factor"-wide block)
    ARRAY<int> remap; // maps ghost cells inside for periodic bc, identity for non-periodic bc
    ARRAY<int> cell_location;
    ARRAY<bool> bdy_element;
    
public:

    const GRID<TV>& grid;
    const int padding;
    const TV_INT size;
    const int flat_size;
    const int coarse_factor;
    const int interface_elements;
    const bool periodic_bc;
    const TV_INT coarse_counts;
    const RANGE<TV_INT> coarse_range;
    
    CELL_DOMAIN_INTERFACE(const GRID<TV>& grid_input,int padding_input,int coarse_factor_input,int interface_elements_input,int periodic_bc_input);

    inline int Flatten(const TV_INT& index) const {return index.Dot(a)+b;}
    inline int Flatten_Diff(const TV_INT& index) const {return index.Dot(a);}
    inline int Get_Flat_Base(int e) const {return flat_base(e);}
    inline int Remap(int i) const {return remap(i);}
    inline int Is_Inside_Cell(int i) const {return cell_location(i)==-1;}
    inline int Is_Outside_Cell(int i) const {return cell_location(i)==1;}
    inline int Is_Boundary_Cell(int i) const {return cell_location(i)==0;}
    inline int Is_Boundary_Element(int i) const {return bdy_element(i);}
    
    void Set_Flat_Base(int start,int end,const TV_INT& index);
    void Initialize();
};

template<class TV>
class CELL_MANAGER:NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:

    const CELL_DOMAIN_INTERFACE<TV>& cdi;
    VECTOR<ARRAY<int>,2> compressed;  // inside and outside; [-1] - inactive; [-2] - active (temporary); [non-negative] - dof number (for active)
    VECTOR<int,2> dofs; // inside and outside; number of dofs

    CELL_MANAGER(const CELL_DOMAIN_INTERFACE<TV>& cdi_input);
    
    inline void Set_Active(int i,int s){compressed[s](i)=-2;compressed[s](cdi.Remap(i))=-2;}
    inline int Get_Index(TV_INT index,int s){return compressed[s](cdi.Flatten(index));}
    
    void Compress_Indices();
};
}
#endif
