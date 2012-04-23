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
    ARRAY<int> remap; // maps ghost cells inside for periodic bc, identity for non-periodic bc
    
public:

    const GRID<TV>& grid;
    const TV_INT size;
    const TV_INT coarse_range;
    const int padding;
    const int flat_size;
    const int coarse_factor;
    const int interface_elements;
    const bool periodic_bc;
    
    CELL_DOMAIN_INTERFACE(const GRID<TV>& grid_input,int padding_input,int coarse_factor_input,int interface_elements_input,int periodic_bc_input):
        grid(grid_input),padding(padding_input),size(grid.counts+2*padding),flat_size(size.Product()),coarse_factor(coarse_factor_input),
        coarse_range(TV_INT()+coarse_factor),interface_elements(interface_elements_input),periodic_bc(periodic_bc_input)
    {
        a(TV::m-1)=1;
        for(int i=TV::m-2;i>=0;i--) a(i)=a(i+1)*size(i+1);
        b=(TV_INT()+padding).Dot(a);
        flat_base.Resize(interface_elements_input);
        Initialize_Remap();
    }

    inline int Flatten(const TV_INT& index){return index.Dot(a)+b;}
    inline int Flatten_Diff(const TV_INT& index){return index.Dot(a);}
    inline int Get_Flat_Base(int e){return flat_base(e);}
    inline int Remap(int i){return remap(i);}
    
    void Set_Flat_Base(int start,int end,const TV_INT& index)
    {int flat=Flatten(index);for(int i=start;i<end;i++) flat_base(i)=flat;}

    void Initialize_Remap()
    {
        remap.Resize(flat_size);
        if(periodic_bc){
            for(int axis=0;axis<TV::m;axis++)
                for(int s=0;s<2;s++){
                    int side=axis*2+s;
                    int sign=s?1:-1;
                    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,padding,GRID<TV>::GHOST_REGION,side);it.Valid();it.Next())
                        remap(Flatten(it.index))=Flatten(it.index+sign*grid.counts(axis)*TV_INT::Axis_Vector(axis));}
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
                int i=Flatten(it.index);
                remap(i)=i;}}
        else for(int i=0;i<flat_size;i++) remap(i)=i;
    }
};

template<class TV>
class CELL_MANAGER
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:

    const CELL_DOMAIN_INTERFACE<TV>& cdi;
    ARRAY<int> cell_mapping[2];  // inside and outside; [-1] - inactive; [-2] - active (temporary); [non-negative] - dof number (for active)
    int dofs[2];

    CELL_MANAGER(const CELL_DOMAIN_INTERFACE<TV>& cdi_input);
    
    inline void Set_Active(int i,int s){cell_mapping[s](i)=-2;cell_mapping[s](cdi->Remap(i))=-2;}
    
    void Create_Cell_Mapping();
};
}
#endif
