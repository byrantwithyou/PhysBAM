//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NODE_ITERATOR
//#####################################################################
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
using namespace PhysBAM;
//#####################################################################
// Function NODE_ITERATOR
//#####################################################################
template<class TV> NODE_ITERATOR<TV>::
NODE_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells,const T_REGION& region_type,const int side)
    :GRID_ITERATOR_BASE<TV>(grid_input)
{
    assert(-1<=side&&side<2*TV::dimension);assert(region_type!=GRID<TV>::BOUNDARY_INTERIOR_REGION);
    TV_INT counts=grid.Numbers_Of_Nodes(); // exact number of nodes (even if MAC)
    RANGE<TV_INT> domain(grid.Node_Indices(number_of_ghost_cells));
    switch(region_type){
        case GRID<TV>::WHOLE_REGION: assert(side<0);Add_Region(domain);break;
        case GRID<TV>::GHOST_REGION: assert(number_of_ghost_cells>0); // ghost region of grid with specified ghost cells
            // TODO(jontg): counts doesn't take into account number_of_ghost_cells, so I believe this to be incorrect.
            if(side<0){ // don't loop over the same cell twice!
                TV_INT max_copy(domain.max_corner);
                for(int axis=TV::dimension-1;axis>=0;axis--){
                    domain.max_corner(axis)=0;
                    Add_Region(domain);
                    domain.max_corner(axis)=max_copy(axis);
                    domain.min_corner(axis)=counts(axis);
                    Add_Region(domain);
                    domain.max_corner(axis)=counts(axis);
                    domain.min_corner(axis)=0;}}
            else{int axis=side/2;if(side&1) domain.min_corner(axis)=counts(axis)+1;else domain.max_corner(axis)=0;Add_Region(domain);}
            break;
        case GRID<TV>::BOUNDARY_REGION: // outer boundary of grid with specified ghost cells
            if(side<0){ // don't loop over the same node twice!
                RANGE<TV_INT> domain_copy(domain);
                for(int axis=TV::dimension-1;axis>=0;axis--){
                    domain.max_corner(axis)=domain.min_corner(axis);
                    Add_Region(domain);
                    domain.max_corner(axis)=domain.min_corner(axis)=domain_copy.max_corner(axis);
                    Add_Region(domain);
                    domain.max_corner(axis)=domain_copy.max_corner(axis)-1;
                    domain.min_corner(axis)=domain_copy.min_corner(axis)+1;}}
            else{int axis=side/2;if(side&1) domain.min_corner(axis)=domain.max_corner(axis);else domain.max_corner(axis)=domain.min_corner(axis);Add_Region(domain);}
            break;
        default: assert(region_type==GRID<TV>::INTERIOR_REGION && side<0);domain.Change_Size(-1);Add_Region(domain);break;}
    Reset();
}
//#####################################################################
// Function NODE_ITERATOR
//#####################################################################
template<class TV> NODE_ITERATOR<TV>::
NODE_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input)
    :GRID_ITERATOR_BASE<TV>(grid_input,region_input)
{
}
//#####################################################################
namespace PhysBAM{
template class NODE_ITERATOR<VECTOR<float,1> >;
template class NODE_ITERATOR<VECTOR<float,2> >;
template class NODE_ITERATOR<VECTOR<float,3> >;
template class NODE_ITERATOR<VECTOR<double,1> >;
template class NODE_ITERATOR<VECTOR<double,2> >;
template class NODE_ITERATOR<VECTOR<double,3> >;
}
