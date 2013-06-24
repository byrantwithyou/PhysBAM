//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CUSTOM_PHI
//#####################################################################
//
//#####################################################################
// Fedkiw - August 8, 2001
//#####################################################################
#ifndef __BOUNDARY_CUSTOM_PHI__
#define __BOUNDARY_CUSTOM_PHI__

#include <Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{

class BOUNDARY_CUSTOM_PHI:public BOUNDARY
{
public:
    BOUNDARY_CUSTOM_PHI() 
                                            :BOUNDARY()
    {}

//#####################################################################
    void Fill_Ghost_Cells(GRID_2D& grid,ARRAY_2D& u,ARRAY_2D& u_ghost,const double dt=0,const double time=0,const int number_of_ghost_cells=3) const;
    void Apply_Boundary_Condition(GRID_2D& grid,ARRAY_2D& u,const double time=0) const{} // do nothing
//#####################################################################
};
}
#endif
