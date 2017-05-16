//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MAC_GRID_PERIODIC
//#####################################################################
#ifndef __BOUNDARY_MAC_GRID_PERIODIC__
#define __BOUNDARY_MAC_GRID_PERIODIC__

#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_MAC_GRID_PERIODIC:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    using BOUNDARY<TV,T2>::Find_Ghost_Regions;
    VECTOR<bool,TV::m> is_periodic;

    BOUNDARY_MAC_GRID_PERIODIC()
    {is_periodic.Fill(true);}

    ~BOUNDARY_MAC_GRID_PERIODIC() = default;

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T2,FACE_INDEX<TV::m> >& u,ARRAY<T2,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const override {} // do nothing
    void Apply_Boundary_Condition_Face(const GRID<TV>& grid,ARRAY<T2,FACE_INDEX<TV::m> >& u,const T time) const override;
//#####################################################################
};
template<class T2,class TV> void
Fill_Ghost_Cells_Periodic(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& u,
    ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& u_ghost,const VECTOR<bool,TV::m>& is_periodic,
    const int number_of_ghost_cells);
template<class T2,class TV> void
Fill_Ghost_Faces_Periodic(const GRID<TV>& grid,const ARRAY<T2,FACE_INDEX<TV::m> >& u,
    ARRAY<T2,FACE_INDEX<TV::m> >& u_ghost,const VECTOR<bool,TV::m>& is_periodic,
    const int number_of_ghost_cells);
template<class T2,class TV> void
Apply_Boundary_Condition_Face_Periodic(const GRID<TV>& grid,ARRAY<T2,FACE_INDEX<TV::m> >& u,
    const VECTOR<bool,TV::m>& is_periodic);
template<class T2,class TV> void
Fill_Ghost_Cells_Periodic_Accum(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& u,
    ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& u_ghost,const VECTOR<bool,TV::m>& is_periodic,
    const int number_of_ghost_cells);
template<class T2,class TV> void
Fill_Ghost_Faces_Periodic_Accum(const GRID<TV>& grid,const ARRAY<T2,FACE_INDEX<TV::m> >& u,
    ARRAY<T2,FACE_INDEX<TV::m> >& u_ghost,const VECTOR<bool,TV::m>& is_periodic,
    const int number_of_ghost_cells);

}
#endif
