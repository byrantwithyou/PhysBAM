//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ARRAYS_UTILITIES__
#define __ARRAYS_UTILITIES__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV,class T2>
class ARRAYS_UTILITIES
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;typedef VECTOR<T2,TV::m> TV_T2;
    typedef ARRAY<T2,TV_INT> T_ARRAYS_DIMENSION_T2;
    typedef ARRAY<TV_T2,TV_INT> T_ARRAYS_DIMENSION_VECTOR_T2;
    typedef ARRAY<T2,FACE_INDEX<TV::m> > T_FACE_ARRAYS_T2;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
public:

//#####################################################################
    static void Make_Ghost_Mask_From_Active_Mask(const GRID<TV>& grid,const ARRAY<bool,TV_INT>& input_mask,ARRAY<bool,TV_INT>& output_mask,const int stencil_width,const int number_of_ghost_cells=0);
    static void Compute_Face_Data_From_Cell_Data(const GRID<TV>& face_grid,T_FACE_ARRAYS_T2& face_array,const T_ARRAYS_DIMENSION_T2& cell_array,const int number_of_ghost_cells=0);
    static void Compute_Gradient_At_Faces_From_Cell_Data(const GRID<TV>& face_grid,T_FACE_ARRAYS_T2& grad_face_array,const T_ARRAYS_DIMENSION_T2& cell_array,const int number_of_ghost_cells=0);
    static void Compute_Gradient_At_Cells_From_Face_Data(const GRID<TV>& face_grid,T_ARRAYS_DIMENSION_VECTOR_T2& grad_cell_array,const T_FACE_ARRAYS_T2& face_array,const int number_of_ghost_cells=0);
    static void Compute_Gradient_At_Cells_From_Cell_Data(const GRID<TV>& grid,T_ARRAYS_DIMENSION_VECTOR_T2& grad_cell_array,const T_ARRAYS_DIMENSION_T2& cell_array,const int number_of_ghost_cells);
    static void Compute_Divergence_At_Cells_From_Face_Data(const GRID<TV>& face_grid,T_ARRAYS_DIMENSION_T2& div_cell_array,const T_FACE_ARRAYS_T2& face_array,const int number_of_ghost_cells=0);
//#####################################################################
};
}
#endif
