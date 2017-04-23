//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL_MPI
//#####################################################################
#ifndef __FLOOD_FILL_MPI__
#define __FLOOD_FILL_MPI__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Parallel_Computation/MPI_GRID.h>
namespace PhysBAM{
template<class T> class MPI_UNIFORM_GRID;

template<class TV>
class FLOOD_FILL_MPI
{
    typedef typename TV::SCALAR T;
    typedef MPI_UNIFORM_GRID<TV> T_MPI_GRID;typedef VECTOR<int,TV::m> TV_INT;
    typedef TV_INT T_INDEX;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
    typedef GRID<TV> T_PARALLEL_GRID;
public:
    const T_MPI_GRID& mpi_grid;
    const GRID<TV>& local_grid;
    const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N;
    int& number_of_regions;
    T_ARRAYS_INT& colors;
    ARRAY<ARRAY<int> >& color_ranks;
    ARRAY<bool>* color_touches_uncolorable;

    FLOOD_FILL_MPI(const T_MPI_GRID& mpi_grid_input,const GRID<TV>& local_grid_input,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N_input,int& number_of_regions_input,T_ARRAYS_INT& colors_input,
        ARRAY<ARRAY<int> >& color_ranks_input,ARRAY<bool>* color_touches_uncolorable);
    FLOOD_FILL_MPI(const FLOOD_FILL_MPI&) = delete;
    void operator=(const FLOOD_FILL_MPI&) = delete;
    virtual ~FLOOD_FILL_MPI();

//#####################################################################
    int Synchronize_Colors();
protected:
    void Find_Global_Colors(ARRAY<bool,VECTOR<int,1> >& color_is_global,const RANGE<TV_INT>&) const;
    template<class T_BOX_HORIZONTAL_INT> void Find_Global_Colors(ARRAY<bool,VECTOR<int,1> >& color_is_global,const T_BOX_HORIZONTAL_INT&) const;
    struct Find_Global_Colors_Helper{template<class T_FACE> static void Apply(const T_MPI_GRID& mpi_grid,const GRID<TV>& local_grid,const ARRAY<bool>& psi_N,const ARRAY<int>& colors,
        ARRAY<bool,VECTOR<int,1> >& color_is_global);};
    void Translate_Local_Colors_To_Global_Colors(const ARRAY<int,VECTOR<int,1> >& color_map,T_ARRAYS_INT& colors_copy,const RANGE<TV_INT>& region,const int global_color_offset) const;
    template<class T_BOX_HORIZONTAL_INT> void Translate_Local_Colors_To_Global_Colors(const ARRAY<int,VECTOR<int,1> >& color_map,ARRAY<int>& colors_copy,const T_BOX_HORIZONTAL_INT& region,
        const int global_color_offset) const;
    void Find_Color_Matches(const ARRAY<int,VECTOR<int,1> >& color_map,UNION_FIND<>& union_find,T_ARRAYS_INT& colors_copy,const RANGE<TV_INT>& region,const int global_color_offset) const;
    template<class T_BOX_HORIZONTAL_INT> void Find_Color_Matches(const ARRAY<int,VECTOR<int,1> >& color_map,UNION_FIND<>& union_find,ARRAY<int>& colors_copy,const T_BOX_HORIZONTAL_INT& region,
        const int global_color_offset) const;
    void Remap_Colors(ARRAY<int,VECTOR<int,1> >& color_map,const RANGE<TV_INT>&);
    template<class T_BOX_HORIZONTAL_INT> void Remap_Colors(ARRAY<int,VECTOR<int,1> >& color_map,const T_BOX_HORIZONTAL_INT&);
//#####################################################################
};
}
#endif
