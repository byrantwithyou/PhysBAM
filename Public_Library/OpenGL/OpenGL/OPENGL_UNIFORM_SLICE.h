//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_UNIFORM_SLICE
//##################################################################### 
#ifndef __OPENGL_UNIFORM_SLICE__
#define __OPENGL_UNIFORM_SLICE__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Math_Tools/max.h>
#include <Tools/Math_Tools/min.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <OpenGL/OpenGL/OPENGL_SLICE.h>

namespace PhysBAM{
template<class T> class OPENGL_WORLD;
template<class T>
class OPENGL_UNIFORM_SLICE:public OPENGL_SLICE
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    OPENGL_WORLD<T>& world;
    GRID<TV> grid;

    OPENGL_SLICE::SLICE_MODE mode;

    int axis;
    int index;

    GLenum clip_plane_id1,clip_plane_id2;

    OPENGL_UNIFORM_SLICE(OPENGL_WORLD<T>& world_input)
        : world(world_input),clip_plane_id1(0),clip_plane_id2(0)
    {
        Initialize(GRID<TV>(TV_INT()+2,RANGE<TV>::Unit_Box()));
    }

    bool Is_Slice_Mode() override
    {
        return mode==NODE_SLICE || mode==CELL_SLICE;
    }

    int Maximum_Slice_Index_In_Current_Mode()
    {
        VECTOR<int,3> grid_cells(grid.numbers_of_cells);
        if(mode==NODE_SLICE) return grid_cells[axis]+1;
        else if(mode==CELL_SLICE) return grid_cells[axis];
        else return -1;
    }

    void Initialize(GRID<TV> grid_input)
    {
        grid=grid_input;
        mode=NO_SLICE;
        axis=2;
        index=grid_input.counts.z/2;
        Update_Clip_Planes();
    }

    void Set_Slice_Mode(OPENGL_SLICE::SLICE_MODE mode_input) override
    {
        mode=mode_input;
        if(Is_Slice_Mode()) {
            // keep current values of axis and slice number as long as they're in range
            axis=clamp(axis,0,2);
            index=clamp(index,0,Maximum_Slice_Index_In_Current_Mode());
        }
        Update_Clip_Planes();
    }

    void Toggle_Slice_Mode() override
    {
        Set_Slice_Mode((OPENGL_SLICE::SLICE_MODE)(((int)mode+1)%3));
    }
    
    void Toggle_Slice_Axis() override
    {
        if(Is_Slice_Mode()) {
            axis=(axis+1)%3;
            index=Maximum_Slice_Index_In_Current_Mode()/2+1;
            Update_Clip_Planes();
        }
    }

    void Increment_Slice() override
    {
        if(Is_Slice_Mode()) {
            index=min(index+1,Maximum_Slice_Index_In_Current_Mode());
            Update_Clip_Planes();
        }
    }

    void Decrement_Slice() override
    {
        if(Is_Slice_Mode()) {
            index=max(index-1,1);
            Update_Clip_Planes();
        }
    }
    
    void Enable_Clip_Planes() override
    {
        glEnable(clip_plane_id1);
        glEnable(clip_plane_id2);
    }

    template<class T2>
    static void Get_Face_Index_Range(const OPENGL_UNIFORM_SLICE<T>* slice,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& array,int face,VECTOR<int,3> &index_start,VECTOR<int,3> &index_end,int scale=1)
    {
        index_start=VECTOR<int,3>(array.domain.min_corner.x,array.domain.min_corner.y,array.domain.min_corner.z);
        index_end=VECTOR<int,3>(array.domain.max_corner.x,array.domain.max_corner.y,array.domain.max_corner.z);
        if(!slice) return;
        else if(slice->mode==CELL_SLICE) {
            if(face==slice->axis) {
                index_start[slice->axis]=max(slice->index/scale,index_start[slice->axis]);
                index_end[slice->axis]=min(slice->index/scale+1,index_end[slice->axis]);
            } else {
                index_start[slice->axis]=max(slice->index/scale,index_start[slice->axis]);
                index_end[slice->axis]=min(slice->index/scale,index_end[slice->axis]);
            }
        } else if(slice->mode==NODE_SLICE) {
            if(face==slice->axis) {
                index_start[slice->axis]=max(slice->index/scale,index_start[slice->axis]);
                index_end[slice->axis]=min(slice->index/scale,index_end[slice->axis]);
            } else {
                index_start=VECTOR<int,3>(1,1,1);
                index_end=VECTOR<int,3>(0,0,0);
            }
        }
    }

    GRID<TV> Get_Slice_Grid() const
    {
        GRID<TV> slice_grid(grid.Get_MAC_Grid());
        slice_grid.domain.min_corner(axis)+=index*grid.dX(axis);
        slice_grid.domain.max_corner(axis)=slice_grid.domain.min_corner(axis)+grid.dX(axis);
        slice_grid.counts(axis)=1;
        slice_grid.numbers_of_cells=slice_grid.counts;
        return slice_grid;
    }

    void Print_Slice_Info(std::ostream& output_stream) override;
private:
    void Update_Clip_Planes() override;
};
}
#endif
