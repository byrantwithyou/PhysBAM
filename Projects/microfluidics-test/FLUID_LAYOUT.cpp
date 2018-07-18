//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <map>
#include "FLUID_LAYOUT.h"

namespace PhysBAM{
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Compute(const PARSE_DATA<TV>& pd)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    used_faces.Resize(grid,1,true,true,{nodof,-1,-1,0});
    used_cells.Resize(grid.Domain_Indices(1),use_init,{nodof,-1,-1,0,-1});
    
    for(auto& i:pd.pipes)
    {
        RANGE<TV_INT> box=pd.Pipe_Full_Range(i);
        for(FACE_RANGE_ITERATOR<TV::m> it(box,RF::ghost);it.Valid();it.Next())
            used_faces(it.face).type=wall;
        for(RANGE_ITERATOR<TV::m> it(box,RI::ghost);it.Valid();it.Next())
            used_cells(it.index).type=wall;
    }

    for(auto& i:pd.pts)
    {
        for(FACE_RANGE_ITERATOR<TV::m> it(i.box,RF::ghost);it.Valid();it.Next())
            used_faces(it.face).type=wall;
        for(FACE_RANGE_ITERATOR<TV::m> it(i.box,RF::skip_outer);it.Valid();it.Next())
            used_faces(it.face).type=fluid;
        for(RANGE_ITERATOR<TV::m> it(i.box);it.Valid();it.Next())
            used_cells(it.index).type=fluid;
    }
    
    for(auto& i:pd.pipes)
    {
        RANGE<TV_INT> box=pd.Pipe_Full_Range(i);
        for(FACE_RANGE_ITERATOR<TV::m> it(box,RF::skip_outer);it.Valid();it.Next())
            used_faces(it.face).type=fluid;
        for(RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next())
            used_cells(it.index).type=fluid;
    }

    int next_block=0;
    for(auto& i:pd.pts)
    {
        RANGE<TV_INT> box=i.box;
        if(i.bc_type==dirichlet)
        {
            for(RANGE_ITERATOR<TV::m> it(box,1,0,RI::ghost|RI::omit_corners,i.bc_side);it.Valid();it.Next())
            {
                used_cells(it.index)={dirichlet,next_block,-1,i.bc_value,-1};
            }
            for(FACE_RANGE_ITERATOR<TV::m> it(box,0,0,RF::ghost,i.bc_side);it.Valid();it.Next())
            {
                used_faces(it.face)={fluid,next_block,-1,0,-1};
            }
        }
        if(i.bc_type==wall)
        {
            for(FACE_RANGE_ITERATOR<TV::m> it(box,0,0,RF::ghost,i.bc_side);it.Valid();it.Next())
            {
                used_faces(it.face).type=i.bc_type;
                used_faces(it.face).bc_value=i.bc_value;
            }
        }
        for(RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next())
        {
            used_cells(it.index).block_id=next_block;
        }
        for(FACE_RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next())
        {
            auto& uf=used_faces(it.face);
            if(uf.type==fluid) uf.block_id=next_block;
        }
        next_block++;
        blocks.Append({0});
    }
    for(auto& i:pd.pipes)
    {
        int dir=pd.Pipe_Dir(i);
        TV_INT dpt=pd.pts(i.y).pt-pd.pts(i.x).pt;
        int num_blocks=dpt.Max()-2*pd.half_width;
        TV_INT pt=pd.pts(i.x).pt;
        RANGE<TV_INT> box=pd.Pipe_Inner_Range(i);
        for(RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next())
        {
            TV_INT da=it.index-pd.pts(i.x).pt;
            int diff=-1;
            for(int i=0;i<TV::m;i++)
                if(da(i)>=pd.half_width)
                    diff=da(i)-pd.half_width;
            PHYSBAM_ASSERT(diff>=0);
            used_cells(it.index).block_id=diff+next_block;
        }
        for(FACE_RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next())
        {
            TV_INT cell=it.face.index;
            TV_INT da=cell-pd.pts(i.x).pt;
            if(it.face.axis==dir)
                for(int i=0;i<TV::m;i++)
                    if(da(i)>=0 && da(i)<pd.half_width){
                        cell(dir)--;
                        break;}
            used_faces(it.face).block_id=used_cells(cell).block_id;
        }
        for(int i=0;i<num_blocks;i++) blocks.Append({0});
        next_block+=num_blocks;
    }
    for(FACE_RANGE_ITERATOR<TV::m> it(used_faces.domain_indices);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        if(f.type==fluid)
        {
            f.block_dof=blocks(f.block_id).num_dofs++;
            f.global_id=dof_map.Append({f.block_id,f.block_dof});
        }
    }
    for(RANGE_ITERATOR<TV::m> it(used_cells.domain);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        if(c.type==fluid)
        {
            c.block_dof=blocks(c.block_id).num_dofs++;
            c.global_id=dof_map.Append({c.block_id,c.block_dof});
        }
    }
}
//#####################################################################
// Function Dump_Layout
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Dump_Layout() const
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        if(f.type==fluid) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));
        else if(f.type==wall) Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,1,1));
        else PHYSBAM_ASSERT(f.type==nodof);
    }
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        if(c.type==fluid) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));
        else if(c.type==dirichlet) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,1));
        else PHYSBAM_ASSERT(c.type==nodof);
    }
}
//#####################################################################
// Function Dump_Dofs
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Dump_Dofs() const
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        std::string s=LOG::sprintf("%i",f.global_id);
        if(f.type==fluid) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
        else if(f.type==wall) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(0,1,1));
        else Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,0,0));
    }
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        std::string s=LOG::sprintf("%i",c.global_id);
        if(c.type==fluid) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
        else if(c.type==dirichlet) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,0,1));
        else Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,0,0));
    }
}
//#####################################################################
// Function Dump_Blocks
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Dump_Blocks() const
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        if(f.type!=fluid)  continue;
        std::string s=LOG::sprintf("%i",f.block_id);
        Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
    }
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        if(c.type!=fluid)  continue;
        std::string s=LOG::sprintf("%i",c.block_id);
        Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
    }
}
template struct FLUID_LAYOUT<VECTOR<double,2> >;
}
