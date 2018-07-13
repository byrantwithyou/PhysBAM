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
    used_faces.Resize(grid,1,true,true,{noslip,-1,-1,0});
    used_cells.Resize(grid.Domain_Indices(1),use_init,{dirichlet,-1,-1,0,-1});
    
    for(auto& i:pd.pts)
    {
        for(FACE_RANGE_ITERATOR<TV::m> it(i.box,RF::skip_outer);it.Valid();it.Next())
            used_faces(it.face).type=fluid;
        for(RANGE_ITERATOR<TV::m> it(i.box);it.Valid();it.Next())
            used_cells(it.index).type=fluid;
    }
    
    for(auto& i:pd.pipes)
    {
        RANGE<TV_INT> box=RANGE<TV_INT>::Combine(pd.pts(i.x).box,pd.pts(i.y).box);
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
            for(RANGE_ITERATOR<TV::m> it(box,1,0,RI::ghost,i.bc_side);it.Valid();it.Next())
            {
                used_cells(it.index)={i.bc_type,next_block,-1,i.bc_value,-1};
            }
            for(FACE_RANGE_ITERATOR<TV::m> it(box,0,0,RF::ghost,i.bc_side);it.Valid();it.Next())
            {
                used_faces(it.face).type=fluid;
            }
        }
        if(i.bc_type==source)
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
        next_block++;
        int block_type=1+16*i.connected_sides;
        if(i.bc_type==dirichlet) block_type+=2*i.bc_side;
        blocks.Append({0,block_type});
    }
    for(auto& i:pd.pipes)
    {
        TV_INT dpt=abs(pd.pts(i.x).pt-pd.pts(i.y).pt);
        int num_blocks=dpt.Max()-2*pd.half_width;
        TV_INT pt=pd.pts(i.x).pt;
        RANGE<TV_INT> box=RANGE<TV_INT>::Combine(pd.pts(i.x).box,pd.pts(i.y).box);
        for(RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next())
        {
            if(pd.pts(i.x).box.Lazy_Inside_Half_Open(it.index)) continue;
            if(pd.pts(i.y).box.Lazy_Inside_Half_Open(it.index)) continue;
            TV_INT da=it.index-pd.pts(i.x).pt;
            int diff=-1;
            for(int i=0;i<TV::m;i++)
            {
                if(da(i)>=pd.half_width) diff=da(i)-pd.half_width;
                if(da(i)<-pd.half_width) diff=-da(i)-pd.half_width-1;
            }
            PHYSBAM_ASSERT(diff>=0);
            used_cells(it.index).block_id=diff+next_block;
        }
        int block_type=2*dpt.Arg_Max();
        for(int i=0;i<num_blocks;i++) blocks.Append({0,block_type});
        next_block+=num_blocks;
    }
    for(FACE_RANGE_ITERATOR<TV::m> it(used_faces.domain_indices);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        if(f.type==fluid)
        {
            f.block_id=used_cells(it.face.index).block_id;
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

    next_block=0;
    std::map<int,int> block_rename;
    for(int i=0;i<blocks.m;i++)
    {
        blocks(i).is_rep=false;
        auto it=block_rename.find(blocks(i).block_type);
        if(it!=block_rename.end()) blocks(i).block_type=it->second;
        else
        {
            block_rename[blocks(i).block_type]=next_block;
            blocks(i).block_type=next_block++;
            blocks(i).is_rep=true;
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
        else if(f.type==noslip) Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,1,1));
        else if(f.type==source) Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,1,0));
        else Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
    }
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        if(c.type==fluid) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));
        else if(c.type==dirichlet) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,1));
        else Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
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
        else if(f.type==noslip) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(0,1,1));
        else if(f.type==source) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(0,1,0));
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
//#####################################################################
// Function Dump_Block_Types
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Dump_Block_Types() const
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        if(f.type!=fluid)  continue;
        std::string s=LOG::sprintf("%i",blocks(f.block_id).block_type);
        Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
    }
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        if(c.type!=fluid)  continue;
        std::string s=LOG::sprintf("%i",blocks(c.block_id).block_type);
        Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
    }
}
template struct FLUID_LAYOUT<VECTOR<double,2> >;
}
