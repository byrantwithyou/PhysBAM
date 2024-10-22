//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Tools/Symbolics/PROGRAM.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
#include <Grid_Tools/Grids/GRID.h>
#include "PARSE_DATA.h"
#include <fstream>
#include <map>
#include <string>
namespace PhysBAM{
//#####################################################################
// Function Inflow_BC_Value
//#####################################################################
template<class TV> typename TV::SCALAR PARSE_DATA<TV>::
Inflow_BC_Value(const TV& X,const typename PARSE_DATA<TV>::VERTEX_DATA& vd,const GRID<TV>& grid) const
{
    T a=half_width*grid.dX(1-vd.bc_side/2);
    TV center=grid.Node(vd.pt);
    T sign=vd.bc_side%2?1:-1;
    center(vd.bc_side/2)=center(vd.bc_side/2)+sign*a;
    T r=(X-center).Magnitude();
    T v=-3.0/4*vd.bc_value/(a*a*a)*(r*r-a*a);
    return v;
}
//#####################################################################
// Function Parse_Input
//#####################################################################
template<class TV> void PARSE_DATA<TV>::
Parse_Input(const std::string& pipe_file)
{
    typedef VECTOR<int,TV::m> TV_INT;
    std::ifstream fin(pipe_file);
    std::string line;
    half_width=1;

    RANGE<TV_INT> pts_box;
    std::map<std::string,int> pts_index;
    std::map<std::string,INDEX_TYPE> pts_bc;
    TV_INT pt;
    int a;
    std::string name1,name2;
    char c;
    while(getline(fin,line))
    {
        std::stringstream ss(line);
        ss>>c;
        if(isspace(c) || c=='#') continue;
        switch(c)
        {
            case 'w':
                ss>>half_width;
                break;
            case 'v':
                ss>>name1>>pt;
                pts_box.Enlarge_To_Include_Point(pt);
                pts_index[name1]=pts.m;
                pts.Append({pt,fluid,-1,0,RANGE<TV_INT>(pt-half_width,pt+half_width)});
                break;
            case 'p':
                {
                    ss>>name1>>name2;
                    int i0=pts_index[name1];
                    int i1=pts_index[name2];
                    TV_INT d=pts(i1).pt-pts(i0).pt;
                    if(d.Sum()<0) std::swap(i0,i1);
                    TV_INT size=abs(d).Sorted();
                    assert(size.x==0 && size.y>half_width*2);
                    pipes.Append({i0,i1});
                }
                break;
            case 's':
                {
                    char* end=0;
                    ss>>name1>>a>>name2;
                    int i0=pts_index[name1];
                    pts(i0).bc_type=wall;
                    pts(i0).bc_side=a;
                    pts(i0).bc_value=strtod(name2.c_str(),&end);
                    if(*end)
                    {
                    }
                }
                break;
            case 'd':
                {
                    char* end=0;
                    ss>>name1>>a>>name2;
                    int i0=pts_index[name1];
                    pts(i0).bc_type=dirichlet;
                    pts(i0).bc_side=a;
                    pts(i0).bc_value=strtod(name2.c_str(),&end);
                }
                break;
            default:
                LOG::printf("PARSE FAIL: %c %s\n",c,ss.str());
        }
    }

    pts_box=pts_box.Thickened(half_width);
    for(auto& i:pts)
    {
        i.pt-=pts_box.min_corner;
        i.box-=pts_box.min_corner;
    }
    box_size=pts_box.Edge_Lengths();
}
//#####################################################################
// Function Pipe_Full_Range
//#####################################################################
template<class TV> auto PARSE_DATA<TV>::
Pipe_Full_Range(const VECTOR<int,2>& p) const -> RANGE<TV_INT>
{
    return RANGE<TV_INT>(pts(p.x).box.min_corner,pts(p.y).box.max_corner);
}
//#####################################################################
// Function Pipe_Inner_Range
//#####################################################################
template<class TV> auto PARSE_DATA<TV>::
Pipe_Inner_Range(const VECTOR<int,2>& p) const -> RANGE<TV_INT>
{
    int dir=Pipe_Dir(p);
    RANGE<TV_INT> A=pts(p.x).box;
    RANGE<TV_INT> B=pts(p.y).box;
    RANGE<TV_INT> box(A.min_corner,B.max_corner);
    box.min_corner(dir)=A.max_corner(dir);
    box.max_corner(dir)=B.min_corner(dir);
    return box;
}
//#####################################################################
// Function Pipe_Dir
//#####################################################################
template<class TV> int PARSE_DATA<TV>::
Pipe_Dir(const VECTOR<int,2>& p) const
{
    return (pts(p.y).pt-pts(p.x).pt).Arg_Max();
}
template struct PARSE_DATA<VECTOR<double,2> >;
template struct PARSE_DATA<VECTOR<double,3> >;
}
