//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Tools/Symbolics/PROGRAM.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
#include <fstream>
#include <map>
#include <string>
#include "PARSE_DATA.h"
namespace PhysBAM{

// template<class TV>
// void Set_BC_Func(const std::string& expr,PARSE_DATA<TV>& pd,typename PARSE_DATA<TV>::VERTEX_DATA& vd)
// {
//     typedef typename TV::SCALAR T;
//     PROGRAM<T> prog;
//     PROGRAM_CONTEXT<T> context;
//     for(int i=0;i<TV::m-1;i++) prog.var_in.Append("rs"[i]);
//     prog.var_out.Append("v");
//     prog.Parse(expr.c_str(),false);
//     prog.Optimize();
//     prog.Finalize();
//     context.Initialize(prog);
    
//     vd.bc_func_value=[prog,context,&pd](VECTOR<T,TV::m-1> x)
//         {
//             for(int i=0;i<TV::m-1;i++) context.data_in(i)=x(i)/(2*pd.half_width)+(T).5;
//             prog.Execute(context);
//             return context.data_out(0);
//         };
// }
// //#####################################################################
// // Function Value
// //#####################################################################
// template<class TV> auto PARSE_DATA<TV>::VERTEX_DATA::
// Value(const VECTOR<int,TV::m>& c) const -> T
// {
//     if(bc_func_value) return bc_func_value(c.Remove_Index(bc_side/2)-pt.Remove_Index(bc_side/2)+(T).5);
//     return bc_value;
// }
// //#####################################################################
// // Function Value
// //#####################################################################
// template<class TV> auto PARSE_DATA<TV>::VERTEX_DATA::
// Value(const FACE_INDEX<TV::m>& f) const -> T
// {
//     if(bc_func_value) return bc_func_value(f.index.Remove_Index(bc_side/2)-pt.Remove_Index(bc_side/2)+(T).5);
//     return bc_value;
// }
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
                    TV_INT d=pts(i0).pt-pts(i1).pt;
                    for(int i=0;i<TV::m;i++)
                    {
                        if(d(i)<0)
                        {
                            pts(i0).connected_sides|=1<<(2*i+1);
                            pts(i1).connected_sides|=1<<(2*i);
                        }
                        if(d(i)>0)
                        {
                            pts(i0).connected_sides|=1<<(2*i);
                            pts(i1).connected_sides|=1<<(2*i+1);
                        }
                    }
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
                    pts(i0).bc_type=source;
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
template struct PARSE_DATA<VECTOR<double,2> >;
}
