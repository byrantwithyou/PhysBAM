//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <fstream>
#include <map>
#include <string>
#include "PARSE_DATA_FEM.h"
namespace PhysBAM{
//#####################################################################
// Function Parse_Input
//#####################################################################
template<class TV> void PARSE_DATA_FEM<TV>::
Parse_Input(const std::string& pipe_file)
{
    std::ifstream fin(pipe_file);
    std::string line;
    half_width=1;

    std::map<std::string,VERTEX_ID> pts_index;
    TV pt;
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
            case 'h':
                ss>>unit_length;
                break;
            case 'v':
                ss>>name1>>pt;
                pts_index[name1]=pts.m;
                pts.Append({pt,nobc,TV(),default_joint,{}});
                break;
            case 'p':
                {
                    ss>>name1>>name2;
                    VERTEX_ID i0=pts_index[name1];
                    VERTEX_ID i1=pts_index[name2];
                    PIPE_ID p=pipes.Append({i0,i1});
                    pts(i0).joints.Append(p);
                    pts(i1).joints.Append(p);
                }
                break;
            case 's':
                {
                    ss>>name1>>pts(pts_index[name1]).bc;
                    pts(pts_index[name1]).bc_type=dirichlet_v;
                }
                break;
            case 't':
                {
                    ss>>name1>>pts(pts_index[name1]).bc;
                    pts(pts_index[name1]).bc_type=traction;
                }
                break;
            case 'j':
                {
                    ss>>name1>>name2;
                    VERTEX_ID i=pts_index[name1];
                    if(name2=="corner")
                        pts(i).joint_type=corner_joint;
                }
                break;
            default:
                LOG::printf("PARSE FAIL: %c %s\n",c,ss.str());
        }
    }
}
template struct PARSE_DATA_FEM<VECTOR<double,2> >;
template struct PARSE_DATA_FEM<VECTOR<double,3> >;
}
