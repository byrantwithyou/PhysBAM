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

    std::map<std::string,int> pts_index;
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
                joints.Set(pts.m,ARRAY<int>());
                pts.Append({pt,nobc,default_joint});
                break;
            case 'p':
                {
                    ss>>name1>>name2;
                    int i0=pts_index[name1];
                    int i1=pts_index[name2];
                    int p=pipes.Append({i0,i1});
                    if(ARRAY<int>* r=joints.Get_Pointer(i0))
                        r->Append(p);
                    else PHYSBAM_ASSERT(false);
                    if(ARRAY<int>* r=joints.Get_Pointer(i1))
                        r->Append(p);
                    else PHYSBAM_ASSERT(false);
                }
                break;
            case 's':
                {
                    ss>>name1;
                    pts(pts_index[name1]).bc_type=dirichlet_v;
                }
                break;
            case 't':
                {
                    ss>>name1;
                    pts(pts_index[name1]).bc_type=traction;
                }
                break;
            default:
                LOG::printf("PARSE FAIL: %c %s\n",c,ss.str());
        }
    }

    for(auto& i:joints)
    {
        if(i.data.m==1)
            pts(i.key).joint_type=end_vertex;
        else if(i.data.m==2 && pts(i.key).joint_type==default_joint)
            pts(i.key).joint_type=arc;
    }
}
template struct PARSE_DATA_FEM<VECTOR<double,2> >;
template struct PARSE_DATA_FEM<VECTOR<double,3> >;
}
