//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_VECTOR.h>
#include <fstream>
#include <map>
#include <string>
#include "PARSE_DATA_FEM.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARSE_DATA_FEM<TV>::
PARSE_DATA_FEM()
{
    bc.Append({dirichlet_v,0,0,0});
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARSE_DATA_FEM<TV>::
~PARSE_DATA_FEM()
{
    delete force;
    for(BC_ID i(0);i<bc.m;i++){
        delete bc(i).velocity;
        delete bc(i).pressure;
        delete bc(i).traction;}
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV PARSE_DATA_FEM<TV>::
Velocity(const TV& X,BC_ID bc_id) const
{
    BC_ID id=analytic_bc!=BC_ID(-1)?analytic_bc:bc_id;
    PHYSBAM_ASSERT(id!=BC_ID(-1));
    if(bc(id).velocity) return bc(id).velocity->v(X,0);
    else return TV();
}
//#####################################################################
// Function Traction
//#####################################################################
template<class TV> TV PARSE_DATA_FEM<TV>::
Traction(const TV& X,const TV& N,T mu,BC_ID bc_id) const
{
    if(analytic_bc!=BC_ID(-1)){
        SYMMETRIC_MATRIX<T,2> stress=bc(analytic_bc).velocity->dX(X,0).Twice_Symmetric_Part()*mu;
        stress-=bc(analytic_bc).pressure->f(X,0);
        return stress*N;}
    else return bc(bc_id).traction->v(X,0);
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> TV PARSE_DATA_FEM<TV>::
Force(const TV& X,T mu) const
{
    if(analytic_bc!=BC_ID(-1)){
        SYMMETRIC_TENSOR<T,0,TV::m> ddU=bc(analytic_bc).velocity->ddX(X,0);
        TV f=bc(analytic_bc).pressure->dX(X,0);
        f-=mu*(Contract<1,2>(ddU)+Contract<0,2>(ddU));
        return f;}
    else return force?force->v(X,0):TV();
}
//#####################################################################
// Function Divergence
//#####################################################################
template<class TV> typename TV::SCALAR PARSE_DATA_FEM<TV>::
Divergence(const TV& X) const
{
    PHYSBAM_ASSERT(analytic_bc!=BC_ID(-1));
    return bc(analytic_bc).velocity->dX(X,0).Trace();
}
//#####################################################################
// Function Pressure
//#####################################################################
template<class TV> typename TV::SCALAR PARSE_DATA_FEM<TV>::
Pressure(const TV& X) const
{
    PHYSBAM_ASSERT(analytic_bc!=BC_ID(-1));
    return bc(analytic_bc).pressure->f(X,0);
}
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
                {
                    VERTEX_DATA vd;
                    ss>>name1>>vd.pt;
                    pts_index[name1]=pts.m;
                    pts.Append(vd);
                }
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
                    std::string expr;
                    ss>>name1>>expr;
                    BC_FUNC bcf;
                    bcf.type=dirichlet_v;
                    bcf.velocity=new ANALYTIC_VECTOR_PROGRAM<TV>(expr.c_str());
                    pts(pts_index[name1]).bc_id=bc.Append(bcf);
                }
                break;
            case 't':
                {
                    std::string expr;
                    ss>>name1>>expr;
                    BC_FUNC bcf;
                    bcf.type=traction;
                    bcf.traction=new ANALYTIC_VECTOR_PROGRAM<TV>(expr.c_str());
                    pts(pts_index[name1]).bc_id=bc.Append(bcf);
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
            case 'U':
                if(analytic_bc==BC_ID(-1))
                    analytic_bc=bc.Append({analytic,0,0,0});
                else delete bc(analytic_bc).velocity;
                bc(analytic_bc).velocity=new ANALYTIC_VECTOR_PROGRAM<TV>(ss.str().c_str()+2);
                break;
            case 'P':
                if(analytic_bc==BC_ID(-1))
                    analytic_bc=bc.Append({analytic,0,0,0});
                else delete bc(analytic_bc).pressure;
                bc(analytic_bc).pressure=new ANALYTIC_SCALAR_PROGRAM<TV>(ss.str().c_str()+2);
                break;
            case 'F':
                delete force;
                force=new ANALYTIC_VECTOR_PROGRAM<TV>(ss.str().c_str()+2);
                break;
            default:
                LOG::printf("PARSE FAIL: %c %s\n",c,ss.str());
        }
    }
}
template struct PARSE_DATA_FEM<VECTOR<double,2> >;
}
