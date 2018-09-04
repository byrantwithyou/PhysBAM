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
    BC_FUNC bcf;
    bcf.type=dirichlet_v;
    bcf.flowrate=0;
    bc.Append(bcf);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARSE_DATA_FEM<TV>::
~PARSE_DATA_FEM()
{
    delete force;
    delete analytic_velocity;
    delete analytic_pressure;
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV PARSE_DATA_FEM<TV>::
Velocity(const TV& X,BC_ID bc_id) const
{
    if(analytic_velocity) return analytic_velocity->v(X,0);
    else if(!bc(bc_id).flowrate) return TV();
    else{
        T a=half_width*unit_length;
        T r=(X-bc(bc_id).pt).Magnitude();
        T v=-3.0/4*bc(bc_id).flowrate/(a*a*a)*(r*r-a*a);
        return v*bc(bc_id).dir;}
}
//#####################################################################
// Function Traction
//#####################################################################
template<class TV> TV PARSE_DATA_FEM<TV>::
Traction(const TV& X,const TV& N,T mu,BC_ID bc_id) const
{
    if(analytic_velocity && analytic_pressure){
        SYMMETRIC_MATRIX<T,2> stress=analytic_velocity->dX(X,0).Twice_Symmetric_Part()*mu;
        stress-=analytic_pressure->f(X,0);
        return stress*N;}
    else return bc(bc_id).traction;
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> TV PARSE_DATA_FEM<TV>::
Force(const TV& X,T mu) const
{
    if(analytic_velocity && analytic_pressure){
        SYMMETRIC_TENSOR<T,0,TV::m> ddU=analytic_velocity->ddX(X,0);
        TV f=analytic_pressure->dX(X,0);
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
    PHYSBAM_ASSERT(analytic_velocity);
    return analytic_velocity->dX(X,0).Trace();
}
//#####################################################################
// Function Pressure
//#####################################################################
template<class TV> typename TV::SCALAR PARSE_DATA_FEM<TV>::
Pressure(const TV& X) const
{
    PHYSBAM_ASSERT(analytic_pressure);
    return analytic_pressure->f(X,0);
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
                    BC_FUNC bcf;
                    bcf.type=dirichlet_v;
                    ss>>name1>>bcf.flowrate;
                    VERTEX_ID self=pts_index[name1];
                    bcf.pt=pts(self).pt;
                    PHYSBAM_ASSERT(pts(self).joints.m==1);
                    PIPE_ID pid=pts(self).joints(0);
                    VERTEX_ID p0=pipes(pid)(0),p1=pipes(pid)(1);
                    if(p0!=self) std::swap(p0,p1);
                    bcf.dir=(pts(p1).pt-pts(p0).pt).Normalized();
                    pts(pts_index[name1]).bc_id=bc.Append(bcf);
                }
                break;
            case 't':
                {
                    BC_FUNC bcf;
                    bcf.type=traction;
                    ss>>name1>>bcf.traction;
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
            case 'f':
                delete force;
                force=new ANALYTIC_VECTOR_PROGRAM<TV>(ss.str().c_str()+2);
                break;
            case 'U':
                delete analytic_velocity;
                analytic_velocity=new ANALYTIC_VECTOR_PROGRAM<TV>(ss.str().c_str()+2);
                break;
            case 'P':
                delete analytic_pressure;
                analytic_pressure=new ANALYTIC_SCALAR_PROGRAM<TV>(ss.str().c_str()+2);
                break;
            default:
                LOG::printf("PARSE FAIL: %c %s\n",c,ss.str());
        }
    }
}
template struct PARSE_DATA_FEM<VECTOR<double,2> >;
}
