//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Projection/CUT_CELL_PROJECTION.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,TV::m> TV_INT;

int main(int argc, char* argv[])
{
    CUT_CELL_PROJECTION<TV> proj;
    // Negative means inside fluid for both surface and object
    std::string surface_str="p=-1;",object_str="p=-1",u_str="u=0;v=0;",p_str="p=0;",bc_str="pppp";
    int ghost=3;
    RANGE<TV> domain=RANGE<TV>::Centered_Box();
    TV_INT size;
    int resolution=32;
    T density=100;
    T dt=.01;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-sf",&surface_str,"prog","free surface level set");
    parse_args.Add("-ob",&object_str,"prog","object level set");
    parse_args.Add("-u",&u_str,"prog","solution velocity");
    parse_args.Add("-p",&p_str,"prog","solution pressure");
    parse_args.Add("-bc",&bc_str,"prog","boundary conditions [fsp]");
    parse_args.Add("-print_matrix",&proj.print_matrix,"print system matrix and vectors");
    parse_args.Add("-print_residual",&proj.print_residual,"print residual after linear solve");
    parse_args.Add_Not("-no_precon",&proj.use_preconditioner,"disable preconditioner");
    parse_args.Add("-test_system",&proj.test_system,"print system matrix and vectors");
    parse_args.Add("-solver_tolerance",&proj.solver_tolerance,"tol","Krylov tolerance");
    parse_args.Add("-solver_iterations",&proj.solver_iterations,"num","Krylov iterations");
    parse_args.Add("-resolution",&resolution,"num","grid resolution");
    parse_args.Add("-lo",&domain.min_corner,"vec","grid domain lower corner");
    parse_args.Add("-hi",&domain.max_corner,"vec","grid domain upper corner");
    parse_args.Add("-rho",&density,"value","fluid density");
    parse_args.Add("-dt",&dt,"value","time step size");
    parse_args.Parse();

    GRID<TV> grid(TV_INT()+resolution,domain,true);
    GRID<TV> phi_grid=grid.Get_Regular_Grid().Get_MAC_Grid_At_Regular_Positions();
    VIEWER_DIR viewer_dir("proj");
    VIEWER_OUTPUT vo(STREAM_TYPE((RW)0),viewer_dir);
    Use_Debug_Particles<TV>();
    vo.Add_Common("grid",grid);

    ARRAY<T,FACE_INDEX<TV::m> > u(grid.Domain_Indices(ghost));
    vo.Add("mac_velocities",u);
    ARRAY<bool,FACE_INDEX<TV::m> > valid_u(grid,ghost);
    ARRAY<T,TV_INT> object_phi(grid.Node_Indices(ghost));
    ARRAY<T,TV_INT> surface_phi(grid.Node_Indices(ghost));
    ARRAY<T,TV_INT> pressure(grid.Domain_Indices(ghost));
    ARRAY<bool,TV_INT> valid_p(grid.Domain_Indices(ghost));
    proj.object_phi=&object_phi;
    proj.surface_phi=&surface_phi;
    proj.valid_u=&valid_u;
    proj.valid_p=&valid_p;
    proj.pressure=&pressure;

    ANALYTIC_VECTOR_PROGRAM<TV> av(u_str);
    ANALYTIC_SCALAR_PROGRAM<TV> ap(p_str);
    ANALYTIC_SCALAR_PROGRAM<TV> as(surface_str);
    ANALYTIC_SCALAR_PROGRAM<TV> ao(object_str);
    proj.bc_v=[&av](const TV& X){return av.v(X,0);};
    proj.bc_p=[&ap](const TV& X){return ap.f(X,0);};

    for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        object_phi(it.index)=ao.f(it.Location(),0);
        surface_phi(it.index)=as.f(it.Location(),0);}

    for(FACE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        TV X=it.Location();
        TV V=av.v(X,0),dp=ap.dX(X,0);
        u(it.face)=(V+dt*dp/density)(it.face.axis);}

    for(int i=0;i<2*TV::m;i++){
        if(bc_str[i]=='p') proj.bc_type(i)=0;
        else if(bc_str[i]=='f') proj.bc_type(i)=1;
        else if(bc_str[i]=='s') proj.bc_type(i)=2;
        else PHYSBAM_FATAL_ERROR();}

    Dump_Levelset(phi_grid,object_phi,VECTOR<T,3>(1,1,0));
    Dump_Levelset(phi_grid,surface_phi,VECTOR<T,3>(1,0,1));
    Flush_Frame("before projection");

    proj.Cut_Cell_Projection(grid,ghost,u,density,dt);

    Dump_Levelset(phi_grid,object_phi,VECTOR<T,3>(1,1,0));
    Dump_Levelset(phi_grid,surface_phi,VECTOR<T,3>(1,0,1));
    Flush_Frame("after projection");
    
    T u_l2=0,u_inf=0,u_max=0;
    int u_cnt=0;
    for(FACE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        if(!valid_u(it.face)) continue;
        TV X=it.Location();
        T v=av.v(X,0)(it.face.axis),w=u(it.face),e=std::abs(v-w);
        Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>("display_size",e);
        u_inf=std::max(u_inf,e);
        u_max=std::max(u_max,std::abs(v));
        u_l2+=e*e;
        u_cnt++;}
    u_l2=sqrt(u_l2/u_cnt);
    LOG::printf("u: %.3f %.3f   %.3f\n",u_inf,u_l2,u_max);
    Dump_Levelset(phi_grid,object_phi,VECTOR<T,3>(1,1,0));
    Dump_Levelset(phi_grid,surface_phi,VECTOR<T,3>(1,0,1));
    Flush_Frame("velocity error");

    T p_l2=0,p_inf=0,p_max=0;
    int p_cnt=0;
    for(CELL_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()){
        if(!valid_p(it.index)) continue;
        TV X=it.Location();
        T p=ap.f(X,0),q=pressure(it.index),e=std::abs(p-q);
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));
        Debug_Particle_Set_Attribute<TV>("display_size",e);
        p_inf=std::max(p_inf,e);
        p_max=std::max(p_max,std::abs(p));
        p_l2+=e*e;
        p_cnt++;}
    p_l2=sqrt(p_l2/p_cnt);
    LOG::printf("p: %.3f %.3f   %.3f\n",p_inf,p_l2,p_max);
    Dump_Levelset(phi_grid,object_phi,VECTOR<T,3>(1,1,0));
    Dump_Levelset(phi_grid,surface_phi,VECTOR<T,3>(1,0,1));
    Flush_Frame("pressure error");

    return 0;
}

