#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <cmath>
#include <cstdio>
#include "HEADER.h"
#include "OBJECTS_COMMON.h"
#include "PARAMETERS_COMMON.h"
#include "SIM_COMMON.h"
using namespace PhysBAM;

template<class T,class TV,int d>
void Invalidate_Outside(const GRID<TV>& grid,int ghost,ARRAY<T,FACE_INDEX<d> >& u,const BOUNDARY_CONDITIONS<TV>& bc,T value)
{
    for(FACE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next()) if(!bc.Inside(grid.Face(it.Full_Index()))) u(it.Full_Index())=value;
}

template<class TV>
void Fill_Ghost_Cells(const GRID<TV>& grid,int ghost,int distance,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,const BOUNDARY_CONDITIONS<TV>& bc)
{
    typedef typename TV::SCALAR T;
    ARRAY<bool,FACE_INDEX<TV::m> > inside(grid.Domain_Indices());
    for(FACE_ITERATOR<TV> it(grid);it.Valid();it.Next())
        inside(it.Full_Index())=bc.phi->Phi(it.Location())<=0;
    EXTRAPOLATION_HIGHER_ORDER<TV,T>(grid,*bc.phi,100,3,distance).Extrapolate_Face(inside,u);
}

template<class TV>
void Initialize_Grid_From_Domains(GRID<TV>& grid,int resolution,const RANGE<TV>& sample_box,const RANGE<TV>& bounding_box,RANGE<VECTOR<int,TV::m> >& sample_domain)
{
    TV dX=sample_box.Edge_Lengths()/resolution;
    TV o=sample_box.min_corner;
    VECTOR<int,TV::m> mn(floor((bounding_box.min_corner-o)/dX)),mx(ceil((bounding_box.max_corner-o)/dX));
    grid.Initialize(mx-mn,RANGE<TV>(TV(mn)*dX+o,TV(mx)*dX+o),true);
    sample_domain.min_corner-=mn;
    sample_domain.max_corner-=mn;
}

template<class TV>
void Prune_Outside_Sample_Points(const GRID<TV>& grid,const BOUNDARY_CONDITIONS<TV>& bc,ACCURACY_INFO<TV::m>& ai)
{
    int k=0;
    for(int i=1;i<ai.cell_samples.m;i++) if(bc.Inside(grid.X(ai.cell_samples(i)))) ai.cell_samples(++k)=ai.cell_samples(i);
    ai.cell_samples.Resize(k);
    k=0;
    for(int i=1;i<ai.face_samples.m;i++) if(bc.Inside(grid.Face(ai.face_samples(i)))) ai.face_samples(++k)=ai.face_samples(i);
    ai.face_samples.Resize(k);
}

template<class TV>
void Add_Viscosity(const OBJECTS_COMMON<TV>& obj,const PARAMETERS_COMMON<typename TV::SCALAR>& param,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,bool use_extrapolation)
{
    if(use_extrapolation){
        Fill_Ghost_Cells(obj.grid,3,3,u,*obj.bc);
        Flush_Frame(u,"fill 1");}
    for(int a=0;a<TV::m;a++) Apply_Viscosity(obj.grid,u,*obj.bc,param.dt,param.time,param.mu,param.rho,a,param.theta_threshold,param.cg_tolerance,param.print_matrix);
    Flush_Frame(u,"after viscosity");
    obj.ai.Print("AFTVIS",u);
}

template<class TV>
void Add_Advection(const OBJECTS_COMMON<TV>& obj,const PARAMETERS_COMMON<typename TV::SCALAR>& param,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,bool use_extrapolation)
{
    typedef typename TV::SCALAR T;
    ADVECTION_HAMILTON_JACOBI_ENO<TV,T> advection;
    BOUNDARY<TV,T> boundary;

    if(use_extrapolation){
        Fill_Ghost_Cells(obj.grid,3,3,u,*obj.bc);
        Flush_Frame(u,"after extrapolation");}
    advection.Update_Advection_Equation_Face(obj.grid,u,u,u,boundary,param.dt,param.time);
    Flush_Frame(u,"after advection");
    for(FACE_ITERATOR<TV> it(obj.grid,4);it.Valid();it.Next()) if(!obj.bc->Inside(it.Location())) u(it.Full_Index())=0;
    Flush_Frame(u,"after clear ghost");
}

template<class TV>
void Project_Pressure(const OBJECTS_COMMON<TV>& obj,const PARAMETERS_COMMON<typename TV::SCALAR>& param,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,bool use_extrapolation,proj_type proj)
{
    if(use_extrapolation){
        Fill_Ghost_Cells(obj.grid,3,3,u,*obj.bc);
        Flush_Frame(u,"fill after viscosity");}
    if(proj==proj_gibou) Project_Incompressibility_Gibou(obj.grid,u,*obj.bc,obj.ai,param.time,param.rho,param.theta_threshold,param.cg_tolerance,param.print_matrix);
    else if(proj==proj_slip)  Project_Incompressibility_Slip(obj.grid,u,*obj.bc,obj.ai,param.time,param.rho,param.theta_threshold,param.cg_tolerance,param.print_matrix);
    else if(proj==proj_default) Project_Incompressibility(obj.grid,u,*obj.bc,obj.ai,param.time,param.rho,param.theta_threshold,param.cg_tolerance,param.print_matrix);
    else PHYSBAM_FATAL_ERROR();
    obj.ai.Print("AFTPROJ",u);
    Flush_Frame(u,"after project");
}

template<class TV>
void First_Order_Step(SIM_COMMON<TV>& sim,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u)
{
    sim.obj.bc->Update_Parameters(sim.param);
    if(sim.use_viscosity) Add_Viscosity(sim.obj,sim.param,u,sim.use_extrapolation);
    if(sim.use_advection) Add_Advection(sim.obj,sim.param,u,sim.use_extrapolation);
    sim.param.time+=sim.param.dt;
    if(sim.use_projection) Project_Pressure(sim.obj,sim.param,u,sim.use_extrapolation,sim.proj_algo);
}

template<class TV>
void Second_Order_RE_Step(SIM_COMMON<TV>& sim,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u2)
{
    typedef typename TV::SCALAR T;
    u2=u;
    First_Order_Step(sim,u2);
    Flush_Frame(u,"lo step");
    sim.param.time-=sim.param.dt;

    sim.param.dt/=2;
    First_Order_Step(sim,u);
    First_Order_Step(sim,u);
    sim.param.dt*=2;
    Flush_Frame(u2,"hi step");

    u.Copy((T)2,u,-(T)1,u2);
}

template<class TV>
void Dump_Error(const SIM_COMMON<TV>& sim,const ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u2)
{
    for(FACE_ITERATOR<TV> it(sim.obj.grid,4);it.Valid();it.Next()){
        if(sim.obj.bc->Inside(it.Location())) u2(it.Full_Index())=u(it.Full_Index())-sim.obj.bc->Analytic_Velocity(it.Location(),sim.param.time)(it.Axis());
        else u2(it.Full_Index())=0;}
    Flush_Frame(u2,"error");
}

template<class T>
ERROR_COLOR_MAP<T>::ERROR_COLOR_MAP(T min_value,T max_value,bool log_scale,bool reverse,bool two_sets)
    :mn(min_value),mx(max_value),use_log(log_scale),reverse_order(reverse),use_two_sets(two_sets),colors(*new INTERPOLATION_CURVE<T,VECTOR<T,3> >)
{
    if(use_log)
    {
        mn=std::log(std::abs(mn));
        mx=std::log(std::abs(mx));
    }
    PHYSBAM_ASSERT(mn<mx);
    T a=mn,d=(mx-mn)/(6*(1+use_two_sets));
    if(reverse_order){a=mx,d=-d;}
    colors.Add_Control_Point(a+0*d,VECTOR<T,3>(1,1,1));
    colors.Add_Control_Point(a+1*d,VECTOR<T,3>(1,0,0));
    colors.Add_Control_Point(a+2*d,VECTOR<T,3>(1,1,0));
    colors.Add_Control_Point(a+3*d,VECTOR<T,3>(0,1,0));
    colors.Add_Control_Point(a+4*d,VECTOR<T,3>(0,1,1));
    colors.Add_Control_Point(a+5*d,VECTOR<T,3>(0,0,1));
    if(use_two_sets){
        colors.Add_Control_Point(a+6*d,VECTOR<T,3>(.3,0,.3));
        colors.Add_Control_Point(a+7*d,VECTOR<T,3>(.3,0,0));
        colors.Add_Control_Point(a+8*d,VECTOR<T,3>(.3,.3,0));
        colors.Add_Control_Point(a+9*d,VECTOR<T,3>(0,.3,0));
        colors.Add_Control_Point(a+10*d,VECTOR<T,3>(0,.3,.3));
        colors.Add_Control_Point(a+11*d,VECTOR<T,3>(0,0,.3));}
    colors.Add_Control_Point(a+6*(1+use_two_sets)*d,VECTOR<T,3>());
}

template<class T>
ERROR_COLOR_MAP<T>::~ERROR_COLOR_MAP()
{
    delete &colors;
}

template<class T>
VECTOR<T,3> ERROR_COLOR_MAP<T>::operator()(T x) const
{
    if(use_log){if(x==0) x=mn;else x=std::log(std::abs(x));}
    return colors.Value(x);
}

template<class TV>
void Dump_Error_Image(const SIM_COMMON<TV>& sim,const ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >& u,VECTOR<int,TV::m> dim)
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    static int id=0;id++;
    LOG::cout<<"ERROR IMAGE ID "<<id<<std::endl;

    GRID<TV> image_grid(dim,sim.obj.grid.domain);
    ARRAY<T,TV_INT> errors(image_grid.Domain_Indices());
    ARRAY<VECTOR<T,3>,TV_INT> color_errors(image_grid.Domain_Indices());
    ERROR_COLOR_MAP<T> em(1e-12,1,true,true,true);

    for(FACE_ITERATOR<TV> it(sim.obj.grid);it.Valid();it.Next()){
        if(!sim.obj.bc->Inside(it.Location())) continue;
        T err=abs(u(it.Full_Index())-sim.obj.bc->Analytic_Velocity(it.Location(),sim.param.time)(it.Axis()));
        TV_INT index;
        image_grid.Cell(it.Location(),index,0);
        if(err>errors(index)) errors(index)=err;}

    for(RANGE_ITERATOR<TV::m> it(errors.domain);it.Valid();it.Next())
        if(T e=errors(it.index)) color_errors(it.index)=em(e);

    OCTAVE_OUTPUT<T>(LOG::sprintf("error-image-%i.txt",id).c_str()).Write("EI",color_errors);
}

template struct ERROR_COLOR_MAP<double>;
template void Add_Advection<VECTOR<double,1> >(OBJECTS_COMMON<VECTOR<double,1> > const&,PARAMETERS_COMMON<double> const&,ARRAY<double,FACE_INDEX<1> >&,bool);
template void Add_Advection<VECTOR<double,2> >(OBJECTS_COMMON<VECTOR<double,2> > const&,PARAMETERS_COMMON<double> const&,ARRAY<double,FACE_INDEX<2> >&,bool);
template void Add_Viscosity<VECTOR<double,1> >(OBJECTS_COMMON<VECTOR<double,1> > const&,PARAMETERS_COMMON<double> const&,ARRAY<double,FACE_INDEX<1> >&,bool);
template void Add_Viscosity<VECTOR<double,2> >(OBJECTS_COMMON<VECTOR<double,2> > const&,PARAMETERS_COMMON<double> const&,ARRAY<double,FACE_INDEX<2> >&,bool);
template void Dump_Error<VECTOR<double,1> >(SIM_COMMON<VECTOR<double,1> > const&,ARRAY<double,FACE_INDEX<1> > const&,ARRAY<double,FACE_INDEX<1> >&);
template void Dump_Error<VECTOR<double,2> >(SIM_COMMON<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<2> > const&,ARRAY<double,FACE_INDEX<2> >&);
template void Dump_Error_Image<VECTOR<double,2> >(SIM_COMMON<VECTOR<double,2> > const&,ARRAY<VECTOR<double,2>::SCALAR,FACE_INDEX<VECTOR<double,2>::m> > const&,VECTOR<int,VECTOR<double,2>::m>);
template void Fill_Ghost_Cells<VECTOR<double,1> >(GRID<VECTOR<double,1> > const&,int,int,ARRAY<double,FACE_INDEX<1> >&,const BOUNDARY_CONDITIONS<VECTOR<double,1> >&);
template void Fill_Ghost_Cells<VECTOR<double,2> >(GRID<VECTOR<double,2> > const&,int,int,ARRAY<double,FACE_INDEX<2> >&,const BOUNDARY_CONDITIONS<VECTOR<double,2> >&);
template void First_Order_Step<VECTOR<double,2> >(SIM_COMMON<VECTOR<double,2> >&,ARRAY<VECTOR<double,2>::SCALAR,FACE_INDEX<VECTOR<double,2>::m> >&);
template void Initialize_Grid_From_Domains<VECTOR<double,1> >(GRID<VECTOR<double,1> >&,int,RANGE<VECTOR<double,1> > const&,RANGE<VECTOR<double,1> > const&,RANGE<VECTOR<int,VECTOR<double,1>::m> >&);
template void Initialize_Grid_From_Domains<VECTOR<double,2> >(GRID<VECTOR<double,2> >&,int,RANGE<VECTOR<double,2> > const&,RANGE<VECTOR<double,2> > const&,RANGE<VECTOR<int,VECTOR<double,2>::m> >&);
template void Project_Pressure<VECTOR<double,1> >(OBJECTS_COMMON<VECTOR<double,1> > const&,PARAMETERS_COMMON<double> const&,ARRAY<double,FACE_INDEX<1> >&,bool,proj_type);
template void Project_Pressure<VECTOR<double,2> >(OBJECTS_COMMON<VECTOR<double,2> > const&,PARAMETERS_COMMON<double> const&,ARRAY<double,FACE_INDEX<2> >&,bool,proj_type);
template void Prune_Outside_Sample_Points<VECTOR<double,1> >(GRID<VECTOR<double,1> > const&,BOUNDARY_CONDITIONS<VECTOR<double,1> > const&,ACCURACY_INFO<VECTOR<double,1>::m>&);
template void Prune_Outside_Sample_Points<VECTOR<double,2> >(GRID<VECTOR<double,2> > const&,BOUNDARY_CONDITIONS<VECTOR<double,2> > const&,ACCURACY_INFO<VECTOR<double,2>::m>&);
template void Second_Order_RE_Step<VECTOR<double,1> >(SIM_COMMON<VECTOR<double,1> >&,ARRAY<double,FACE_INDEX<1> >&,ARRAY<double,FACE_INDEX<1> >&);
template void Second_Order_RE_Step<VECTOR<double,2> >(SIM_COMMON<VECTOR<double,2> >&,ARRAY<double,FACE_INDEX<2> >&,ARRAY<double,FACE_INDEX<2> >&);
