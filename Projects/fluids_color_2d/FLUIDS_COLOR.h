//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR__
#define __FLUIDS_COLOR__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_BOX.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ELLIPSOID.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_LINE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_NEST.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ROTATE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SCALE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SPHERE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_TRANSLATE.h>
#include <PhysBAM_Geometry/Analytic_Tests/ANALYTIC_LEVELSET_VORTEX.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include "FLUIDS_COLOR_BASE.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace PhysBAM{

template<class TV>
struct ANALYTIC_VELOCITY_VORTEX:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    T nu,rho;
    ANALYTIC_VELOCITY_VORTEX(T mu,T rho): nu(mu/rho),rho(rho){}
    virtual TV u(const TV& X,T t) const
    {return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y))*exp(-2*nu*t);}
    virtual MATRIX<T,TV::m> du(const TV& X,T t) const
    {T c=cos(X.x)*cos(X.y),s=sin(X.x)*sin(X.y);return MATRIX<T,TV::m>(c,s,-s,-c)*exp(-2*nu*t);}
    virtual T p(const TV& X,T t) const {return (T).25*rho*(cos(2*X.x)+cos(2*X.y))*exp(-4*nu*t);}
    virtual TV F(const TV& X,T t) const {return TV();}
};
    
template<class TV>
struct ANALYTIC_VELOCITY_VORTEX_AND_SHIFT:public ANALYTIC_VELOCITY<TV>
{
    typedef typename TV::SCALAR T;
    T nu,rho;
    TV au;
    T const_p;
    ANALYTIC_VELOCITY_VORTEX_AND_SHIFT(T mu,T rho,TV v): nu(mu/rho),rho(rho),au(v),const_p(0){}
    virtual TV u(const TV& X,T t) const
    {return TV(sin(X.x)*cos(X.y),-cos(X.x)*sin(X.y))*exp(-2*nu*t)+au;}
    virtual MATRIX<T,2> du(const TV& X,T t) const
    {T c=cos(X.x)*cos(X.y),s=sin(X.x)*sin(X.y);return MATRIX<T,2>(c,s,-s,-c)*exp(-2*nu*t);}
    virtual T p(const TV& X,T t) const {return const_p+(T).25*rho*(cos(2*X.x)+cos(2*X.y))*exp(-4*nu*t);}
    virtual TV F(const TV& X,T t) const {return TV();}
};

template<class T>
class FLUIDS_COLOR<VECTOR<T,2> >:public FLUIDS_COLOR_BASE<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    typedef FLUIDS_COLOR_BASE<TV> BASE;

public:
    using BASE::grid;using BASE::use_level_set_method;using BASE::use_p_null_mode;using BASE::test_number;
    using BASE::analytic_velocity;using BASE::m;using BASE::s;using BASE::kg;using BASE::resolution;
    using BASE::analytic_levelset;using BASE::Large_Phi;using BASE::mu0;using BASE::mu1;using BASE::rho0;
    using BASE::rho1;using BASE::bc_type;using BASE::SLIP;using BASE::DIRICHLET;using BASE::NEUMANN;
    using BASE::unit_rho;using BASE::unit_mu;using BASE::unit_st;using BASE::surface_tension;
    using BASE::override_rho0;using BASE::override_rho1;using BASE::override_mu0;using BASE::override_mu1;
    using BASE::test_analytic_diff;using BASE::Initialize_Common_Example;using BASE::After_Initialize_Example;

    T epsilon,radius;
    int mode;

    FLUIDS_COLOR(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
        :FLUIDS_COLOR_BASE<TV>(stream_type,parse_args),epsilon((T).1),radius((T).05),mode(2)
    {
        parse_args.Add("-mode",&mode,"mode","Oscillation mode for surface tension test");
        parse_args.Add("-radius",&radius,"radius","Radius mode for surface tension test");
        parse_args.Add("-epsilon",&epsilon,"eps","Epsilon for surface tension test");
        parse_args.Parse();

        if(!Initialize_Common_Example())
            Initialize_Example();

        After_Initialize_Example();
    }

    ~FLUIDS_COLOR()
    {
    }

    void Initialize_Example()
    {
        switch(test_number){
            case 1:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX<TV>(mu0/unit_mu,rho0/unit_rho));
                use_p_null_mode=true;
                break;
            case 2:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX<TV>(mu0/unit_mu,rho0/unit_rho));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 5:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX<TV>(mu0/unit_mu,rho0/unit_rho));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 6:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+(T).5,VECTOR<T,1>(1),rho0/unit_rho));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 8:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE<TV>(new ANALYTIC_VELOCITY_VORTEX<TV>(mu0/unit_mu,rho0/unit_rho),TV((T).2,(T).5)));
                use_p_null_mode=true;
                break;
            case 10:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(TV()+1));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                break;
            case 16:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                {
                    TV vel((T).2,(T).5);
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4),vel);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE<TV>(new ANALYTIC_VELOCITY_VORTEX<TV>(mu0/unit_mu,rho0/unit_rho),vel));
                    if(bc_type!=NEUMANN) use_p_null_mode=true;
                }
                break;
            case 19:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                {
                    analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_SCALE<TV>(new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4),-2),VECTOR<T,1>(4)),TV(9,.4));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_TRANSLATE<TV>(new ANALYTIC_VELOCITY_VORTEX<TV>(mu0/unit_mu,rho0/unit_rho),TV(.5,-.2)));
                    if(bc_type!=NEUMANN) use_p_null_mode=true;
                }
                break;
            case 23:
                struct ANALYTIC_LEVELSET_MODE:public ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>
                {
                    T e,r;
                    int mode;
                    
                    ANALYTIC_LEVELSET_MODE(T e,T r,int mode,int c_i,int c_o): ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>(c_i,c_o),e(e),r(r),mode(mode) {}
                    virtual ~ANALYTIC_LEVELSET_MODE() {}
                    
                    virtual T f(const TV& X,T t) const {return X.Magnitude()-r*(1+e*cos(mode*atan2(X.y,X.x)));}
                    virtual TV df(const TV& X,T t) const {TV N=X;T d=N.Normalize();return N+r*mode*e*sin(mode*atan2(X.y,X.x))/d*N.Orthogonal_Vector();}
                    virtual MATRIX<T,TV::m> ddf(const TV& X,T t) const
                    {
                        TV N=X;
                        T d=N.Normalize(),m_th=mode*atan2(X.y,X.x),cs=cos(m_th),sn=sin(m_th);
                        TV O=N.Orthogonal_Vector();
                        MATRIX<T,TV::m> OO=MATRIX<T,TV::m>::Outer_Product(O,O),ON=MATRIX<T,TV::m>::Outer_Product(O*2,N).Symmetric_Part();
                        return (1/d+r*e*sqr(mode/d)*cs)*OO-r*e*mode*sn/sqr(d)*ON;
                    }
                    virtual TV Closest_Point_Estimate(const TV& X,T t) const {return X.Normalized()*r;}
                };

                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*(T).1*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_MODE(epsilon,radius,mode,0,1);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(TV()));
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST<TV>(TV()));
                surface_tension=(T)0.07197*unit_st;
                use_p_null_mode=true;
                use_level_set_method=true;
                if(!override_rho0) rho0=1000*unit_rho;
                if(!override_rho1) rho1=(T)1.1839*unit_rho;
                if(!override_mu0) mu0=0*unit_mu;
                if(!override_mu1) mu1=0*unit_mu;
//                analytic_initial_only=true;

                if(test_analytic_diff){
                    ANALYTIC_LEVELSET_MODE al(epsilon,radius,mode,0,1);
                    RANDOM_NUMBERS<T> rand;
                    T e=1e-6,t=rand.Get_Uniform_Number(0,1);
                    TV X=rand.Get_Uniform_Vector(grid.domain),dX;
                    rand.Fill_Uniform(dX,-e,e);
                    T f0=al.f(X,t),f1=al.f(X+dX,t);
                    TV df0=al.df(X,t),df1=al.df(X+dX,t);
                    MATRIX<T,TV::m> ddf0=al.ddf(X,t)/s,ddf1=al.ddf(X+dX,t)/s;
                    T err=abs((df0+df1).Dot(dX)/2-(f1-f0))/e;
                    T errd=((ddf0+ddf1)*dX/2-(df1-df0)).Magnitude()/e;
                    LOG::cout<<"analytic diff test f "<<err<<"  "<<errd<<std::endl;}
                break;
            case 27:{//Like test 25 but with an ellipse instead of the smaller circle
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                
                ANALYTIC_LEVELSET<TV>* ab=new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_ELLIPSOID<TV>(TV(.5,.3),TV(.15,.1),1,0),VECTOR<T,1>(1),TV()+(T).5);
                ANALYTIC_LEVELSET<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),-4,-4);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).42,0,1)))->Add(ab)->Add(cd);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+(T).5,VECTOR<T,1>(1),rho0/unit_rho));
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION<TV>(TV()+(T).5,VECTOR<T,1>(1),rho0/unit_rho));
                if(bc_type!=NEUMANN) use_p_null_mode=true;
                //use_level_set_method=true;
                
                break;}
            case 29:{
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*pi*m,true);
                TV vel((T).2,(T).5);
                analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3*pi,1,0),vel);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX_AND_SHIFT<TV>(mu0/unit_mu,rho0/unit_rho,vel));
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX_AND_SHIFT<TV>(mu0/unit_mu,rho0/unit_rho,vel));
                use_p_null_mode=true;
                use_level_set_method=true;
                break;}
            default: PHYSBAM_FATAL_ERROR("Missing test number");}
    }

//#####################################################################
};
}

#endif
