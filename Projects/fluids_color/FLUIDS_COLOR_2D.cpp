//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ELLIPSOID.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_LINE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_MOVING_ELLIPSOID.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_NEST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ROTATE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SCALE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SPHERE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_TRANSLATE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_VORTEX.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include "FLUIDS_COLOR_2D.h"
namespace PhysBAM{
//#####################################################################
// Function Force_Value
//#####################################################################
template<class T,class A> enable_if_t<!is_same<T,A>::value,T> Force_Value(A v){return v.x;}
template<class T> T Force_Value(T v){return v;}
//#####################################################################
// Constructor
//#####################################################################
template<class T> FLUIDS_COLOR<VECTOR<T,2> >::
FLUIDS_COLOR(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :FLUIDS_COLOR_BASE<TV>(stream_type_input,parse_args),epsilon((T).1),radius((T).05),mode(2)
{
    parse_args.Add("-mode",&mode,"mode","Oscillation mode for surface tension test");
    parse_args.Add("-radius",&radius,"radius","Radius mode for surface tension test");
    parse_args.Add("-epsilon",&epsilon,"eps","Epsilon for surface tension test");
    parse_args.Parse();

    if(!Initialize_Common_Example())
        Initialize_Example();

    After_Initialize_Example();
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FLUIDS_COLOR<VECTOR<T,2> >::
~FLUIDS_COLOR()
{
}
////////////////////////////////////////////////////////////////////////////
//
//   Tests:
//   1. Periodic decaying vortex, no interface
//   2. Translating and decaying vortex, interface
//   3. In BASE class
//   4. In BASE class
//   5. Vortex with sphere boundary
//   6. Rotation with vortex stream boundary
//   7. In BASE class
//   8. Translating and decaying vortex
//   9. In BASE class
//   10. Constant velocity with vortex stream boundary
//   11-15. In BASE class
//   16. Translating and decaying vortex, boundary
//   17-18. In BASE class
//   19. Distorted vortex with boundary
//   20-22. In BASE class
//   23.
//   24-26. In BASE class
//   27. 
//   28. In BASE class
//   29.
//   30.
//
//   101. Decaying vortex, identical to Test 2. Example 1 in the JCP paper.
//   102. Translating and decaying vortex, similar to Test 8. JCP Example 2.
//   103. Different velocity fields with embedded Neumann bc. JCP Example 3.
//   104. In BASE class. JCP Example 4 (2d).
//   105. In BASE class. JCP Example 5 (3d).
//   106. In BASE class. JCP Example 6, also Test 9.
//   107. In BASE class. JCP Example 7, also Test 22.
//
//   
////////////////////////////////////////////////////////////////////////////
template<class T> void FLUIDS_COLOR<VECTOR<T,2> >::
Initialize_Example()
{
    SYMMETRIC_MATRIX<T,TV::m> id(SYMMETRIC_MATRIX<T,TV::m>()+1);
    switch(test_number){
        case 1:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            Add_Vortex(mu0/unit_mu,rho0/unit_rho);
            use_p_null_mode=true;
            break;
        case 2:
        case 101:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4);
            Add_Vortex(mu0/unit_mu,rho0/unit_rho);
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            break;
        case 5:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3,0,-4);
            Add_Vortex(mu0/unit_mu,rho0/unit_rho);
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            break;
        case 6:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4);
            auto rot=[=](auto X,auto t){return MATRIX<T,2>::Cross_Product_Matrix(VECTOR<T,1>(1))*(X-(TV()+(T).5));};
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(this->rho0/this->unit_rho)*rot(X,t).Magnitude_Squared();});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            break;}
        case 8:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            Add_Vortex(mu0/unit_mu,rho0/unit_rho,TV((T).2,(T).5));
            use_p_null_mode=true;
            break;
        case 10:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4);
            Add_Velocity([](auto X,auto t){return TV()+1;});
            Add_Pressure([](auto X,auto t){return 0;});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            break;
        case 16:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
            {
                TV vel((T).2,(T).5);
                analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4),vel);
                Add_Vortex(mu0/unit_mu,rho0/unit_rho,vel);
                if(bc_type!=NEUMANN) use_p_null_mode=true;
            }
            break;
        case 19:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
            {
                analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_ROTATE<TV>(new ANALYTIC_LEVELSET_SCALE<TV>(new ANALYTIC_LEVELSET_VORTEX<TV>((T).2,0,-4),-2),VECTOR<T,1>(4)),TV(9,.4));
                Add_Vortex(mu0/unit_mu,rho0/unit_rho,TV(.5,-.2));
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
            Add_Velocity([](auto X,auto t){return TV();});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Velocity([](auto X,auto t){return TV();});
            Add_Pressure([](auto X,auto t){return 0;});
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
            auto rot=[=](auto X,auto t){return MATRIX<T,2>::Cross_Product_Matrix(VECTOR<T,1>(1))*(X-(TV()+(T).5));};
            Add_Velocity(rot);
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(this->rho0/this->unit_rho)*rot(X,t).Magnitude_Squared();});
            Add_Pressure([=](auto X,auto t){return (T).5*(this->rho1/this->unit_rho)*rot(X,t).Magnitude_Squared();});
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            //use_level_set_method=true;
                
            break;}

        case 29:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*pi*m,true);
            TV vel((T).0,(T).0);
            T lambda = (T)0;
            analytic_levelset=new ANALYTIC_LEVELSET_TRANSLATE<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T).5,(T).3*pi,1,0),vel);
            Add_Vortex(mu0/unit_mu,rho0/unit_rho,vel,lambda);
            Add_Vortex(mu1/unit_mu,rho1/unit_rho,vel,lambda);
            use_p_null_mode=true;
            use_level_set_method=true;
            break;}
        case 30:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*pi*m,true);
            TV vel((T).5,(T).3);
            T lambda = (T)1;
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,-4);
            Add_Vortex(mu0/unit_mu,rho0/unit_rho,vel,lambda);
            use_p_null_mode=true;
            use_level_set_method=true;
            break;}
        case 32:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            {
                T radius=(T).6,curvature=(TV::m-1)/radius;
                if(!override_surface_tension) surface_tension=(T)1*unit_st;
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),radius,0,1);
                Add_Velocity([=](auto X,auto t){auto q=X.Magnitude_Squared()-sqr(radius);return sqr(q)*(T)12*(MATRIX<T,2>(0,-1,1,0)*X);});
                Add_Pressure([=](auto X,auto t){return (this->surface_tension*curvature)/this->unit_p;});
                use_p_null_mode=true;
                use_level_set_method=true;
                use_advection=false;
            }
            break;
        case 33:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            {
                T outer_radius=(T).9,ellipse_mean_radius=(T).5;
                    
                auto a=[=](auto t)
                    {
                        return (T)1*(mode==0)+(T)1.1*(mode==1)+(t+1)*(mode==2)+exp(t)*(mode==3)+((T)1+(T).2*exp(-t))*(mode==4);
                    };
                auto da=[=](auto t)
                    {
                        return t*(mode==2)+exp(t)*(mode==3)-(T).2*exp(-t)*(mode==4);
                    };

                MATRIX<T,TV::m> sm(1,0,0,-1);
                if(!override_surface_tension) surface_tension=(T)1*unit_st;
                ANALYTIC_LEVELSET<TV>* ab=new ANALYTIC_LEVELSET_MOVING_ELLIPSOID<TV>(TV(),[=](T t){T r=a(t);return ellipse_mean_radius*TV(r,1/r);},0,1);
                ANALYTIC_LEVELSET<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),-4,-4);
                analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),outer_radius,0,1)))->Add(ab)->Add(cd);

                T st=surface_tension/unit_st;
                T mu_j=(mu1-mu0)/unit_mu;

                Add_Velocity([=](auto X,auto t){return da(t)/a(t)*(sm*X);});
                Add_Velocity([=](auto X,auto t){return da(t)/a(t)*(sm*X);});
                Add_Pressure([=](auto X,auto t){return 0;});
                Add_Pressure([=](auto X,auto t)
                    {
                        auto x=X.Dot(TV(1,0));
                        auto y=X.Dot(TV(0,1));
                        auto NN=(TV(1,0)*x+TV(0,1)*sqr(sqr(a(t)))*y);
                        auto N=NN/NN.Magnitude();
                        auto a4=sqr(sqr(a(t)));
                        auto x2=sqr(x);
                        auto z=a4*sqr(y);
                        auto e=a4/cube(sqrt(a4*z+x2));
                        auto K=e*(x2+z);
                        return -st*K-2*mu_j*N.Dot(da(t)/a(t)*sm*N);
                    });

                if(bc_type!=NEUMANN) use_p_null_mode=true;
                use_level_set_method=true;
            }
            break;
        case 102:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV((T)1.1*(T)pi,(T)pi),(T).6*(T)pi,0,1);
            TV vel((T).2,(T).5);
            Add_Vortex(mu0/unit_mu,rho0/unit_rho,vel);
            Add_Vortex(mu1/unit_mu,rho1/unit_rho,vel);
            use_level_set_method=true;
            use_p_null_mode=true;
            break;}
        case 103:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*(T)pi*m,true);
            ANALYTIC_LEVELSET<TV>* ab=new ANALYTIC_LEVELSET_SPHERE<TV>(TV::Axis_Vector(0)*.2*pi,(T).2*(T)pi,NEUMANN,0);
            ANALYTIC_LEVELSET<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),1,1);
            analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),(T).8*(T)pi,0,1)))->Add(ab)->Add(cd);
            MATRIX<T,TV::m> du0;for(int i=0;i<TV::m;i++)du0(i,i)=-1;du0(1,1)+=TV::m;
            Add_Velocity([=](auto X,auto t){return du0*(X-TV::Axis_Vector(0)*.2*pi);});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Vortex(mu1/unit_mu,rho1/unit_rho);
            use_discontinuous_velocity=true;
            break;}
        case 111:{
            grid.Initialize(TV_INT(resolution,2*resolution),RANGE<TV>(TV(-1,0),TV(1,4))*m,true);
            T top = (T)3.5,bottom=(T).5;
            TV bubble_center((T)0,(T)1.5),bubble_radius((T).5,(T).2);
            ANALYTIC_LEVELSET<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),DIRICHLET,DIRICHLET);

            ANALYTIC_LEVELSET_SIGNED<TV>* ef=new ANALYTIC_LEVELSET_ELLIPSOID<TV>(bubble_center,bubble_radius,1,0);
            ANALYTIC_LEVELSET<TV>* gh=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(1)*bottom,TV::Axis_Vector(1),0,1)))->Add(cd)->Add(ef);
            analytic_levelset = (new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(1)*top,TV::Axis_Vector(1),0,1)))->Add(gh)->Add(cd);
            gravity=TV::Axis_Vector(0)*(-(T)9.8)*m/s/s;
            Add_Velocity([](auto X,auto t){return TV();});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Velocity([](auto X,auto t){return TV();});
            Add_Pressure([](auto X,auto t){return 0;});
            analytic_initial_only=true;
            use_level_set_method=true;
            use_p_null_mode=true;

            break;}
        case 112:{
            grid.Initialize(TV_INT(2*resolution,resolution),RANGE<TV>(TV(0,0),TV(12,6))*m,true);
            T left = (T)1,right=(T)11;
            TV cylinder_center((T)3.25,(T)3);//,bubble_radius((T).5,(T).2);
            T cylinder_radius(.5);
            ANALYTIC_LEVELSET<TV>* ab=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),DIRICHLET,DIRICHLET);
            ANALYTIC_LEVELSET<TV>* ef=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),NEUMANN,NEUMANN);                
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_SPHERE<TV>(cylinder_center,cylinder_radius,SLIP,0);
            ANALYTIC_LEVELSET<TV>* gh=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*left,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
            analytic_levelset = (new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*right,TV::Axis_Vector(0),0,1)))->Add(gh)->Add(ef);
            ANALYTIC_LEVELSET<TV>* mn=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*1.3,TV::Axis_Vector(0),0,1);                
            Add_Velocity([=](auto X,auto t){int c=0;mn->phi(Force_Value<TV>(X),Force_Value<T>(t),c);return TV(1-c,0);});
            Add_Pressure([](auto X,auto t){return 0;});
            use_level_set_method=true;
            use_p_null_mode=false;                
            break;}
        case 253:{
            TV a;a(0)=rho0/unit_rho;
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),(T).6,0,-4);
            Add_Velocity([](auto X,auto t){return TV()+1;});
            Add_Pressure([](auto X,auto t){return 0;});
            Add_Polymer_Stress([=](auto X,auto t){return id*(X.Dot(a)+(T)1);});
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }
        case 254:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),(T).6,0,-4);
            auto rot=[=](auto X,auto t){return MATRIX<T,2>::Cross_Product_Matrix(VECTOR<T,1>(1))*(X-(TV()+(T).5));};
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(this->rho0/this->unit_rho)*rot(X,t).Magnitude_Squared();});
            T k=rho0/unit_rho;
            Add_Polymer_Stress([=](auto X,auto t){return id*(k*X.Magnitude_Squared()+(T)1);});
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }
        case 256:{
            SYMMETRIC_MATRIX<T,TV::m> sms(rho0/unit_rho,0,-rho0/unit_rho),smc(0,rho0/unit_rho,0);
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),(T).6,0,-4);
            auto rot=[=](auto X,auto t){return MATRIX<T,2>::Cross_Product_Matrix(VECTOR<T,1>(1))*(X-(TV()+(T).5));};
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(this->rho0/this->unit_rho)*rot(X,t).Magnitude_Squared();});
            Add_Polymer_Stress([=](auto X,auto t){auto q=X.Magnitude_Squared()-t*2;return id+sms*sin(q)+smc*cos(q);});
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }
        case 257:{
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Centered_Box()*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV(),(T).6,0,-4);
            auto rot=[=](auto X,auto t){return MATRIX<T,2>::Cross_Product_Matrix(VECTOR<T,1>(1))*(X-(TV()+(T).5));};
            Add_Velocity(rot);
            Add_Pressure([=](auto X,auto t){return (T).5*(this->rho0/this->unit_rho)*rot(X,t).Magnitude_Squared();});
            analytic_polymer_stress.Append(new ANALYTIC_SYMMETRIC_MATRIX_QUADRATIC<TV>(rho0/unit_rho));
            if(!override_beta0) polymer_stress_coefficient.Append(1.2*unit_p);
            if(!override_inv_Wi0) inv_Wi.Append(1.3/s);
            if(bc_type!=NEUMANN) use_p_null_mode=true;
            use_polymer_stress=true;
            break;
        }

        default: PHYSBAM_FATAL_ERROR("Missing test number");}
}
//#####################################################################
// Function Add_Vortex
//#####################################################################
template<class T> void FLUIDS_COLOR<VECTOR<T,2> >::
Add_Vortex(T mu,T rho,TV trans_vel,T new_l)
{
    T nu=mu/rho;
    Add_Velocity(
        [=](auto X,auto t)
        {
            auto Z=X-trans_vel*t;
            auto x=Z.Dot(TV(1,0));
            auto y=Z.Dot(TV(0,1));
            return (TV(1,0)*sin(x)*cos(y)-TV(0,1)*cos(x)*sin(y))*exp(-2*nu*t)*cos(new_l*t)+trans_vel;
        });
    Add_Pressure(
        [=](auto X,auto t)
        {
            auto Z=X-trans_vel*t;
            auto x=Z.Dot(TV(1,0));
            auto y=Z.Dot(TV(0,1));
            return (T).25*rho*(cos(x*2)+cos(y*2))*exp(-4*nu*t)*sqr(cos(new_l*t));
        });
}
//#####################################################################
template class FLUIDS_COLOR<VECTOR<float,2> >;
template class FLUIDS_COLOR<VECTOR<double,2> >;
}
