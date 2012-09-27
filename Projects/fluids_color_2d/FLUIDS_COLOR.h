//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR__
#define __FLUIDS_COLOR__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>

namespace PhysBAM{


template<class TV>
class ANALYTIC_IMPLICIT_SURFACE_LEVELSET
{
    typedef typename TV::SCALAR T;
public:

    T tolerance;
    ANALYTIC_IMPLICIT_SURFACE_LEVELSET(): tolerance((T)1e-25) {}

    virtual T f(const TV& X) const=0;
    virtual TV df(const TV& X) const=0;
    virtual MATRIX<T,TV::m> ddf(const TV& X) const=0;
    virtual TV Closest_Point_Estimate(const TV& X) const {return X;}

    VECTOR<T,TV::m+1> Find_Closest_Point(const TV& X) const
    {
        TV w=Closest_Point_Estimate(X);
        VECTOR<T,TV::m+1> z(w.Append(0));
        for(int i=0;i<100;i++){
            TV Z(z.Remove_Index(TV::m)),dg=df(Z);
            T L=z(TV::m),g=f(Z);
            VECTOR<T,TV::m+1> G=((Z-X)*2+L*dg).Append(g),Hcol=dg.Append(0);
            MATRIX<T,TV::m+1> H;
            H.Set_Submatrix(0,0,L*ddf(Z)+2);
            H.Set_Row(TV::m,Hcol);
            H.Set_Column(TV::m,Hcol);
            z-=H.Solve_Linear_System(G);
            if(G.Magnitude_Squared()<1e-25) break;}
        return z;
    }

    T Phi(const TV& X) const
    {return (X-Find_Closest_Point(X).Remove_Index(TV::m)).Magnitude()*sign(f(X));}

    TV Normal(const TV& X) const
    {return (X-Find_Closest_Point(X).Remove_Index(TV::m)).Normalized();}
};

template<class TV>
class FLUIDS_COLOR:public PLS_FC_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef PLS_FC_EXAMPLE<TV> BASE;

public:
    using BASE::grid;using BASE::output_directory;using BASE::domain_boundary;using BASE::face_velocities;
    using BASE::particle_levelset_evolution;using BASE::write_substeps_level;using BASE::restart;using BASE::last_frame;
    using BASE::dt;using BASE::levelset_color;using BASE::mu;using BASE::rho;using BASE::dump_matrix;using BASE::number_of_colors;
    using BASE::use_advection;using BASE::use_reduced_advection;using BASE::omit_solve;

    enum WORKAROUND{SLIP=-3,DIRICHLET=-2,NEUMANN=-1}; // From CELL_DOMAIN_INTERFACE_COLOR

    int test_number;
    int resolution;
    int stored_last_frame;
    bool user_last_frame;
    T mu0,mu1;
    T rho0,rho1;
    T m,s,kg;
    int bc_type;
    bool bc_n,bc_d,bc_s;
    bool test_analytic_diff;
    bool no_advection;
    int refine;

    struct ANALYTIC_VELOCITY
    {
        virtual TV u(const TV& X,T t) const=0;
        virtual MATRIX<T,2> du(const TV& X,T t) const=0;
        virtual T p(const TV& X,T t) const=0;
    };
    ARRAY<ANALYTIC_VELOCITY*> analytic_velocity;

    struct ANALYTIC_LEVELSET
    {
        virtual T phi(const TV& X,T t,int& c) const=0;
        virtual TV N(const TV& X,T t) const=0;
    }* analytic_levelset;

    FLUIDS_COLOR(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
        :PLS_FC_EXAMPLE<TV>(stream_type),test_number(0),resolution(32),stored_last_frame(0),user_last_frame(false),mu0(1),mu1(2),
        rho0(1),rho1(2),m(1),s(1),kg(1),bc_n(false),bc_d(false),bc_s(false),no_advection(false),refine(1)
    {
        last_frame=16;
        parse_args.Set_Extra_Arguments(1,"<test-number>");
        parse_args.Add("-restart",&restart,"frame","restart frame");
        parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
        parse_args.Add("-substep",&write_substeps_level,"level","output-substep level");
        parse_args.Add("-dt",&dt,"dt","time step size to use for simulation");
        parse_args.Add("-steps",&this->time_steps_per_frame,"steps","number of time steps per frame");
        parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
        parse_args.Add("-dump_matrix",&dump_matrix,"dump out system and rhs");
        parse_args.Add("-mu0",&mu0,"viscosity","viscosity for first fluid region");
        parse_args.Add("-mu1",&mu1,"viscosity","viscosity for second fluid region");
        parse_args.Add("-rho0",&rho0,"density","density for first fluid region");
        parse_args.Add("-rho1",&rho1,"density","density for second fluid region");
        parse_args.Add("-m",&m,"scale","meter scale");
        parse_args.Add("-s",&s,"scale","second scale");
        parse_args.Add("-kg",&kg,"scale","kilogram scale");
        parse_args.Add("-bc_n",&bc_n,"use Neumann boundary conditions");
        parse_args.Add("-bc_d",&bc_d,"use Dirichlet boundary conditions");
        parse_args.Add("-bc_s",&bc_s,"use slip boundary conditions");
        parse_args.Add("-test_diff",&test_analytic_diff,"test analytic derivatives");
        parse_args.Add("-no_advect",&no_advection,"Disable advection");
        parse_args.Add("-no_solve",&omit_solve,"Disable visocity and pressure solve");
        parse_args.Add("-reduced_advect",&use_reduced_advection,"Peform reduced advection");
        parse_args.Add("-refine",&refine,"num","Refine space/time by this factor");
        parse_args.Parse();
        if(!STRING_UTILITIES::String_To_Value(parse_args.Extra_Arg(0),test_number)) throw VALUE_ERROR("The argument is not an integer.");

        resolution*=refine;
        dt/=refine;
        this->time_steps_per_frame*=refine;
        stored_last_frame=last_frame;
        mu0*=kg/s;
        mu1*=kg/s;
        rho0*=kg/sqr(m);
        rho1*=kg/sqr(m);
        dt*=s;
        PHYSBAM_ASSERT(bc_n+bc_d+bc_s<2);
        bc_type=bc_n?NEUMANN:(bc_s?SLIP:DIRICHLET);
        use_advection=!no_advection;

        analytic_levelset=0;

        output_directory=STRING_UTILITIES::string_sprintf("Test_%d",test_number);
        switch(test_number){
            case 0:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_PERIODIC;
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(TV()+1));
                break;
            case 1:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_PERIODIC;
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg,TV()));
                break;
            case 2:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_VORTEX((T).2);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg,TV()));
                break;
            case 3:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CIRCLE(TV()+(T).5,(T).3);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION(TV()+(T).5,rho0));
                break;
            case 4:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CIRCLE(TV()+(T).5,(T).3);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(TV()+1));
                break;
            case 5:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CIRCLE(TV()+(T).5,(T).3);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg,TV()));
                break;
            case 6:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_VORTEX((T).2);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION(TV()+(T).5,rho0));
                break;
            case 7:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_PERIODIC;
                analytic_velocity.Append(new ANALYTIC_VELOCITY_RAREFACTION);
                omit_solve=true;
                break;
            case 8:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_PERIODIC;
                analytic_velocity.Append(new ANALYTIC_VELOCITY_VORTEX(mu0*s/kg,rho0*sqr(m)/kg,TV((T).2,(T).5)));
                break;
            case 9:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                {
                    T x0=(T).2,x1=(T).5,x2=(T).8,v0=1,v2=-1;
                    T v1=(v0*mu0/(x1-x0)+v2*mu1/(x2-x1))/(mu0/(x1-x0)+mu1/(x2-x1));
                    MATRIX<T,2> du0(0,(v1-v0)/(x1-x0),0,0),du1(0,(v2-v1)/(x2-x1),0,0);
                    analytic_levelset=(new ANALYTIC_LEVELSET_BANDED(x0,x2,DIRICHLET,DIRICHLET))->Add(x1);
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE(TV(x1,0),TV(0,v1),du0));
                    analytic_velocity.Append(new ANALYTIC_VELOCITY_AFFINE(TV(x1,0),TV(0,v1),du1));
                }
                break;
            case 10:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(T)pi*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_VORTEX((T).2);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(TV()+1));
                break;
            case 11:
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_CIRCLE(TV()+(T).5,(T).3);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_ROTATION(TV((T).6,(T).8),rho0));
                break;
            default: PHYSBAM_FATAL_ERROR("Missing test number");}

        if(analytic_velocity.m) number_of_colors=analytic_velocity.m;
    }

    void Write_Output_Files(const int frame)
    {BASE::Write_Output_Files(frame);}

    void Initialize()
    {
        if(analytic_levelset && analytic_velocity.m) Analytic_Test();
        else
            switch(test_number){
                default: PHYSBAM_FATAL_ERROR("Missing test number");}

        PHYSBAM_ASSERT(rho.m==number_of_colors && mu.m==number_of_colors);
        if(user_last_frame) last_frame=stored_last_frame;
    }

    void Analytic_Test()
    {
        mu.Append(mu0);
        rho.Append(rho0);
        for(int i=1;i<number_of_colors;i++){
            mu.Append(mu1);
            rho.Append(rho1);}
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
            int c=-4;
            analytic_levelset->phi(it.Location()/m,0,c);
            if(c<0) continue;
            face_velocities(c)(it.Full_Index())=analytic_velocity(c)->u(it.Location()/m,0)(it.Axis())*m/s;}
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=-4;
            T p=analytic_levelset->phi(it.Location()/m,0,c)*m;
            levelset_color.phi(it.index)=abs(p);
            levelset_color.color(it.index)=c==-4?bc_type:c;}

        if(test_analytic_diff){
            RANDOM_NUMBERS<T> rand;
            TV X,dX;
            T e=1e-6,t=rand.Get_Uniform_Number(0,1);
            rand.Fill_Uniform(dX,-e,e);
            int c=-4;
            do{
                X=rand.Get_Uniform_Vector(grid.domain);
                analytic_levelset->phi(X/m,0,c);}
            while(c<0);
            TV u0=analytic_velocity(c)->u(X/m,t)*m/s,u1=analytic_velocity(c)->u((X+dX)/m,t)*m/s;
            MATRIX<T,2> du0=analytic_velocity(c)->du(X/m,t)/s,du1=analytic_velocity(c)->du((X+dX)/m,t)/s;
            T erru=((du0+du1)*dX/2-(u1-u0)).Magnitude()/e;
            int c0,c1;
            T l0=analytic_levelset->phi(X/m,t,c0)*m,l1=analytic_levelset->phi((X+dX)/m,t,c1)*m;
            if(c0>=0) l0=-l0;
            if(c1>=0) l1=-l1;
            TV dl0=analytic_levelset->N(X/m,t),dl1=analytic_levelset->N((X+dX)/m,t);
            T errl=abs((dl0+dl1).Dot(dX)/2-(l1-l0))/e;
            LOG::cout<<"analytic diff test "<<erru<<"  "<<errl<<std::endl;}
    }

    struct ANALYTIC_VELOCITY_CONST:public ANALYTIC_VELOCITY
    {
        TV au;
        ANALYTIC_VELOCITY_CONST(TV v): au(v){}
        virtual TV u(const TV& X,T t) const {return au;}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>();}
        virtual T p(const TV& X,T t) const {return 0;}
    };

    struct ANALYTIC_VELOCITY_ROTATION:public ANALYTIC_VELOCITY
    {
        TV c;
        T rho;
        ANALYTIC_VELOCITY_ROTATION(TV cc,T rho): c(cc),rho(rho){}
        virtual TV u(const TV& X,T t) const {return (X-c).Orthogonal_Vector();}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>(0,1,-1,0);}
        virtual T p(const TV& X,T t) const {return (T).5*rho*(X-c).Magnitude_Squared();}
    };

    struct ANALYTIC_VELOCITY_VORTEX:public ANALYTIC_VELOCITY
    {
        T nu,rho;
        TV trans;
        ANALYTIC_VELOCITY_VORTEX(T mu,T rho,TV t): nu(mu/rho),rho(rho),trans(t){}
        virtual TV u(const TV& X,T t) const
        {TV Z=X-t*trans;return TV(sin(Z.x)*cos(Z.y),-cos(Z.x)*sin(Z.y))*exp(-2*nu*t)+trans;}
        virtual MATRIX<T,2> du(const TV& X,T t) const
        {TV Z=X-t*trans;T c=cos(Z.x)*cos(Z.y),s=sin(Z.x)*sin(Z.y);return MATRIX<T,2>(c,s,-s,-c)*exp(-2*nu*t);}
        virtual T p(const TV& X,T t) const {TV Z=X-t*trans;return (T).25*rho*(cos(2*Z.x)+cos(2*Z.y))*exp(-4*nu*t);}
    };

    struct ANALYTIC_VELOCITY_RAREFACTION:public ANALYTIC_VELOCITY
    {
        ANALYTIC_VELOCITY_RAREFACTION(){}
        virtual TV u(const TV& X,T t) const {return TV(X.x,0)/(t+1);}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>(1,0,0,0)/(t+1);}
        virtual T p(const TV& X,T t) const {return 0;}
    };

    struct ANALYTIC_VELOCITY_ROTATION_DECAY:public ANALYTIC_VELOCITY
    {
        ANALYTIC_VELOCITY_ROTATION_DECAY(){}
        virtual TV u(const TV& X,T t) const {return X.Orthogonal_Vector()*2*exp(-X.Magnitude_Squared());}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return MATRIX<T,2>(-2*X.x*X.y,1-2*sqr(X.y),-1+2*sqr(X.x),2*X.x*X.y)*2*exp(-X.Magnitude_Squared());}
        virtual T p(const TV& X,T t) const {return 0;}
    };

    struct ANALYTIC_VELOCITY_AFFINE:public ANALYTIC_VELOCITY
    {
        TV v0;
        MATRIX<T,2> du0;
        ANALYTIC_VELOCITY_AFFINE(const TV& x0,const TV& v0,const MATRIX<T,2>& du0): v0(v0-du0*x0),du0(du0) {}
        virtual TV u(const TV& X,T t) const {return du0*X+v0;}
        virtual MATRIX<T,2> du(const TV& X,T t) const {return du0;}
        virtual T p(const TV& X,T t) const {return 0;}
    };

    struct ANALYTIC_LEVELSET_PERIODIC:public ANALYTIC_LEVELSET
    {
        virtual T phi(const TV& X,T t,int& c) const {c=0;return 1;}
        virtual TV N(const TV& X,T t) const {return TV(1,0);}
    };

    struct ANALYTIC_LEVELSET_CIRCLE:public ANALYTIC_LEVELSET
    {
        TV cen;
        T r;
        ANALYTIC_LEVELSET_CIRCLE(TV cc,T rr): cen(cc),r(rr){}
        virtual T phi(const TV& X,T t,int& c) const {T p=(X-cen).Magnitude()-r;c=p>=0?-4:0;return abs(p);}
        virtual TV N(const TV& X,T t) const {return (X-cen).Normalized();}
    };

    struct ANALYTIC_LEVELSET_VORTEX:public ANALYTIC_LEVELSET
    {
        T k;

        struct VORTEX_IMPLICIT_SURFACE:public ANALYTIC_IMPLICIT_SURFACE_LEVELSET<TV>
        {
            T k;
            virtual T f(const TV& X) const {return k-sin(X.x)*sin(X.y);}
            virtual TV df(const TV& X) const {return -TV(cos(X.x)*sin(X.y),sin(X.x)*cos(X.y));}
            virtual MATRIX<T,TV::m> ddf(const TV& X) const {T A=sin(X.x)*sin(X.y),B=-cos(X.x)*cos(X.y);return MATRIX<T,TV::m>(A,B,B,A);}
            virtual TV Closest_Point_Estimate(const TV& X) const {return (X-pi/2).Normalized()+pi/2;}
        } vis;

        ANALYTIC_LEVELSET_VORTEX(T kk): k(kk){vis.k=k;}

        virtual T phi(const TV& X,T t,int& c) const
        {T p=vis.Phi(X);c=p>0?-4:0;return abs(p);}

        virtual TV N(const TV& X,T t) const
        {return vis.Normal(X)*sign(vis.f(X));}
    };

    struct ANALYTIC_LEVELSET_BANDED:public ANALYTIC_LEVELSET
    {
        int bc0,bc1;
        T x0,x1;
        ARRAY<T> x;
        ANALYTIC_LEVELSET_BANDED(T x0,T x1,int bc0,int bc1): bc0(bc0),bc1(bc1),x0(x0),x1(x1) {}
        ANALYTIC_LEVELSET_BANDED* Add(T nx) {x.Append(nx);return this;}
        virtual T phi(const TV& X,T t,int& c) const
        {
            T p0=X.x-x0;
            if(p0<=0){c=bc0;return -p0;}
            for(int i=0;i<x.m;i++){
                T p1=X.x-x(i);
                if(p1<=0){c=i;return min(p0,-p1);}
                p0=p1;}
            T p1=X.x-x1;
            if(p1<=0){c=x.m;return min(p0,-p1);}
            c=DIRICHLET;
            return bc1;
        }
        virtual TV N(const TV& X,T t) const
        {
            T p0=X.x-x0;
            if(p0<=0) return TV(-1,0);
            for(int i=0;i<x.m;i++){
                T p1=X.x-x(i);
                if(p1<=0) return -TV(p0<-p1?1:-1,0);
                p0=p1;}
            T p1=X.x-x1;
            if(p1<=0) return -TV(p0<-p1?1:-1,0);
            return TV(1,0);
        }
    };

    void Begin_Time_Step(const T time) PHYSBAM_OVERRIDE {}

    void End_Time_Step(const T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m){
            T max_error=0,a=0,b=0;
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
                int c=levelset_color.Color(it.Location());
                if(c<0) continue;
                T A=face_velocities(c)(it.Full_Index()),B=analytic_velocity(c)->u(it.Location()/m,time/s)(it.Axis())*m/s;
                a=max(a,abs(A));
                b=max(b,abs(B));
                max_error=max(max_error,abs(A-B));
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,abs(A-B));}
            LOG::cout<<"max_error "<<max_error<<"  "<<a<<"  "<<b<<std::endl;}
    }

    TV Dirichlet_Boundary_Condition(const TV& X,int bc_color,int fluid_color,T time) PHYSBAM_OVERRIDE
    {
        Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
        if(analytic_velocity.m) return analytic_velocity(fluid_color)->u(X/m,time/s)*m/s;
        return TV();
    }

    TV Neumann_Boundary_Condition(const TV& X,int bc_color,int fluid_color,T time) PHYSBAM_OVERRIDE
    {
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));
        if(analytic_velocity.m && analytic_levelset){
            MATRIX<T,2> du=analytic_velocity(fluid_color)->du(X/m,time/s)/s;
            TV n=analytic_levelset->N(X/m,time/s);
            T p=analytic_velocity(fluid_color)->p(X/m,time/s)*kg/(s*s*m);
            return (du+du.Transposed())*n*mu(fluid_color)-p*n;}
        return TV();
    }

    TV Jump_Interface_Condition(const TV& X,int bc_color,int fluid_color,T time) PHYSBAM_OVERRIDE
    {
        return TV();
    }

//#####################################################################
};
}

#endif
