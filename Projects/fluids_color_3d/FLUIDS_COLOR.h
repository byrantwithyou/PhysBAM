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
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace PhysBAM{

template<class TV>
class FLUIDS_COLOR:public PLS_FC_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef PLS_FC_EXAMPLE<TV> BASE;

public:
    using BASE::grid;using BASE::output_directory;using BASE::domain_boundary;using BASE::face_velocities;
    using BASE::write_substeps_level;using BASE::restart;using BASE::last_frame;
    using BASE::dt;using BASE::levelset_color;using BASE::mu;using BASE::rho;using BASE::dump_matrix;using BASE::number_of_colors;
    using BASE::use_advection;using BASE::use_reduced_advection;using BASE::omit_solve;using BASE::use_discontinuous_velocity;
    using BASE::time_steps_per_frame;using BASE::use_p_null_mode;using BASE::Fill_Levelsets_From_Levelset_Color;

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
        virtual MATRIX<T,3> du(const TV& X,T t) const=0;
        virtual T p(const TV& X,T t) const=0;
        virtual TV F(const TV& X,T t) const=0;
    };
    ARRAY<ANALYTIC_VELOCITY*> analytic_velocity;

    struct ANALYTIC_LEVELSET
    {
        virtual T phi(const TV& X,T t,int& c) const=0;
        virtual TV N(const TV& X,T t,int c) const=0;
    }* analytic_levelset;

    FLUIDS_COLOR(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
        :PLS_FC_EXAMPLE<TV>(stream_type),test_number(0),resolution(32),stored_last_frame(0),user_last_frame(false),mu0(1),mu1(2),
        rho0(1),rho1(2),m(1),s(1),kg(1),bc_n(false),bc_d(false),bc_s(false),no_advection(false),refine(1)
    {
        last_frame=16;
        int number_of_threads=1;
        parse_args.Extra(&test_number,"example number","example number to run");
        parse_args.Add("-restart",&restart,"frame","restart frame");
        parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
        parse_args.Add("-substep",&write_substeps_level,"level","output-substep level");
        parse_args.Add("-dt",&dt,"dt","time step size to use for simulation");
        parse_args.Add("-steps",&time_steps_per_frame,"steps","number of time steps per frame");
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
        parse_args.Add("-null_p",&use_p_null_mode,"Assume pressure null mode and project it out");
        parse_args.Add("-threads",&number_of_threads,"threads","Number of threads");
        parse_args.Parse();

#ifdef USE_OPENMP
        omp_set_num_threads(number_of_threads);
#pragma omp parallel
#pragma omp single
        {
            if(omp_get_num_threads()!=number_of_threads) PHYSBAM_FATAL_ERROR();
            LOG::cout<<"Running on "<<number_of_threads<<" threads"<<std::endl;
        }
#endif
        
        resolution*=refine;
        dt/=refine;
        time_steps_per_frame*=refine;
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
                grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*m,true);
                analytic_levelset=new ANALYTIC_LEVELSET_SPHERE(TV()+(T).5,(T).3);
                analytic_velocity.Append(new ANALYTIC_VELOCITY_CONST(TV()+1));
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

    void Set_Level_Set(T time)
    {
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,1);it.Valid();it.Next()){
            int c=-4;
            T p=analytic_levelset->phi(it.Location()/m,time,c)*m;
            levelset_color.phi(it.index)=abs(p);
            levelset_color.color(it.index)=c==-4?bc_type:c;}
        Fill_Levelsets_From_Levelset_Color();
    }

    void Get_Initial_Velocities()
    {
        if(analytic_levelset && analytic_velocity.m)
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
                int c=face_color(it.Full_Index());
                if(c<0) continue;
                face_velocities(c)(it.Full_Index())=analytic_velocity(c)->u(it.Location()/m,0)(it.Axis())*m/s;}
    }

    void Analytic_Test()
    {
        mu.Append(mu0);
        rho.Append(rho0);
        for(int i=1;i<number_of_colors;i++){
            mu.Append(mu1);
            rho.Append(rho1);}

        Set_Level_Set(0);

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
            MATRIX<T,3> du0=analytic_velocity(c)->du(X/m,t)/s,du1=analytic_velocity(c)->du((X+dX)/m,t)/s;
            T erru=((du0+du1)*dX/2-(u1-u0)).Magnitude()/e;
            int c0,c1;
            T l0=analytic_levelset->phi(X/m,t,c0)*m,l1=analytic_levelset->phi((X+dX)/m,t,c1)*m;
            if(c0>=0) l0=-l0;
            if(c1>=0) l1=-l1;
            TV dl0=analytic_levelset->N(X/m,t,c),dl1=analytic_levelset->N((X+dX)/m,t,c);
            T errl=abs((dl0+dl1).Dot(dX)/2-(l1-l0))/e;
            LOG::cout<<"analytic diff test "<<erru<<"  "<<errl<<std::endl;}
    }

    struct ANALYTIC_VELOCITY_CONST:public ANALYTIC_VELOCITY
    {
        TV au;
        ANALYTIC_VELOCITY_CONST(TV v): au(v){}
        virtual TV u(const TV& X,T t) const {return au;}
        virtual MATRIX<T,3> du(const TV& X,T t) const {return MATRIX<T,3>();}
        virtual T p(const TV& X,T t) const {return 0;}
        virtual TV F(const TV& X,T t) const {return TV();}
    };

    struct ANALYTIC_VELOCITY_AFFINE:public ANALYTIC_VELOCITY
    {
        TV v0;
        MATRIX<T,3> du0;
        T rho;
        ANALYTIC_VELOCITY_AFFINE(const TV& x0,const TV& v0,const MATRIX<T,3>& du0,T rho): v0(v0-du0*x0),du0(du0) {}
        virtual TV u(const TV& X,T t) const {return du0*X+v0;}
        virtual MATRIX<T,3> du(const TV& X,T t) const {return du0;}
        virtual T p(const TV& X,T t) const {return 0;}
        virtual TV F(const TV& X,T t) const {return rho*du0*(du0*X+v0);}
    };

    struct ANALYTIC_VELOCITY_TRANSLATE:public ANALYTIC_VELOCITY
    {
        ANALYTIC_VELOCITY* av;
        TV vel;

        ANALYTIC_VELOCITY_TRANSLATE(ANALYTIC_VELOCITY* av,const TV& vel): av(av),vel(vel) {}
        ~ANALYTIC_VELOCITY_TRANSLATE() {delete av;}
        virtual TV u(const TV& X,T t) const {return av->u(X-vel*t,t)+vel;}
        virtual MATRIX<T,3> du(const TV& X,T t) const {return av->du(X-vel*t,t);}
        virtual T p(const TV& X,T t) const {return av->p(X-vel*t,t);}
        virtual TV F(const TV& X,T t) const {return av->F(X-vel*t,t);}
    };

    struct ANALYTIC_LEVELSET_PERIODIC:public ANALYTIC_LEVELSET
    {
        virtual T phi(const TV& X,T t,int& c) const {c=0;return 1;}
        virtual TV N(const TV& X,T t,int c) const {return TV(1,0,0);}
    };

    struct ANALYTIC_LEVELSET_SPHERE:public ANALYTIC_LEVELSET
    {
        TV cen;
        T r;
        ANALYTIC_LEVELSET_SPHERE(TV cc,T rr): cen(cc),r(rr){}
        virtual T phi(const TV& X,T t,int& c) const {T p=(X-cen).Magnitude()-r;c=p>=0?-4:0;return abs(p);}
        virtual TV N(const TV& X,T t,int c) const {return (X-cen).Normalized();}
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
            c=bc1;
            return p1;
        }
        virtual TV N(const TV& X,T t,int c) const
        {
            if(c>0 && c<x.m)
            {
                T p0=X.x-x(c-1),p1=X.x-x(c);
                return TV(p0<-p1?-1:1,0);
            }
            if(c==0)
            {
                T p0=X.x-x0,p1=X.x-x(0);
                return TV(p0<-p1?-1:1,0);
            }
            T p0=X.x-x.Last(),p1=X.x-x1;
            return TV(p0<-p1?-1:1,0);
        }
    };

    struct ANALYTIC_LEVELSET_CONCENTRIC:public ANALYTIC_LEVELSET
    {
        TV center;
        int bc0,bc1;
        T r0,r1;
        ARRAY<T> array;
        ANALYTIC_LEVELSET_CONCENTRIC(TV c,T r0,T r1,int bc0,int bc1): center(c),bc0(bc0),bc1(bc1),r0(r0),r1(r1) {}
        ANALYTIC_LEVELSET_CONCENTRIC* Add(T nr) {array.Append(nr);return this;}
        virtual T phi(const TV& X,T t,int& c) const
        {
            T r=(X-center).Magnitude(),p0=r-r0;
            if(p0<=0){c=bc0;return -p0;}
            for(int i=0;i<array.m;i++){
                T p1=r-array(i);
                if(p1<=0){c=i;return min(p0,-p1);}
                p0=p1;}
            T p1=r-r1;
            if(p1<=0){c=array.m;return min(p0,-p1);}
            c=bc1;
            return p1;
        }
        virtual TV N(const TV& X,T t,int c) const
        {
            TV n=X-center;
            T r=n.Normalize();
            if(c>0 && c<array.m)
            {
                T p0=r-array(c-1),p1=r-array(c);
                return p0<-p1?-n:n;
            }
            if(c==0)
            {
                T p0=r-r0,p1=r-array(0);
                return p0<-p1?-n:n;
            }
            T p0=r-array.Last(),p1=r-r1;
            return p0<-p1?-n:n;
        }
    };

    struct ANALYTIC_LEVELSET_TRANSLATE:public ANALYTIC_LEVELSET
    {
        ANALYTIC_LEVELSET* al;
        TV vel;

        ANALYTIC_LEVELSET_TRANSLATE(ANALYTIC_LEVELSET* al,const TV& vel): al(al),vel(vel) {}
        ~ANALYTIC_LEVELSET_TRANSLATE() {delete al;}
        virtual T phi(const TV& X,T t,int& c) const {return al->phi(X-vel*t,t,c);}
        virtual TV N(const TV& X,T t,int c) const {return al->N(X-vel*t,t,c);}
    };

    struct ANALYTIC_LEVELSET_SCALE:public ANALYTIC_LEVELSET
    {
        ANALYTIC_LEVELSET* al;
        T scale;

        ANALYTIC_LEVELSET_SCALE(ANALYTIC_LEVELSET* al,T scale): al(al),scale(scale) {}
        ~ANALYTIC_LEVELSET_SCALE() {delete al;}
        virtual T phi(const TV& X,T t,int& c) const {return (1+t*scale)*al->phi(X/(1+t*scale),t,c);}
        virtual TV N(const TV& X,T t,int c) const {return al->N(X/(1+t*scale),t,c);}
    };

    void Begin_Time_Step(const T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset)
            Set_Level_Set(time);
    }

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
            MATRIX<T,3> du=analytic_velocity(fluid_color)->du(X/m,time/s)/s;
            TV n=analytic_levelset->N(X/m,time/s,fluid_color);
            T p=analytic_velocity(fluid_color)->p(X/m,time/s)*kg/(s*s*m);
            return (du+du.Transposed())*n*mu(fluid_color)-p*n;}
        return TV();
    }

    TV Jump_Interface_Condition(const TV& X,int color0,int color1,T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset){
            MATRIX<T,3> du0=analytic_velocity(color0)->du(X/m,time/s)/s,du1=analytic_velocity(color1)->du(X/m,time/s)/s;
            T p0=analytic_velocity(color0)->p(X/m,time/s)*kg/(s*s*m),p1=analytic_velocity(color1)->p(X/m,time/s)*kg/(s*s*m);
            MATRIX<T,3> stress0=(du0+du0.Transposed())*mu(color0)-p0,stress1=(du1+du1.Transposed())*mu(color1)-p1;
            TV n=analytic_levelset->N(X/m,time/s,color1);
            return (stress1-stress0)*n;}
        return TV();
    }

    TV Volume_Force(const TV& X,int color,T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset)
            return analytic_velocity(color)->F(X/m,time/s)*kg/(m*s*s);
        return TV();
    }

    TV Velocity_Jump(const TV& X,int color0,int color1,T time) PHYSBAM_OVERRIDE
    {
        if(analytic_velocity.m && analytic_levelset)
            return (analytic_velocity(color1)->u(X/m,time/s)-analytic_velocity(color0)->u(X/m,time/s))*m/s;
        return TV();
    }

//#####################################################################
};
}

#endif
