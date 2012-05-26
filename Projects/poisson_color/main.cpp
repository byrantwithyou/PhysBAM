//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Krylov_Solvers/MINRES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Finite_Elements/ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_COLOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>

using namespace PhysBAM;

typedef float RW;
std::string output_directory;

template<class TV>
GRID<TV>* Global_Grid(GRID<TV>* grid_in=0)
{
    static GRID<TV>* grid=0;
    GRID<TV>* old_grid=grid;
    if(grid_in) grid=grid_in;
    return old_grid;
}

typedef VECTOR<double,3> TV3;
TV3 color_map[3]={TV3(0,0.7,0),TV3(0.85,0.5,0),TV3(0,0.4,1)};

template<class TV> struct ANALYTIC_TEST;

//#################################################################################################################################################
// Debug Particles ################################################################################################################################
//#################################################################################################################################################

template<class TV> DEBUG_PARTICLES<TV>& Get_Debug_Particles()
{
    static DEBUG_PARTICLES<TV> debug_particles;
    return debug_particles;
}

template<class T, class TV>
void Dump_Frame(const ARRAY<T,FACE_INDEX<TV::m> >& u,const char* title)
{
    static int frame=0;
    char buff[100];
    sprintf(buff, "%s/%i", output_directory.c_str(), frame);
    FILE_UTILITIES::Create_Directory(buff);
    FILE_UTILITIES::Write_To_File<RW>((std::string)buff+"/mac_velocities.gz",u);
    if(title) FILE_UTILITIES::Write_To_Text_File((std::string)buff+"/frame_title",title);
    Get_Debug_Particles<TV>().Write_Debug_Particles(STREAM_TYPE((RW())),output_directory,frame);
    frame++;
}

template<class T,class TV>
void Flush_Frame(const char* title)
{
    Dump_Frame<T,TV>(ARRAY<T,FACE_INDEX<TV::m> >(*Global_Grid<TV>()),title);
}

template<class T,class TV>
void Dump_Interface(const INTERFACE_POISSON_SYSTEM_COLOR<TV>& ips)
{
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::SURFACE_ELEMENT SURFACE_ELEMENT;

    for(int i=0;i<ips.cdi->surface_mesh.m;i++){
        SURFACE_ELEMENT& V=ips.cdi->surface_mesh(i);
        if((V.z)>=0){
            if((V.z)>=0) Add_Debug_Object(V.x.X-V.x.Normal()*(T).03*ips.grid.dX.Min(),color_map[V.z]);
            if("#############################") Add_Debug_Object(V.x.X+V.x.Normal()*(T).03*ips.grid.dX.Min(),color_map[V.y]);}
        else if((V.y)>=0) Add_Debug_Object(V.x.X-V.x.Normal()*(T).03*ips.grid.dX.Min(),color_map[V.y]);}
}

template<class T,class TV>
void Dump_System(const INTERFACE_POISSON_SYSTEM_COLOR<TV>& ips,ANALYTIC_TEST<TV>& at)
{
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::SURFACE_ELEMENT SURFACE_ELEMENT;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::T_FACE T_FACE;
    
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(ips.grid);it.Valid();it.Next())
        Add_Debug_Particle(it.Location(),color_map[at.phi_color(it.Location())]);
    Flush_Frame<T,TV>("level set");

    Dump_Interface<T,TV>(ips);
    Flush_Frame<T,TV>("surfaces");
    
    char buff[100];
    for(int c=0;c<ips.cdi->colors;c++){
        Dump_Interface<T,TV>(ips);
        sprintf(buff,"dofs u%d",c);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(ips.grid);it.Valid();it.Next()){
            int index=ips.cm_u->Get_Index(it.index,c);
            if(index>=0){
                bool duplicated=false;
                for(int c_other=0;c_other<ips.cdi->colors;c_other++){
                    if(c_other==c) continue;
                    if(ips.cm_u->Get_Index(it.index,c_other)>=0) duplicated=true;}
                if(duplicated) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,1));
                else Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));}}                        
        Flush_Frame<T,TV>(buff);}

    for(int i=0;i<ips.cdi->surface_mesh.m;i++){
        SURFACE_ELEMENT& V=ips.cdi->surface_mesh(i);
        Add_Debug_Particle(V.x.Center(),VECTOR<T,3>(0,0.1,0.5));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,at.f_surface(V.x.Center(),V.y,V.z)*T_FACE::Normal(V.x.X));}
    Flush_Frame<T,TV>("surface forces");
    
    Dump_Interface<T,TV>(ips);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(ips.grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        T f_volume=at.f_volume(it.Location(),c);
        Add_Debug_Particle(it.Location(),f_volume==0?VECTOR<T,3>(0.25,0.25,0.25):(f_volume>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0)));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,f_volume);}
        Flush_Frame<T,TV>("volumetric forces");
}

template<class T,class TV>
void Dump_Vector(const INTERFACE_POISSON_SYSTEM_COLOR<TV>& ips,const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>& v,const char* title)
{
    char buff[100];
    for(int c=0;c<ips.cdi->colors;c++){
        Dump_Interface<T,TV>(ips);
        sprintf(buff,"%s u%d",title,c);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(ips.grid);it.Valid();it.Next()){
            int k=ips.cm_u->Get_Index(it.index,c);
            if(k>=0){
                T u_value=v.u(c)(k);
                Add_Debug_Particle(it.Location(),u_value==0?VECTOR<T,3>(0.25,0.25,0.25):(u_value>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0)));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,u_value);}}
        Flush_Frame<T,TV>(buff);}
}

template<class T,class TV>
void Dump_Vector(const INTERFACE_POISSON_SYSTEM_COLOR<TV>& ips,const ARRAY<T,VECTOR<int,TV::m> >& u,const char* title)
{
    INTERPOLATED_COLOR_MAP<T> cm;
    cm.Initialize_Colors(1e-12,1,true,true,true);

    Dump_Interface<T,TV>(ips);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(ips.grid);it.Valid();it.Next()){
        T u_value=u(it.index);
        Add_Debug_Particle(it.Location(),u_value==0?VECTOR<T,3>(0.25,0.25,0.25):(u_value>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0)));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,u_value);}
    Flush_Frame<T,TV>(title);
}

//#################################################################################################################################################
// Analytic Test ##################################################################################################################################
//#################################################################################################################################################

template<class TV>
struct ANALYTIC_TEST: public ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>
{
    typedef typename TV::SCALAR T;
    using ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>::kg;
    using ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>::m;
    using ANALYTIC_BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>::s;

    bool wrap;
    ARRAY<T> mu;

    virtual void Initialize()=0;
    virtual T phi_value(const TV& X)=0;
    virtual int phi_color(const TV& X)=0;
    virtual T u(const TV& X,int color)=0;
    virtual T f_volume(const TV& X,int color)=0;

    T u(const TV& X){return u(X,phi_color(X));}
};

template<class TV>
void Analytic_Test(GRID<TV>& grid,ANALYTIC_TEST<TV>& at,int max_iter,bool use_preconditioner,bool null,bool dump_matrix,bool debug_particles)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    ARRAY<T,TV_INT> phi_value(grid.Node_Indices());
    ARRAY<int,TV_INT> phi_color(grid.Node_Indices());

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        phi_value(it.index)=at.phi_value(it.Location());
        phi_color(it.index)=at.phi_color(it.Location());}
    
    INTERFACE_POISSON_SYSTEM_COLOR<TV> ips(grid,phi_value,phi_color);
    ips.use_preconditioner=use_preconditioner;
    ips.Set_Matrix(at.mu,at.wrap,&at);

    printf("\n");
    for(int c=0;c<ips.cdi->colors;c++) printf("u%d [%i]\t",c,ips.cm_u->dofs(c));printf("\n");
    printf("q [%i] ",ips.cdi->constraint_base_scalar);
    printf("\n");

    INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV> rhs,sol;
    {
        ARRAY<ARRAY<T,TV_INT> > f_volume;
        ARRAY<ARRAY<T,TV_INT> > u;
        
        f_volume.Resize(ips.cdi->colors);
        u.Resize(ips.cdi->colors);
        
        for(int c=0;c<ips.cdi->colors;c++) f_volume(c).Resize(grid.Domain_Indices());
        for(int c=0;c<ips.cdi->colors;c++) u(c).Resize(grid.Domain_Indices());
        
        for(int c=0;c<ips.cdi->colors;c++)
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
                f_volume(c)(it.index)=at.f_volume(it.Location(),c);
                u(c)(it.index)=at.u(it.Location(),c);}
        
        ips.Set_RHS(rhs,f_volume,u);
        ips.Resize_Vector(sol);
    }

    MINRES<T> mr;
    KRYLOV_SOLVER<T>* solver=&mr;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    solver->Solve(ips,sol,rhs,vectors,1e-10,0,max_iter);
    
    ips.Multiply(sol,*vectors(0));
    *vectors(0)-=rhs;
    LOG::cout<<"Residual: "<<ips.Convergence_Norm(*vectors(0))<<std::endl;

    ips.Multiply(ips.null_u,*vectors(0));
    LOG::cout<<"null u "<<ips.Convergence_Norm(*vectors(0))<<std::endl;

    ARRAY<T,TV_INT> exact_u,numer_u,error_u;

    numer_u.Resize(ips.grid.Domain_Indices());
    exact_u.Resize(ips.grid.Domain_Indices());
    error_u.Resize(ips.grid.Domain_Indices());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        int k=ips.cm_u->Get_Index(it.index,c);
        assert(k>=0);
        numer_u(it.index)=sol.u(c)(k);}

    T avg_u=0;
    int cnt_u=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        exact_u(it.index)=at.u(it.Location());
        error_u(it.index)=numer_u(it.index)-exact_u(it.index);
        avg_u+=error_u(it.index);
        cnt_u++;}
    avg_u/=cnt_u;

    T error_u_linf=0,error_u_l2=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        error_u(it.index)-=avg_u;
        error_u_linf=max(error_u_linf,abs(error_u(it.index)));
        error_u_l2+=sqr(error_u(it.index));}
    error_u_l2=sqrt(error_u_l2/cnt_u);

    LOG::cout<<ips.grid.counts<<" U error:   linf "<<error_u_linf<<"   l2 "<<error_u_l2<<std::endl<<std::endl;

    if(debug_particles){
        Dump_System<T,TV>(ips,at);
        Dump_Vector<T,TV>(ips,sol,"solution");
        Dump_Vector<T,TV>(ips,error_u,"error");}
        
    if(dump_matrix) OCTAVE_OUTPUT<T>("M.txt").Write("M",ips,*vectors(0),*vectors(1));
}

//#################################################################################################################################################
// Integration Test ###############################################################################################################################
//#################################################################################################################################################

template<class TV>
void Integration_Test(int argc,char* argv[],PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    Get_Debug_Particles<TV>().debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);

    int test_number;
    if(parse_args.Num_Extra_Args()<1){LOG::cerr<<"Test number is required."<<std::endl; exit(-1);}
    if(!STRING_UTILITIES::String_To_Value(parse_args.Extra_Arg(0),test_number)) throw VALUE_ERROR("The argument is not an integer.");

    ANALYTIC_TEST<TV>* test=0;

    switch(test_number){
        case 0:{ // One color, periodic. No interface, no forces, u=0.
            struct ANALYTIC_TEST_0:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);}
                virtual T phi_value(const TV& X){return (T)1;}
                virtual int phi_color(const TV& X){return T();}
                virtual T u(const TV& X,int color){return T();}
                virtual T f_volume(const TV& X,int color){return T();}
                virtual T f_surface(const TV& X,int color0,int color1){return T();}
            };
            test=new ANALYTIC_TEST_0;
            break;}
        case 1:{ // One color, periodic. No interface, u=sin(2*pi*x), f=(2*pi)^2*sin(2*pi*x).
            struct ANALYTIC_TEST_1:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);}
                virtual T phi_value(const TV& X){return (T)1;}
                virtual int phi_color(const TV& X){return T();}
                virtual T u(const TV& X,int color){return sin(2*M_PI*X.x);}
                virtual T f_volume(const TV& X,int color){return sqr(2*M_PI)*sin(2*M_PI*X.x)*mu(color);}
                virtual T f_surface(const TV& X,int color0,int color1){return T();}
            };
            test=new ANALYTIC_TEST_1;
            break;}
        case 2:{ // Two colors, periodic. Linear on [0,1/3],[1/3,2/3],[2/3,1], no volumetric forces.
            struct ANALYTIC_TEST_2:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);}
                virtual T phi_value(const TV& X){return abs(-m/(T)6+abs(X.x-0.5*m));}
                virtual int phi_color(const TV& X){return (-m/(T)6+abs(X.x-0.5*m))<0;}
                virtual T u(const TV& X,int color){return (color?(2*X.x-m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T f_volume(const TV& X,int color){return T();}
                virtual T f_surface(const TV& X,int color0,int color1){return ((X.x>0.5*m)?(T)(-1):(T)1)*(2*mu(1)+mu(0))/s;}
            };
            test=new ANALYTIC_TEST_2;
            break;}
        case 3:{ // Two colors, periodic. u=x^2 for r<R, zero elsewhere.
            struct ANALYTIC_TEST_3:public ANALYTIC_TEST<TV>
            {
                T r;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);r=m/M_PI;}
                virtual T phi_value(const TV& X){return abs((X-0.5*m).Magnitude()-r);}
                virtual int phi_color(const TV& X){return ((X-0.5*m).Magnitude()-r)<0;}
                virtual T u(const TV& X,int color){return (X-0.5*m).Magnitude_Squared()*color;}
                virtual T f_volume(const TV& X,int color){return -(TV::m)*2*mu(color)*color;}
                virtual T f_surface(const TV& X,int color0,int color1){return (X-0.5*m).Magnitude()*(-2)*mu(1);}
            };
            test=new ANALYTIC_TEST_3;
            break;}
        case 4:{ // Two colors, periodic. u=exp(-x^2) for r<R, zero elsewhere.
            struct ANALYTIC_TEST_4:public ANALYTIC_TEST<TV>
            {
                T r,m2,m4;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);r=m/M_PI;m2=sqr(m);m4=sqr(m2);}
                virtual T phi_value(const TV& X){return abs((X-0.5*m).Magnitude()-r);}
                virtual int phi_color(const TV& X){return ((X-0.5*m).Magnitude()-r)<0;}
                virtual T u(const TV& X,int color){return exp(-(X-0.5*m).Magnitude_Squared()/m2)*color;}
                virtual T f_volume(const TV& X,int color){T x2=(X-0.5*m).Magnitude_Squared(); return exp(-x2/m2)*(2*TV::m/m2-x2*4/m4)*mu(color)*color;}
                virtual T f_surface(const TV& X,int color0,int color1){T x2=(X-0.5*m).Magnitude_Squared(); return exp(-x2/m2)*sqrt(x2)*2*mu(1)/m2;}
            };
            test=new ANALYTIC_TEST_4;
            break;}
        case 5:{ // Three colors, periodic. u=a for r<R and x>0, u=b for r<R and x<0, zero elsewhere.
            struct ANALYTIC_TEST_5:public ANALYTIC_TEST<TV>
            {
                T r;
                TV n;
                VECTOR<T,3> a;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize()
                {
                    wrap=true;mu.Append(1);mu.Append(2);mu.Append(3);
                    r=m/M_PI;a(0)=0;a(1)=5;a(2)=7;
                    for(int i=0;i<TV::m;i++) n(i)=i+M_PI/(i+M_PI);n.Normalize();
                }
                virtual T phi_value(const TV& X){TV x=X-0.5*m;T s=x.Magnitude()-r;return (s<0)?min(abs(s),abs(x.Dot(n))):abs(s);}
                virtual int phi_color(const TV& X){TV x=X-0.5*m;return (x.Magnitude()-r)<0?((x.Dot(n)<0)?2:1):0;}
                virtual T u(const TV& X,int color){return a(color);}
                virtual T f_volume(const TV& X,int color){return T();}
                virtual T f_surface(const TV& X,int color0,int color1){return T();}
            };
            test=new ANALYTIC_TEST_5;
            break;}
        case 6:{ // Three colors, periodic. u=a*x^2 for r<R and x*n>0, u=b*x^2 for r<R and x*n<0, zero elsewhere.
            struct ANALYTIC_TEST_6:public ANALYTIC_TEST<TV>
            {
                T r;
                TV n;
                VECTOR<T,3> a;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize()
                {
                    wrap=true;mu.Append(1);mu.Append(2);mu.Append(3);
                    r=m/M_PI;a(0)=0;a(1)=5;a(2)=7;
                    for(int i=0;i<TV::m;i++) n(i)=i+M_PI/(i+M_PI);n.Normalize();
                }
                virtual T phi_value(const TV& X){TV x=X-0.5*m;T s=x.Magnitude()-r;return (s<0)?min(abs(s),abs(x.Dot(n))):abs(s);}
                virtual int phi_color(const TV& X){TV x=X-0.5*m;return (x.Magnitude()-r)<0?((x.Dot(n)<0)?2:1):0;}
                virtual T u(const TV& X,int color){return (X-0.5*m).Magnitude_Squared()*a(color);}
                virtual T f_volume(const TV& X,int color){return -(TV::m)*2*mu(color)*a(color);}
                virtual T f_surface(const TV& X,int color0,int color1){if(color0==0) return (X-0.5*m).Magnitude()*(-2)*mu(color1)*a(color1);else return T();}
            };
            test=new ANALYTIC_TEST_6;
            break;}
        case 7:{ // Three colors, periodic. u=a*exp(-x^2) for r<R and x*n>0, u=b*exp(-x^2) for r<R and x*n<0, zero elsewhere.
            struct ANALYTIC_TEST_7:public ANALYTIC_TEST<TV>
            {
                T r,a1,a2,m2,m4;
                TV n;
                VECTOR<T,3> a;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize()
                {
                    wrap=true;mu.Append(1);mu.Append(2);mu.Append(3);
                    r=m/M_PI;a(0)=0;a(1)=5;a(2)=7;m2=sqr(m);m4=sqr(m2);
                    for(int i=0;i<TV::m;i++) n(i)=i+M_PI/(i+M_PI);n.Normalize();
                }
                virtual T phi_value(const TV& X){TV x=X-0.5*m;T s=x.Magnitude()-r;return (s<0)?min(abs(s),abs(x.Dot(n))):abs(s);}
                virtual int phi_color(const TV& X){TV x=X-0.5*m;return (x.Magnitude()-r)<0?((x.Dot(n)<0)?2:1):0;}
                virtual T u(const TV& X,int color){return exp(-(X-0.5*m).Magnitude_Squared()/m2)*a(color);}
                virtual T f_volume(const TV& X,int color){T x2=(X-0.5*m).Magnitude_Squared(); return exp(-x2/m2)*(2*TV::m/m2-x2*4/m4)*mu(color)*a(color);}
                virtual T f_surface(const TV& X,int color0,int color1){T x2=(X-0.5*m).Magnitude_Squared();if(color0==0) return exp(-x2/m2)*2*mu(color1)*a(color1)*sqrt(x2)/m2;else return T();}
            };
            test=new ANALYTIC_TEST_7;
            break;}
        default:{
        LOG::cerr<<"Unknown test number."<<std::endl; exit(-1); break;}}

    output_directory=parse_args.Get_String_Value("-o");
    test->m=parse_args.Get_Double_Value("-m");
    test->s=parse_args.Get_Double_Value("-sec");
    test->kg=parse_args.Get_Double_Value("-kg");
    int res=parse_args.Get_Integer_Value("-resolution");
    int max_iter=parse_args.Get_Integer_Value("-max_iter");
    bool use_preconditioner=parse_args.Get_Option_Value("-use_preconditioner");
    bool null=parse_args.Get_Option_Value("-null");
    bool dump_matrix=parse_args.Get_Option_Value("-dump_matrix");
    bool debug_particles=parse_args.Get_Option_Value("-debug_particles");
    test->Initialize();

    TV_INT counts=TV_INT()+res;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1)*test->m,true);

    Global_Grid(&grid);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    Analytic_Test(grid,*test,max_iter,use_preconditioner,null,dump_matrix,debug_particles);
}

//#################################################################################################################################################
// Main ###########################################################################################################################################
//#################################################################################################################################################

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    PARSE_ARGS parse_args;
    parse_args.Set_Extra_Arguments(-1,"<example number>");
    parse_args.Add_String_Argument("-o","output","output directory");
    parse_args.Add_Double_Argument("-m",1,"meter scale");
    parse_args.Add_Double_Argument("-sec",1,"second scale");
    parse_args.Add_Double_Argument("-kg",1,"kilogram scale");
    parse_args.Add_Integer_Argument("-test",1,"test number");
    parse_args.Add_Integer_Argument("-resolution",4,"resolution");
    parse_args.Add_Integer_Argument("-max_iter",1000000,"max number of interations");
    parse_args.Add_Option_Argument("-use_preconditioner","use Jacobi preconditioner");
    parse_args.Add_Option_Argument("-3d","use 3D");
    parse_args.Add_Option_Argument("-null","find extra null modes of the matrix");
    parse_args.Add_Option_Argument("-dump_matrix","dump system matrix");
    parse_args.Add_Option_Argument("-debug_particles","dump debug particles");
    parse_args.Parse(argc,argv);

    if(parse_args.Get_Option_Value("-3d"))
        Integration_Test<VECTOR<double,3> >(argc,argv,parse_args);
    else
        Integration_Test<VECTOR<double,2> >(argc,argv,parse_args);

    return 0;
}
