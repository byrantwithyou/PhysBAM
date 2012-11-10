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
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/VOLUME_FORCE_SCALAR_COLOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <iostream>
#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace PhysBAM;

typedef float RW;
std::string output_directory="output";

template<class TV>
GRID<TV>* Global_Grid(GRID<TV>* grid_in=0)
{
    static GRID<TV>* grid=0;
    GRID<TV>* old_grid=grid;
    if(grid_in) grid=grid_in;
    return old_grid;
}

typedef VECTOR<double,3> TV3;
TV3 color_map[4]={TV3(0,0.7,0),TV3(0.8,0.8,0),TV3(0,0.4,1),TV3(0.8,0.2,0)};

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
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;

    for(int i=0;i<ips.cdi->surface_mesh.m;i++){
        INTERFACE_ELEMENT& V=ips.cdi->surface_mesh(i);
        if(V.color_pair.y>=0){
            if(V.color_pair.y>=0) Add_Debug_Object(V.face.X-V.face.Normal()*(T).003*ips.grid.dX.Min(),color_map[V.color_pair.y]);
            if("Alexey was here") Add_Debug_Object(V.face.X+V.face.Normal()*(T).003*ips.grid.dX.Min(),color_map[V.color_pair.x]);}
        else if(V.color_pair.x>=0) Add_Debug_Object(V.face.X-V.face.Normal()*(T).003*ips.grid.dX.Min(),color_map[V.color_pair.x]);}
}

template<class T,class TV>
void Dump_System(const INTERFACE_POISSON_SYSTEM_COLOR<TV>& ips,ANALYTIC_TEST<TV>& at)
{
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
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
        INTERFACE_ELEMENT& V=ips.cdi->surface_mesh(i);
        Add_Debug_Particle(V.face.Center(),VECTOR<T,3>(0,0.1,0.5));
        T k=0;
        if(V.color_pair.x>=0) k=at.j_surface(V.face.Center(),V.color_pair.x,V.color_pair.y);
        else if(V.color_pair.x==-1) k=at.n_surface(V.face.Center(),V.color_pair.x,V.color_pair.y);
        else if(V.color_pair.x==-2) k=at.d_surface(V.face.Center(),V.color_pair.x,V.color_pair.y);
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,k*T_FACE::Normal(V.face.X));}
    Flush_Frame<T,TV>("surface forces");
    
    Dump_Interface<T,TV>(ips);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(ips.grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        if(c>=0){
            T f_volume=at.f_volume(it.Location(),c);
            Add_Debug_Particle(it.Location(),f_volume==0?VECTOR<T,3>(0.25,0.25,0.25):(f_volume>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0)));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,f_volume);}}
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
void Analytic_Test(GRID<TV>& grid,ANALYTIC_TEST<TV>& at,int max_iter,bool use_preconditioner,bool null,bool dump_matrix,bool dump_geometry,bool debug_particles,bool double_fine)
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
    ips.Set_Matrix(at.mu,at.wrap,&at,double_fine);

    printf("\n");
    for(int c=0;c<ips.cdi->colors;c++) printf("u%d [%i]\t",c,ips.cm_u->dofs(c));printf("\n");
    printf("q [%i] ",ips.cdi->constraint_base_scalar);
    printf("\n");

    INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV> rhs,sol;

    struct VOLUME_FORCE_SCALAR_COLOR_LOCAL: public VOLUME_FORCE_SCALAR_COLOR<TV>
    {
        ANALYTIC_TEST<TV>* at;
        virtual T F(const TV& X,int color){return at->f_volume(X,color);}
    } vfscl;
    vfscl.at=&at;
    ips.Set_RHS(rhs,&vfscl);
    ips.Resize_Vector(sol);

    MINRES<T> mr;
    KRYLOV_SOLVER<T>* solver=&mr;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    LOG::SCOPE("System solve");

    solver->Solve(ips,sol,rhs,vectors,1e-10,0,max_iter);
    
    ips.Multiply(sol,*vectors(0));
    *vectors(0)-=rhs;
    LOG::cout<<"Residual: "<<ips.Convergence_Norm(*vectors(0))<<std::endl;

    ips.Multiply(ips.null_u,*vectors(0));
    LOG::cout<<"D constraints: "<<(ips.cdi->dc_present?"yes":"no")<<std::endl;
    LOG::cout<<"N constraints: "<<(ips.cdi->nc_present?"yes":"no")<<std::endl;
    if(ips.cdi->dc_present) LOG::cout<<"No null modes"<<std::endl;
    else LOG::cout<<"Null mode residual: "<<ips.Convergence_Norm(*vectors(0))<<std::endl;

    ARRAY<T,TV_INT> exact_u,numer_u,error_u;

    numer_u.Resize(ips.grid.Domain_Indices());
    exact_u.Resize(ips.grid.Domain_Indices());
    error_u.Resize(ips.grid.Domain_Indices());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        if(c>=0){
            int k=ips.cm_u->Get_Index(it.index,c);
            assert(k>=0);
            numer_u(it.index)=sol.u(c)(k);}}

    T avg_u=0;
    int cnt_u=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        if(c>=0){
            exact_u(it.index)=at.u(it.Location());
            error_u(it.index)=numer_u(it.index)-exact_u(it.index);
            avg_u+=error_u(it.index);
            cnt_u++;}}
    avg_u/=cnt_u;

    T error_u_linf=0,error_u_l2=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        if(c>=0){
            error_u(it.index)-=avg_u;
            error_u_linf=max(error_u_linf,abs(error_u(it.index)));
            error_u_l2+=sqr(error_u(it.index));}}
    error_u_l2=sqrt(error_u_l2/cnt_u);

    LOG::cout<<"error "<<error_u_linf<<" "<<error_u_l2<<std::endl<<std::endl;

    if(debug_particles){
        Dump_System<T,TV>(ips,at);
        Dump_Vector<T,TV>(ips,sol,"solution");
        Dump_Vector<T,TV>(ips,error_u,"error");
        if(null&&ips.Nullspace_Check(rhs)){
            OCTAVE_OUTPUT<T>("n.txt").Write("n",rhs);
            ips.Multiply(rhs,*vectors(0));
            LOG::cout<<"Extra nullspace found: "<<sqrt(ips.Inner_Product(*vectors(0),*vectors(0)))<<std::endl;
            rhs*=1/rhs.Max_Abs();
            Dump_Vector<T,TV>(ips,rhs,"extra null mode");}}
    
    if(dump_geometry){
        std::ofstream fout("geometry.txt");
        fout<<TV::m<<std::endl;
        for(int i=0;i<TV::m;i++) fout<<grid.counts(i)<<((i==TV::m-1)?"\n":" ");
        for(int i=0;i<TV::m;i++) fout<<grid.domain.min_corner(i)<<((i==TV::m-1)?"\n":" ");
        for(int i=0;i<TV::m;i++) fout<<grid.domain.max_corner(i)<<((i==TV::m-1)?"\n":" ");
        for(int k=0;k<ips.cdi->surface_mesh.m;k++){
            typename CELL_DOMAIN_INTERFACE_COLOR<TV>::INTERFACE_ELEMENT& se=ips.cdi->surface_mesh(k);
            for(int i=0;i<TV::m;i++) for(int j=0;j<TV::m;j++) fout<<se.face.X(i)(j)<<" ";
            fout<<se.color_pair.x<<" "<<se.color_pair.y<<std::endl;}
        fout.close();}
        
    if(dump_matrix) OCTAVE_OUTPUT<T>("M.txt").Write("M",ips,*vectors(0),*vectors(1));
    vectors.Delete_Pointers_And_Clean_Memory();
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

    T m=1,s=1,kg=1;
    int threads=1;
    int test_number=1,resolution=4,max_iter=1000000;
    bool use_preconditioner=false,use_test=false,null=false,dump_matrix=false,debug_particles=false,double_fine=false,dump_geometry=false,opt_arg=false;
    parse_args.Extra_Optional(&test_number,&opt_arg,"example number","example number to run");
    parse_args.Add("-o",&output_directory,"output","output directory");
    parse_args.Add("-m",&m,"unit","meter scale");
    parse_args.Add("-s",&s,"unit","second scale");
    parse_args.Add("-kg",&kg,"unit","kilogram scale");
    parse_args.Add("-test",&test_number,&use_test,"number","test number");
    parse_args.Add("-resolution",&resolution,"res","resolution");
    parse_args.Add("-threads",&threads,"threads","number of threads");
    parse_args.Add("-use_preconditioner",&use_preconditioner,"Use Jacobi preconditioner");
    parse_args.Add("-max_iter",&max_iter,"iter","max number of interations");
    parse_args.Add("-null",&null,"find extra null modes of the matrix");
    parse_args.Add("-dump_matrix",&dump_matrix,"dump system matrix");
    parse_args.Add("-debug_particles",&debug_particles,"dump debug particles");
    parse_args.Add("-double_fine",&double_fine,"set level set exactly on double fine grid");
    parse_args.Add("-dump_geometry",&dump_geometry,"dump grid info and interface");
    parse_args.Parse();

#ifdef USE_OPENMP
    omp_set_num_threads(threads);
#pragma omp parallel
#pragma omp master
    {
        PHYSBAM_ASSERT(threads==omp_get_num_threads());
        LOG::cout<<"Running on "<<threads<<" threads"<<std::endl;
    }
#endif

    if(!use_test && !opt_arg){
        LOG::cerr<<"Test number is required."<<std::endl;
        exit(-1);}

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
                virtual T j_surface(const TV& X,int color0,int color1){return T();}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
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
                virtual T j_surface(const TV& X,int color0,int color1){return T();}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
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
                virtual T j_surface(const TV& X,int color0,int color1){return ((X.x>0.5*m)?(T)(-1):(T)1)*(2*mu(1)+mu(0))/s;}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
            };
            test=new ANALYTIC_TEST_2;
            break;}
        case 3:{ // Two colors, periodic. u=||x||^2 for r<R, zero elsewhere.
            struct ANALYTIC_TEST_3:public ANALYTIC_TEST<TV>
            {
                T r;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){
                    wrap=true;mu.Append(1);mu.Append(2);r=m/M_PI;
                    this->use_discontinuous_scalar_field=true;}
                virtual T phi_value(const TV& X){return abs((X-0.5*m).Magnitude()-r);}
                virtual int phi_color(const TV& X){return ((X-0.5*m).Magnitude()-r)<0;}
                virtual T u(const TV& X,int color){return (X-0.5*m).Magnitude_Squared()*color;}
                virtual T f_volume(const TV& X,int color){return -(TV::m)*2*mu(color)*color;}
                virtual T j_surface(const TV& X,int color0,int color1){return (X-0.5*m).Magnitude()*(-2)*mu(1);}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
            };
            test=new ANALYTIC_TEST_3;
            break;}
        case 4:{ // Two colors, periodic. u=exp(-x^2) for r<R, zero elsewhere.
            struct ANALYTIC_TEST_4:public ANALYTIC_TEST<TV>
            {
                T r,m2,m4;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){
                    wrap=true;mu.Append(1);mu.Append(2);r=m/M_PI;m2=sqr(m);m4=sqr(m2);
                    this->use_discontinuous_scalar_field=true;}
                virtual T phi_value(const TV& X){return abs((X-0.5*m).Magnitude()-r);}
                virtual int phi_color(const TV& X){return ((X-0.5*m).Magnitude()-r)<0;}
                virtual T u(const TV& X,int color){return exp(-(X-0.5*m).Magnitude_Squared()/m2)*color;}
                virtual T f_volume(const TV& X,int color){T x2=(X-0.5*m).Magnitude_Squared(); return exp(-x2/m2)*(2*TV::m/m2-x2*4/m4)*mu(color)*color;}
                virtual T j_surface(const TV& X,int color0,int color1){T x2=(X-0.5*m).Magnitude_Squared(); return exp(-x2/m2)*sqrt(x2)*2*mu(1)/m2;}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
            };
            test=new ANALYTIC_TEST_4;
            break;}
        case 5:{ // Three colors, periodic. Stripes in x 0:[0,a], 1:[a,b], 2:[b,c], 0:[c,1].
            struct ANALYTIC_TEST_5:public ANALYTIC_TEST<TV>
            {
                T a,b,c;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize()
                {
                    wrap=true;mu.Append(1);mu.Append(2);mu.Append(3);
                    a=m/6;b=m*5/12;c=m*5/6;
                    this->use_discontinuous_scalar_field=true;
                }
                virtual T phi_value(const TV& X)
                {
                    if(X.x<a) return abs(X.x-a);
                    if(X.x<b) return (b-a)/2-abs(X.x-(b+a)/2);
                    if(X.x<c) return (c-b)/2-abs(X.x-(c+b)/2);
                    return abs(X.x-c);
                }
                virtual int phi_color(const TV& X)
                {
                    if(X.x<a) return 0;
                    if(X.x<b) return 1;
                    if(X.x<c) return 2;
                    return 0;
                }
                virtual T u(const TV& X,int color){
                    switch (color){
                        case 0: return 0;
                        case 1: return exp(X.x/m);
                        case 2: return sin(2*M_PI*X.y/m);
                        default: PHYSBAM_FATAL_ERROR();}
                }
                virtual T f_volume(const TV& X,int color)
                {
                    switch (color){
                        case 0: return 0;
                        case 1: return -exp(X.x/m)*mu(color)/sqr(m);
                        case 2: return sqr(2*M_PI)*sin(2*M_PI*X.y/m)*mu(color)/sqr(m);
                        default: PHYSBAM_FATAL_ERROR();}
                }
                virtual T j_surface(const TV& X,int color0,int color1)
                {
                    if(color0==0 && color1==1) return exp(X.x/m)*mu(1)/m;
                    if(color0==0 && color1==2) return 0;
                    if(color0==1 && color1==2) return -exp(X.x/m)*mu(1)/m;
                    PHYSBAM_FATAL_ERROR();
                }
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
            };
            test=new ANALYTIC_TEST_5;
            break;}
        case 6:{ // Three colors (one dirichlet or neumann), periodic. Stripes in x 0:[0,a], 1:[a,b], 2:[b,c], 0:[c,1].
            struct ANALYTIC_TEST_6:public ANALYTIC_TEST<TV>
            {
                T a,b,c;
                int constraint; // [-1] - Neumann, [-2] - Dirichlet
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize()
                {
                    wrap=true;mu.Append(1);mu.Append(2);
                    a=m/6;b=m*5/12;c=m*5/6;
                    constraint=-1;
                    this->use_discontinuous_scalar_field=true;
                }
                virtual T phi_value(const TV& X)
                {
                    if(X.x<a) return abs(X.x-a);
                    if(X.x<b) return (b-a)/2-abs(X.x-(b+a)/2);
                    if(X.x<c) return (c-b)/2-abs(X.x-(c+b)/2);
                    return abs(X.x-c);
                }
                virtual int phi_color(const TV& X)
                {
                    if(X.x<a) return constraint;
                    if(X.x<b) return 0;
                    if(X.x<c) return 1;
                    return constraint;
                }
                virtual T u(const TV& X,int color){
                    switch (color){
                        case 0: return exp(X.x/m);
                        case 1: return sin(2*M_PI*X.y/m);
                        default: PHYSBAM_FATAL_ERROR();}
                }
                virtual T f_volume(const TV& X,int color)
                {
                    switch (color){
                        case 0: return -exp(X.x/m)*mu(color)/sqr(m);
                        case 1: return sqr(2*M_PI)*sin(2*M_PI*X.y/m)*mu(color)/sqr(m);
                        default: PHYSBAM_FATAL_ERROR();}
                }
                virtual T j_surface(const TV& X,int color0,int color1)
                {
                    return -exp(X.x/m)*mu(0)/m;
                }
                virtual T n_surface(const TV& X,int color0,int color1)
                {
                    PHYSBAM_ASSERT(constraint==-1);
                    if(color1==0) return exp(X.x/m)*mu(0)/m;
                    if(color1==1) return 0;
                    PHYSBAM_FATAL_ERROR();
                }
                virtual T d_surface(const TV& X,int color0,int color1)
                {
                    PHYSBAM_ASSERT(constraint==-2);
                    if(color1==0) return exp(X.x/m)*mu(0)/m;
                    if(color1==1) return 0;
                    PHYSBAM_FATAL_ERROR();
                }
            };
            test=new ANALYTIC_TEST_6;
            break;}
        case 7:{ // Three colors, periodic. u=a for r<R and x>0, u=b for r<R and x<0, zero elsewhere.
            struct ANALYTIC_TEST_7:public ANALYTIC_TEST<TV>
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
                    this->use_discontinuous_scalar_field=true;
                }
                virtual T phi_value(const TV& X){TV x=X-0.5*m;T s=x.Magnitude()-r;return (s<0)?min(abs(s),abs(x.Dot(n))):abs(s);}
                virtual int phi_color(const TV& X){TV x=X-0.5*m;return (x.Magnitude()-r)<0?((x.Dot(n)<0)?2:1):0;}
                virtual T u(const TV& X,int color){return a(color);}
                virtual T f_volume(const TV& X,int color){return T();}
                virtual T j_surface(const TV& X,int color0,int color1){return T();}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
            };
            test=new ANALYTIC_TEST_7;
            break;}
        case 8:{ // Three colors, periodic. u=a*x^2 for r<R and x*n>0, u=b*x^2 for r<R and x*n<0, zero elsewhere.
            struct ANALYTIC_TEST_8:public ANALYTIC_TEST<TV>
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
                    this->use_discontinuous_scalar_field=true;
                }
                virtual T phi_value(const TV& X){TV x=X-0.5*m;T s=x.Magnitude()-r;return (s<0)?min(abs(s),abs(x.Dot(n))):abs(s);}
                virtual int phi_color(const TV& X){TV x=X-0.5*m;return (x.Magnitude()-r)<0?((x.Dot(n)<0)?2:1):0;}
                virtual T u(const TV& X,int color){return (X-0.5*m).Magnitude_Squared()*a(color);}
                virtual T f_volume(const TV& X,int color){return -(TV::m)*2*mu(color)*a(color);}
                virtual T j_surface(const TV& X,int color0,int color1){if(color0==0) return (X-0.5*m).Magnitude()*(-2)*mu(color1)*a(color1);else return T();}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
            };
            test=new ANALYTIC_TEST_8;
            break;}
        case 9:{ // Three colors, periodic. u=a*exp(-x^2) for r<R and x*n>0, u=b*exp(-x^2) for r<R and x*n<0, zero elsewhere.
            struct ANALYTIC_TEST_9:public ANALYTIC_TEST<TV>
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
                    this->use_discontinuous_scalar_field=true;
                }
                virtual T phi_value(const TV& X){TV x=X-0.5*m;T s=x.Magnitude()-r;return (s<0)?min(abs(s),abs(x.Dot(n))):abs(s);}
                virtual int phi_color(const TV& X){TV x=X-0.5*m;return (x.Magnitude()-r)<0?((x.Dot(n)<0)?2:1):0;}
                virtual T u(const TV& X,int color){return exp(-(X-0.5*m).Magnitude_Squared()/m2)*a(color);}
                virtual T f_volume(const TV& X,int color){T x2=(X-0.5*m).Magnitude_Squared(); return exp(-x2/m2)*(2*TV::m/m2-x2*4/m4)*mu(color)*a(color);}
                virtual T j_surface(const TV& X,int color0,int color1){T x2=(X-0.5*m).Magnitude_Squared();if(color0==0) return exp(-x2/m2)*2*mu(color1)*a(color1)*sqrt(x2)/m2;else return T();}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
            };
            test=new ANALYTIC_TEST_9;
            break;}
        case 10:{ // Three colors (dirichlet/neumann outside), periodic. u=a*exp(-x^2) for r<R and x*n>0, u=b*exp(-x^2) for r<R and x*n<0, zero elsewhere.
            struct ANALYTIC_TEST_10:public ANALYTIC_TEST<TV>
            {
                T r,a1,a2,m2,m4;
                TV n;
                VECTOR<T,3> a;
                int constraint;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize()
                {
                    wrap=true;mu.Append(1);mu.Append(2);
                    r=m/M_PI;a(0)=0;a(1)=5;a(2)=7;m2=sqr(m);m4=sqr(m2);
                    for(int i=0;i<TV::m;i++) n(i)=i+M_PI/(i+M_PI);n.Normalize();
                    constraint=-2;
                    this->use_discontinuous_scalar_field=true;
                }
                virtual T phi_value(const TV& X){TV x=X-0.5*m;T s=x.Magnitude()-r;return (s<0)?min(abs(s),abs(x.Dot(n))):abs(s);}
                virtual int phi_color(const TV& X){TV x=X-0.5*m;return (x.Magnitude()-r)<0?((x.Dot(n)<0)?1:0):constraint;}
                virtual T u(const TV& X,int color){return exp(-(X-0.5*m).Magnitude_Squared()/m2)*a(color);}
                virtual T f_volume(const TV& X,int color){T x2=(X-0.5*m).Magnitude_Squared(); return exp(-x2/m2)*(2*TV::m/m2-x2*4/m4)*mu(color)*a(color);}
                virtual T j_surface(const TV& X,int color0,int color1){return T();}
                virtual T d_surface(const TV& X,int color0,int color1){PHYSBAM_ASSERT(constraint==-2); return u(X,color1);}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
            };
            test=new ANALYTIC_TEST_10;
            break;}
        case 11:{ // Four colors (dirichlet outside). Three equal bubbles.
            struct ANALYTIC_TEST_11:public ANALYTIC_TEST<TV>
            {
                T r;
                TV n; // rotate angle with respect to e_y
                TV a; // shift
                VECTOR<TV,3> centers;
                VECTOR<TV,3> normals;
                VECTOR<VECTOR<int,3>,3> sectors;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize()
                {
                    wrap=true;mu.Append(1);mu.Append(2);mu.Append(3);
                    r=(T)1/(2*M_PI-1);
                    n=TV::Axis_Vector(1);a=TV()+M_PI/200;
                    centers(0)=TV::Axis_Vector(1);
                    centers(1).x=(T)sqrt(3)/2;centers(1).y=-(T)1/2;
                    centers(2).x=-(T)sqrt(3)/2;centers(2).y=-(T)1/2;
                    for(int i=0;i<3;i++){
                        normals(i).x=-centers(i).y;
                        normals(i).y=centers(i).x;
                        centers(i)*=r;
                        for(int j=0;j<3;j++) sectors(i)(j)=(i+j)%3;}
                    this->use_discontinuous_scalar_field=true;
                }
                virtual TV Transform(const TV& X){return X-0.5+a;}
                virtual T phi_value(const TV& X)
                {
                    TV x=Transform(X);
                    int i;
                    for(i=0;i<2;i++) if(x.Dot(normals(sectors(i)(1)))>=0 && x.Dot(normals(sectors(i)(2)))<0) break;
                    T d=(x-centers(i)).Magnitude();
                    if(d>r && x.Magnitude()>r/100) return d-r;
                    else return min(abs(d-r),abs(x.Dot(normals(sectors(i)(1)))),abs(x.Dot(normals(sectors(i)(2)))));
                }
                virtual int phi_color(const TV& X){
                    TV x=Transform(X);
                    int i;
                    for(i=0;i<2;i++) if(x.Dot(normals(sectors(i)(1)))>=0 && x.Dot(normals(sectors(i)(2)))<0) break;
                    T d=(x-centers(i)).Magnitude();
                    if(d>r && x.Magnitude()>r/100) return -2;
                    else return i;
                }
                virtual T u(const TV& X,int color)
                {
                    switch(color){
                        case 0: return sin(X.x);
                        case 1: return cos(X.y);
                        case 2: return X.Magnitude_Squared();
                        default: PHYSBAM_FATAL_ERROR();}
                }
                virtual TV grad_u(const TV& X,int color)
                {
                    switch(color){
                        case 0: return TV::Axis_Vector(0)*cos(X.x)*mu(color);
                        case 1: return -TV::Axis_Vector(1)*sin(X.y)*mu(color);
                        case 2: return (TV()+2)*X*mu(color);
                        default: PHYSBAM_FATAL_ERROR();}
                }
                virtual T f_volume(const TV& X,int color)
                {
                    switch(color){
                        case 0: return sin(X.x)*mu(color);
                        case 1: return cos(X.y)*mu(color);
                        case 2: return -2*TV::m*mu(color);
                        default: PHYSBAM_FATAL_ERROR();}
                }
                virtual T j_surface(const TV& X,int color0,int color1)
                {
                    TV x=Transform(X);
                    if(color0==-2){return -grad_u(X,color1).Dot((x-centers(color1)).Normalized());}
                    if(color0==0 && color1==1) return normals(2).Dot(grad_u(X,1)-grad_u(X,0));
                    if(color0==1 && color1==2) return normals(0).Dot(grad_u(X,2)-grad_u(X,1));
                    if(color0==0 && color1==2) return normals(1).Dot(grad_u(X,0)-grad_u(X,2));
                    PHYSBAM_FATAL_ERROR();
                }
                virtual T d_surface(const TV& X,int color0,int color1){return u(X,color1);}
                virtual T n_surface(const TV& X,int color0,int color1){PHYSBAM_FATAL_ERROR();}
            };
            test=new ANALYTIC_TEST_11;
            break;}
        default:{
        LOG::cerr<<"Unknown test number."<<std::endl; exit(-1); break;}}

    test->m=m;
    test->s=s;
    test->kg=kg;
    test->Initialize();

    TV_INT counts=TV_INT()+resolution;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1)*test->m,true);

    Global_Grid(&grid);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    Analytic_Test(grid,*test,max_iter,use_preconditioner,null,dump_matrix,dump_geometry,debug_particles,double_fine);
    LOG::Finish_Logging();
    delete test;
}

//#################################################################################################################################################
// Main ###########################################################################################################################################
//#################################################################################################################################################

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"Use 3D");
    parse_args.Parse(true);
    
    if(use_3d)
        Integration_Test<VECTOR<double,3> >(argc,argv,parse_args);
    else
        Integration_Test<VECTOR<double,2> >(argc,argv,parse_args);

    return 0;
}
