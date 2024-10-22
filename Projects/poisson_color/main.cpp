//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_LINE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_NEST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SPHERE.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_COLOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/REINITIALIZATION.h>
#include "ANALYTIC_POISSON_TEST.h"
#include <iostream>
#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace PhysBAM;

typedef float RW;

typedef VECTOR<double,3> TV3;
TV3 color_map[7]={TV3(.5,.5,.5),TV3(.8,.8,.8),TV3(1,1,1),TV3(0,0.7,0),TV3(0.8,0.8,0),TV3(0,0.4,1),TV3(0.8,0.2,0)};

//#################################################################################################################################################
// Debug Particles ################################################################################################################################
//#################################################################################################################################################

template<class T,class TV>
void Dump_Interface(const INTERFACE_POISSON_SYSTEM_COLOR<TV>& ips)
{
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::CELL_ELEMENTS CELL_ELEMENTS;
    typedef VECTOR<int,TV::m> TV_INT;

    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(ips.cdi->index_to_cell_elements);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<INTERFACE_ELEMENT>& interface_elements=cell_elements.interface;
        for(int i=0;i<interface_elements.m;i++){
            const INTERFACE_ELEMENT& V=interface_elements(i);
            if(V.color_pair.y>=0){
                if(V.color_pair.y>=0) Add_Debug_Object(V.face.X-V.face.Normal()*(T).003*ips.grid.dX.Min(),color_map[V.color_pair.y+3]);
                if("Alexey was here") Add_Debug_Object(V.face.X+V.face.Normal()*(T).003*ips.grid.dX.Min(),color_map[V.color_pair.x+3]);}
            else if(V.color_pair.x>=0) Add_Debug_Object(V.face.X-V.face.Normal()*(T).003*ips.grid.dX.Min(),color_map[V.color_pair.x+3]);}}
}

template<class T,class TV>
void Dump_System(const INTERFACE_POISSON_SYSTEM_COLOR<TV>& ips,ANALYTIC_POISSON_TEST<TV>& at)
{
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::CELL_ELEMENTS CELL_ELEMENTS;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::T_FACE T_FACE;
    typedef VECTOR<int,TV::m> TV_INT;

    Dump_Interface<T,TV>(ips);
    Flush_Frame("surfaces");

    for(CELL_ITERATOR<TV> it(ips.grid);it.Valid();it.Next()){
        int c=0;
        at.analytic_levelset->phi(it.Location(),0,c);
        Add_Debug_Particle(it.Location(),color_map[c+3]);}
    Flush_Frame("level set");
    
    char buff[100];
    for(int c=0;c<ips.cdi->colors;c++){
        Dump_Interface<T,TV>(ips);
        sprintf(buff,"dofs u%d",c);
        for(CELL_ITERATOR<TV> it(ips.grid);it.Valid();it.Next()){
            int index=ips.cm_u->Get_Index(it.index,c);
            if(index>=0){
                bool duplicated=false;
                for(int c_other=0;c_other<ips.cdi->colors;c_other++){
                    if(c_other==c) continue;
                    if(ips.cm_u->Get_Index(it.index,c_other)>=0) duplicated=true;}
                if(duplicated) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,1));
                else Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));}}                        
        Flush_Frame(buff);}

    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(ips.cdi->index_to_cell_elements);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<INTERFACE_ELEMENT>& interface_elements=cell_elements.interface;
        for(int i=0;i<interface_elements.m;i++){
            const INTERFACE_ELEMENT& V=interface_elements(i);
            Add_Debug_Particle(V.face.Center(),VECTOR<T,3>(0,0.1,0.5));
            T k=0;
            if(V.color_pair.x>=-1) k=at.j_surface(V.face.Center(),V.color_pair.x,V.color_pair.y);
            else if(V.color_pair.x==-2) k=at.u_jump(V.face.Center(),V.color_pair.x,V.color_pair.y);
            Debug_Particle_Set_Attribute<TV>("V",k*T_FACE::Normal(V.face.X));}}
    
    Dump_Interface<T,TV>(ips);
    for(CELL_ITERATOR<TV> it(ips.grid);it.Valid();it.Next()){
        int c=0;
        at.analytic_levelset->phi(it.Location(),0,c);
        if(c>=0){
            T f_volume=-at.mu(c)*at.analytic_solution(c)->L(it.Location(),0);
            Add_Debug_Particle(it.Location(),f_volume==0?VECTOR<T,3>(0.25,0.25,0.25):(f_volume>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0)));
            Debug_Particle_Set_Attribute<TV>("display_size",f_volume);}}
        Flush_Frame("volumetric forces");
}

template<class T,class TV>
void Dump_Vector(const INTERFACE_POISSON_SYSTEM_COLOR<TV>& ips,const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>& v,const char* title)
{
    char buff[100];
    for(int c=0;c<ips.cdi->colors;c++){
        Dump_Interface<T,TV>(ips);
        sprintf(buff,"%s u%d",title,c);
        for(CELL_ITERATOR<TV> it(ips.grid);it.Valid();it.Next()){
            int k=ips.cm_u->Get_Index(it.index,c);
            if(k>=0){
                T u_value=v.u(c)(k);
                Add_Debug_Particle(it.Location(),u_value==0?VECTOR<T,3>(0.25,0.25,0.25):(u_value>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0)));
                Debug_Particle_Set_Attribute<TV>("display_size",u_value);}}
        Flush_Frame(buff);}
}

template<class T,class TV>
void Dump_Vector(const INTERFACE_POISSON_SYSTEM_COLOR<TV>& ips,const ARRAY<T,VECTOR<int,TV::m> >& u,const char* title)
{
    INTERPOLATED_COLOR_MAP<T> cm;
    cm.Initialize_Colors(1e-12,1,true,true,true);

    Dump_Interface<T,TV>(ips);
    for(CELL_ITERATOR<TV> it(ips.grid);it.Valid();it.Next()){
        T u_value=u(it.index);
        Add_Debug_Particle(it.Location(),u_value==0?VECTOR<T,3>(0.25,0.25,0.25):(u_value>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0)));
        Debug_Particle_Set_Attribute<TV>("display_size",u_value);}
    Flush_Frame(title);
}

//#################################################################################################################################################
// Analytic Test ##################################################################################################################################
//#################################################################################################################################################

template<class TV>
void Analytic_Test(GRID<TV>& grid,ANALYTIC_POISSON_TEST<TV>& at,int max_iter,bool use_preconditioner,bool null,bool dump_matrix,bool debug_particles)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    ARRAY<ARRAY<T,TV_INT> > color_phi(at.analytic_solution.m+3);
    for(int i=0;i<color_phi.m;i++)
        color_phi(i).Resize(grid.Node_Indices(3));
    // ARRAY<T,TV_INT> phi_value(grid.Node_Indices());
    // ARRAY<int,TV_INT> phi_color(grid.Node_Indices());

    for(int c=0;c<color_phi.m;c++){
        for(NODE_ITERATOR<TV> it(grid,3);it.Valid();it.Next())
            color_phi(c)(it.index)=at.analytic_levelset->dist(it.Location(),0,c-3);

        LEVELSET<TV> ls(grid,color_phi(c),3);
        Reinitialize(ls,20,(T)0,(T)grid.dX.Max()*5,(T)100,(T).9,3,5,1);
        for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
            T p=color_phi(c)(it.index);
            Add_Debug_Particle(it.Location(),color_map[c]/(p>0?2:1));
            Debug_Particle_Set_Attribute<TV>("display_size",p);}
        Flush_Frame("reinitialized level set");
        LEVELSET<TV> levelset(grid,color_phi(c),1);
        for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
            Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
            Debug_Particle_Set_Attribute<TV>("V",levelset.Normal(it.Location()));}
        Flush_Frame("normals");}

    INTERFACE_POISSON_SYSTEM_COLOR<TV> ips(grid,color_phi);
    ips.use_preconditioner=use_preconditioner;
    ips.Set_Matrix(at.mu,true,
        [&](const TV& X,int color0,int color1){return at.u_jump(X,color0,color1);},
        [&](const TV& X,int color0,int color1){return at.j_surface(X,color0,color1);});
    printf("\n");
    for(int c=0;c<ips.cdi->colors;c++) printf("u%d [%i]\t",c,ips.cm_u->dofs(c));
    printf("\n");
    printf("q [%i] ",ips.cdi->constraint_base_scalar);
    printf("\n");

    INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV> rhs,sol;

    ips.Set_RHS(rhs,[&](const TV& X,int color){return -at.mu(color)*at.analytic_solution(color)->L(X,0);});
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
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        int c=0;
        at.analytic_levelset->phi(it.Location(),0,c);
        if(c>=0){
            int k=ips.cm_u->Get_Index(it.index,c);
            assert(k>=0);
            numer_u(it.index)=sol.u(c)(k);}}

    T avg_u=0;
    int cnt_u=0;
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        int c=0;
        at.analytic_levelset->phi(it.Location(),0,c);
        if(c>=0){
            exact_u(it.index)=at.analytic_solution(c)->f(it.Location(),0);
            error_u(it.index)=numer_u(it.index)-exact_u(it.index);
            avg_u+=error_u(it.index);
            cnt_u++;}}
    avg_u/=cnt_u;

    T error_u_linf=0,error_u_l2=0;
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        int c=0;
        at.analytic_levelset->phi(it.Location(),0,c);
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
        if(null){
            ARRAY<KRYLOV_VECTOR_BASE<T>*> array;
            ips.Compute_Nullspace(rhs,array,1);
            if(array.m){
                OCTAVE_OUTPUT<T>("n.txt").Write("n",*array(0));
                ips.Multiply(*array(0),*vectors(0));
                LOG::cout<<"Extra nullspace found: "<<sqrt(ips.Inner_Product(*vectors(0),*vectors(0)))<<std::endl;
                *array(0)*=1/static_cast<INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>*>(array(0))->Max_Abs();
                Dump_Vector<T,TV>(ips,*static_cast<INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>*>(array(0)),"extra null mode");}}}
    
    if(dump_matrix) OCTAVE_OUTPUT<T>("M.txt").Write("M",ips,*vectors(0));
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

    T m=1,s=1,kg=1;
    int threads=1;
    int test_number=1,resolution=4,max_iter=1000000;
    bool use_preconditioner=false,use_test=false,null=false,dump_matrix=false,debug_particles=false,opt_arg=false;
    bool test_analytic_diff=false;
    VIEWER_DIR viewer_dir("output");
    parse_args.Extra_Optional(&test_number,&opt_arg,"example number","example number to run");
    parse_args.Add("-o",&viewer_dir.output_directory,"output","output directory");
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
    parse_args.Add("-test_diff",&test_analytic_diff,"test analytic derivatives");
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

    ANALYTIC_POISSON_TEST<TV> test;

    auto Add_Exp_Sin_Cos=[=,&test](T s,T c,TV v,TV e)
        {
            auto f=Make_Analytic_Scalar<TV>([=](auto X,auto t){return (s*sin(X.Dot(v))+c*cos(X.Dot(v)))*exp(X.Dot(e));});
            test.analytic_solution.Append(f);
        };

    auto Add_Exp_Quad=[=,&test](T k,const MATRIX<T,TV::m>& M,TV a,T b)
        {
            auto f=Make_Analytic_Scalar<TV>([=](auto X,auto t){return k*exp(X.Dot(M*X+a)+b);});
            test.analytic_solution.Append(f);
        };

    switch(test_number){
        case 0:{ // One color, periodic. No interface, no forces, u=0.
            test.mu.Append(1);
            test.analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-1,0,-4);
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return 0;}));
            break;}
        case 1:{ // One color, periodic. No interface, u=sin(2*pi*x), f=(2*pi)^2*sin(2*pi*x).
            test.mu.Append(1);
            test.analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-1,0,-4);
            Add_Exp_Sin_Cos(1,0,TV::Axis_Vector(0)*(2*pi),TV());
            break;}
        case 2:{ // Two colors, periodic. Linear on [0,1/3],[1/3,2/3],[2/3,1], no volumetric forces.
            test.mu.Append(1);
            test.mu.Append(2);
            ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)/3,TV::Axis_Vector(0),0,1);
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-ANALYTIC_LEVELSET<TV>::Large_Phi(),0,-4);
            test.analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*((T)2/3),TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return (Value(X(0))>0.5)?(1-X(0)):(-X(0));}));
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return X(0)*2-1;}));
            break;}
        case 3:{ // Two colors, periodic. u=||x||^2 for r<R, zero elsewhere.
            test.mu.Append(1);
            test.mu.Append(2);
            test.analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T)0.5,1/(T)pi,1,0);
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return 0;}));
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return X.Magnitude_Squared()+X.Dot(TV()-1)+(T).25*TV::m;}));
            break;}
        case 4:{ // Two colors, periodic. u=exp(-x^2) for r<R, zero elsewhere.
            test.mu.Append(1);
            test.mu.Append(2);
            test.analytic_levelset=new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+(T)0.5,1/(T)pi,1,0);
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return 0;}));
            Add_Exp_Quad(1,MATRIX<T,TV::m>()-1,TV()+1,-(T).25*TV::m);
            break;}
        case 5:{ // Three colors, periodic. Stripes in x 0:[0,a], 1:[a,b], 2:[b,c], 0:[c,1].
            test.mu.Append(1);
            test.mu.Append(2);
            test.mu.Append(3);
            T a=(T)1/6,b=(T)5/12,c=(T)5/6;
            ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*a,TV::Axis_Vector(0),0,1);
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*c,TV::Axis_Vector(0),2,0);
            test.analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*b,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return 0;}));
            Add_Exp_Sin_Cos(0,1,TV(),TV::Axis_Vector(0));
            Add_Exp_Sin_Cos(1,0,TV::Axis_Vector(1)*(2*(T)pi),TV());
            break;}
        case 6:{ // Three colors (one dirichlet or neumann), periodic. Stripes in x 0:[0,a], 1:[a,b], 2:[b,c], 0:[c,1].
            T a=(T)1/6,b=(T)5/12,c=(T)5/6,constraint=-1;
            test.mu.Append(1);
            test.mu.Append(2);
            ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*a,TV::Axis_Vector(0),constraint,0);
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*c,TV::Axis_Vector(0),1,constraint);
            test.analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_LINE<TV>(TV::Axis_Vector(0)*b,TV::Axis_Vector(0),0,1)))->Add(ab)->Add(cd);
            Add_Exp_Sin_Cos(0,1,TV(),TV::Axis_Vector(0));
            Add_Exp_Sin_Cos(1,0,TV::Axis_Vector(1)*(2*(T)pi),TV());
            break;}
        case 7:{ // Three colors, periodic. u=a for r<R and x>0, u=b for r<R and x<0, zero elsewhere.
            T r=1/(T)pi,a0=0,a1=5,a2=7;
            TV n;
            for(int i=0;i<TV::m;i++) n(i)=i+(T)pi/(i+(T)pi);
            n.Normalize();
            test.mu.Append(1);
            test.mu.Append(2);
            test.mu.Append(3);
            ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV()+(T).5,n,2,1);
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-ANALYTIC_LEVELSET<TV>::Large_Phi(),0,-4);
            test.analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+.5,r,0,1)))->Add(ab)->Add(cd);
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return a0;}));
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return a1;}));
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return a2;}));
            break;}
        case 8:{ // Three colors, periodic. u=a*x^2 for r<R and x*n>0, u=b*x^2 for r<R and x*n<0, zero elsewhere.
            T r=1/(T)pi,a0=0,a1=5,a2=7;
            TV n;
            for(int i=0;i<TV::m;i++) n(i)=i+(T)pi/(i+(T)pi);
            n.Normalize();
            test.mu.Append(1);
            test.mu.Append(2);
            test.mu.Append(3);
            ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV()+(T).5,n,2,1);
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-ANALYTIC_LEVELSET<TV>::Large_Phi(),0,-4);
            test.analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+.5,r,0,1)))->Add(ab)->Add(cd);
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return a0*X.Magnitude_Squared()+X.Dot(TV()-a0)+(T).25*a0*TV::m;}));
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return a1*X.Magnitude_Squared()+X.Dot(TV()-a1)+(T).25*a1*TV::m;}));
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return a2*X.Magnitude_Squared()+X.Dot(TV()-a2)+(T).25*a2*TV::m;}));
            break;}
        case 9:{ // Three colors, periodic. u=a*exp(-x^2) for r<R and x*n>0, u=b*exp(-x^2) for r<R and x*n<0, zero elsewhere.
            T r=1/(T)pi,a0=0,a1=5,a2=7;
            TV n;
            for(int i=0;i<TV::m;i++) n(i)=i+(T)pi/(i+(T)pi);
            n.Normalize();
            test.mu.Append(1);
            test.mu.Append(2);
            test.mu.Append(3);
            ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV()+(T).5,n,2,1);
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-ANALYTIC_LEVELSET<TV>::Large_Phi(),0,-4);
            test.analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+.5,r,0,1)))->Add(ab)->Add(cd);
            Add_Exp_Quad(a0,MATRIX<T,TV::m>()-1,TV()+1,-(T).25*TV::m*1);
            Add_Exp_Quad(a1,MATRIX<T,TV::m>()-1,TV()+1,-(T).25*TV::m*1);
            Add_Exp_Quad(a2,MATRIX<T,TV::m>()-1,TV()+1,-(T).25*TV::m*1);
            break;}
        case 10:{ // Three colors (dirichlet/neumann outside), periodic. u=a*exp(-x^2) for r<R and x*n>0, u=b*exp(-x^2) for r<R and x*n<0, zero elsewhere.
            T r=1/(T)pi,a0=0,a1=5;
            TV n;
            int constraint=-2;
            for(int i=0;i<TV::m;i++) n(i)=i+(T)pi/(i+(T)pi);
            n.Normalize();
            test.mu.Append(1);
            test.mu.Append(2);
            ANALYTIC_LEVELSET_SIGNED<TV>* ab=new ANALYTIC_LEVELSET_LINE<TV>(TV()+(T).5,n,1,0);
            ANALYTIC_LEVELSET_SIGNED<TV>* cd=new ANALYTIC_LEVELSET_CONST<TV>(-ANALYTIC_LEVELSET<TV>::Large_Phi(),constraint,-4);
            test.analytic_levelset=(new ANALYTIC_LEVELSET_NEST<TV>(new ANALYTIC_LEVELSET_SPHERE<TV>(TV()+.5,r,0,1)))->Add(ab)->Add(cd);
            Add_Exp_Quad(a0,MATRIX<T,TV::m>()-1,TV()+1,-(T).25*TV::m*1);
            Add_Exp_Quad(a1,MATRIX<T,TV::m>()-1,TV()+1,-(T).25*TV::m*1);
            break;}
        case 11:{ // Four colors (dirichlet outside). Three equal bubbles.
            test.mu.Append(1);
            test.mu.Append(2);
            test.mu.Append(3);
            test.mu.Append(4);
            struct ANALYTIC_POISSON_LEVELSET_11:public ANALYTIC_LEVELSET<TV>
            {
                T r;
                TV n; // rotate angle with respect to e_y
                TV a; // shift
                VECTOR<TV,3> centers;
                VECTOR<TV,3> normals;
                ANALYTIC_POISSON_LEVELSET_11()
                {
                    r=(T)1/(2*pi-1);
                    n=TV::Axis_Vector(1);
                    a=TV()+pi/200;
                    centers(0)=TV::Axis_Vector(1);
                    centers(1).x=(T)sqrt(3)/2;centers(1).y=-(T)1/2;
                    centers(2).x=-(T)sqrt(3)/2;centers(2).y=-(T)1/2;
                    for(int i=0;i<3;i++){
                        normals(i).x=-centers(i).y;
                        normals(i).y=centers(i).x;
                        centers(i)*=r;}
                }
                T phi(const TV& X,T t,int& c) const
                {
                    TV x=X-0.5+a;
                    int i;
                    for(i=0;i<2;i++) if(x.Dot(normals((i+1)%3))>=0 && x.Dot(normals((i+2)%3))<0) break;
                    T d=(x-centers(i)).Magnitude();
                    if(d>r && x.Magnitude()>r/100){
                        c=3;
                        return d-r;}
                    c=i;
                    return min(abs(d-r),abs(x.Dot(normals((i+1)%3))),abs(x.Dot(normals((i+2)%3))));
                }
                TV N(const TV& X,T t,int c) const
                {
                    // Rough approximation
                    TV x=X-0.5+a;
                    if(c==3){
                        VECTOR<T,3> D((x-centers(0)).Magnitude()-r,(x-centers(1)).Magnitude()-r,(x-centers(2)).Magnitude()-r);
                        return -(x-centers(D.Arg_Min())).Normalized();}
                    if(c<0) return TV::Axis_Vector(0);
                    T x_dot_n0=-x.Dot(normals((c+1)%3)),x_dot_n1=x.Dot(normals((c+2)%3));
                    T d=(x-centers(c)).Magnitude()-r;
                    if(d>x_dot_n0 && d>x_dot_n1) return (x-centers(c)).Normalized();
                    if(x_dot_n0>x_dot_n1) return -normals((c+1)%3);
                    return normals((c+2)%3);
                }
                T dist(const TV& X,T t,int c) const
                {
                    // Rough approximation
                    TV x=X-0.5+a;
                    if(c==3){
                        VECTOR<T,3> D((x-centers(0)).Magnitude()-r,(x-centers(1)).Magnitude()-r,(x-centers(2)).Magnitude()-r);
                        return -D.Min();}
                    if(c<0) return ANALYTIC_LEVELSET<TV>::Large_Phi();
                    T x_dot_n0=-x.Dot(normals((c+1)%3)),x_dot_n1=x.Dot(normals((c+2)%3));
                    T d=(x-centers(c)).Magnitude()-r;
                    return max(max(x_dot_n0,x_dot_n1),d);
                }
            };
            test.analytic_levelset=new ANALYTIC_POISSON_LEVELSET_11;
            Add_Exp_Sin_Cos(1,0,TV::Axis_Vector(0),TV());
            Add_Exp_Sin_Cos(0,1,TV::Axis_Vector(1),TV());
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return X.Magnitude_Squared();}));
            test.analytic_solution.Append(Make_Analytic_Scalar<TV>([=](auto X,auto t){return 0;}));
            break;}
        default:{
            LOG::cerr<<"Unknown test number."<<std::endl; exit(-1); break;}}

    TV_INT counts=TV_INT()+resolution;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1),true);
    GRID<TV> grid2(counts*2,RANGE<TV>(TV(),TV()+1),true);

    VIEWER_OUTPUT vo(STREAM_TYPE((RW)0),viewer_dir);
    Use_Debug_Particles<TV>();
    Get_Debug_Particles<TV>().edge_separation=(T).02;

    if(test_analytic_diff) test.Test(grid.domain);

    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        int c=-4;
        T p=test.analytic_levelset->phi(it.Location(),0,c);
        Add_Debug_Particle(it.Location(),color_map[c+3]);
        Debug_Particle_Set_Attribute<TV>("display_size",p);}
    Flush_Frame("level set");
    for(int c=-3;c<test.analytic_solution.m;c++){
        for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
            T p=test.analytic_levelset->dist(it.Location(),0,c);
            Add_Debug_Particle(it.Location(),color_map[c+3]/(p>0?2:1));
            Debug_Particle_Set_Attribute<TV>("display_size",abs(p));
            Debug_Particle_Set_Attribute<TV>("V",test.analytic_levelset->N(it.Location(),0,c));}
        Flush_Frame("normals");}

    Analytic_Test(grid,test,max_iter,use_preconditioner,null,dump_matrix,debug_particles);
    Flush_Frame("finish");
    LOG::Finish_Logging();
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
