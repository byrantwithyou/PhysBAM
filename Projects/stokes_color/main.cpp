//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/MINRES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Finite_Elements/ANALYTIC_BOUNDARY_CONDITIONS_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

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

template<class TV> struct ANALYTIC_TEST;

// TODO - Insert dump debug particles routines

//#################################################################################################################################################
// Analytic Test ##################################################################################################################################
//#################################################################################################################################################

template<class TV>
struct ANALYTIC_TEST: public ANALYTIC_BOUNDARY_CONDITIONS_COLOR<TV>
{
    typedef typename TV::SCALAR T;
    using ANALYTIC_BOUNDARY_CONDITIONS_COLOR<TV>::kg;
    using ANALYTIC_BOUNDARY_CONDITIONS_COLOR<TV>::m;
    using ANALYTIC_BOUNDARY_CONDITIONS_COLOR<TV>::s;

    bool wrap;
    ARRAY<T> mu;

    virtual void Initialize()=0;
    virtual T phi_value(const TV& X)=0;
    virtual int phi_color(const TV& X)=0;
    virtual TV u(const TV& X,int color)=0;
    virtual T p(const TV& X)=0;
    virtual TV f_volume(const TV& X,int color)=0;

    TV u(const TV& X){return u(X,phi_color(X));}
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
    

    INTERFACE_STOKES_SYSTEM_COLOR<TV> iss(grid,phi_value,phi_color);
    iss.use_preconditioner=use_preconditioner;
    iss.Set_Matrix(at.mu,at.wrap,&at);

    printf("\n");
    for(int i=0;i<TV::m;i++){for(int c=0;c<iss.cdi->colors;c++) printf("%c%d [%i]\t","uvw"[i],c,iss.cm_u(i)->dofs(c));printf("\n");}
    for(int c=0;c<iss.cdi->colors;c++) printf("p%d [%i]\t",c,iss.cm_p->dofs(c));printf("\n");
    printf("qn [%i]\t",iss.cdi->constraint_base_normal);
    printf("qt [%i] ",iss.cdi->constraint_base_tangent);
    printf("\n");

    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> rhs,sol;
    {
        ARRAY<ARRAY<TV,TV_INT> > f_volume;
        ARRAY<ARRAY<T,FACE_INDEX<TV::m> > > u;
        
        f_volume.Resize(iss.cdi->colors);
        u.Resize(iss.cdi->colors);
        
        for(int c=0;c<iss.cdi->colors;c++) f_volume(c).Resize(grid.Domain_Indices());
        for(int c=0;c<iss.cdi->colors;c++) u(c).Resize(grid);
        
        for(int c=0;c<iss.cdi->colors;c++)
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
                f_volume(c)(it.index)=at.f_volume(it.Location(),c);
        
        for(int c=0;c<iss.cdi->colors;c++)
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
                FACE_INDEX<TV::m> face(it.Full_Index()); 
                u(c)(face)=at.u(it.Location(),c)(face.axis);}
        
        iss.Set_RHS(rhs,f_volume,u);
        iss.Resize_Vector(sol);
    }

    MINRES<T> mr;
    KRYLOV_SOLVER<T>* solver=&mr;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    solver->Solve(iss,sol,rhs,vectors,1e-10,0,max_iter);
    
    iss.Multiply(sol,*vectors(0));
    *vectors(0)-=rhs;
    LOG::cout<<"Residual: "<<iss.Convergence_Norm(*vectors(0))<<std::endl;

    for(int i=0;i<TV::m;i++){
        iss.Multiply(iss.null_u(i),*vectors(0));
        LOG::cout<<"null u["<<i<<"] "<<iss.Convergence_Norm(*vectors(0))<<std::endl;}
    iss.Multiply(iss.null_p,*vectors(0));
    LOG::cout<<"null p "<<" "<<iss.Convergence_Norm(*vectors(0))<<std::endl;

    ARRAY<T,FACE_INDEX<TV::m> > exact_u,numer_u,error_u;
    ARRAY<T,TV_INT> exact_p,numer_p,error_p;

    numer_u.Resize(iss.grid);
    exact_u.Resize(iss.grid);
    error_u.Resize(iss.grid);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        int i=it.Axis();
        int c=at.phi_color(it.Location());
        int k=iss.cm_u(it.Axis())->Get_Index(it.index,c);
        assert(k>=0);
        numer_u(it.Full_Index())=sol.u(i)(c)(k);}

    numer_p.Resize(iss.grid.Domain_Indices());
    exact_p.Resize(iss.grid.Domain_Indices());
    error_p.Resize(iss.grid.Domain_Indices());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        int k=iss.cm_p->Get_Index(it.index,c);
        assert(k>=0);
        numer_p(it.index)=sol.p(c)(k);}

    TV avg_u;
    TV_INT cnt_u;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face(it.Full_Index()); 
        exact_u(face)=at.u(it.Location())(face.axis);
        error_u(face)=numer_u(face)-exact_u(face);
        avg_u(face.axis)+=error_u(face);
        cnt_u(face.axis)++;}
    avg_u/=(TV)cnt_u;

    TV error_u_linf,error_u_l2;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face(it.Full_Index()); 
        error_u(face)-=avg_u(face.axis);
        error_u_linf(face.axis)=max(error_u_linf(face.axis),abs(error_u(face)));
        error_u_l2(face.axis)+=sqr(error_u(face));}
    error_u_l2=sqrt(error_u_l2/(TV)cnt_u);

    LOG::cout<<iss.grid.counts<<" U error:   linf="<<error_u_linf<<"   l2="<<error_u_l2<<std::endl;

    T avg_p=0;
    int cnt_p=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        exact_p(it.index)=at.p(it.Location());
        error_p(it.index)=numer_p(it.index)-exact_p(it.index);
        avg_p+=error_p(it.index);
        cnt_p++;}
    avg_p/=cnt_p;

    T error_p_linf=0,error_p_l2=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        error_p(it.index)-=avg_p;
        error_p_linf=max(error_p_linf,abs(error_p(it.index)));
        error_p_l2+=sqr(error_p(it.index));}
    error_p_l2=sqrt(error_p_l2/cnt_p);

    LOG::cout<<iss.grid.counts<<" P error:   linf "<<error_p_linf<<"   l2 "<<error_p_l2<<std::endl<<std::endl;

    /*if(debug_particles){
        Dump_System<T,TV>(iss,at);
        Dump_Vector<T,TV>(iss,sol,"solution");
        Dump_u_p(iss,error_u,error_p,"error");
        Dump_u_p(iss.grid,error_u,error_p,"color mapped error");
        if(null&&iss.Nullspace_Check(rhs)){
            OCTAVE_OUTPUT<T>("n.txt").Write("n",rhs);
            iss.Multiply(rhs,*vectors(0));
            LOG::cout<<"nullspace found: "<<sqrt(iss.Inner_Product(*vectors(0),*vectors(0)))<<std::endl;
            rhs*=1/rhs.Max_Abs();
            Dump_Vector2<T,TV>(iss,rhs,"extra null mode");}}*/
        
    if(dump_matrix) OCTAVE_OUTPUT<T>("M.txt").Write("M",iss,*vectors(0),*vectors(1));
}

//#################################################################################################################################################
// Integration Test ###############################################################################################################################
//#################################################################################################################################################

template<class TV>
void Integration_Test(int argc,char* argv[],PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    // Get_Debug_Particles<TV>().debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);

    int test_number;
    if(parse_args.Num_Extra_Args()<1){LOG::cerr<<"Test number is required."<<std::endl; exit(-1);}
    if(!STRING_UTILITIES::String_To_Value(parse_args.Extra_Arg(0),test_number)) throw VALUE_ERROR("The argument is not an integer.");

    ANALYTIC_TEST<TV>* test=0;

    switch(test_number){
        case 0:{
            struct ANALYTIC_TEST_0:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);}
                virtual T phi_value(const TV& X){return abs(-m/(T)6+abs(X.x-0.5*m));}
                virtual int phi_color(const TV& X){return (-m/(T)6+abs(X.x-0.5*m))<0;}
                virtual TV u(const TV& X,int color){return TV::Axis_Vector(1)*(color?(2*X.x-m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return T();}
                virtual TV f_volume(const TV& X,int color){return TV();}
                virtual TV f_surface(const TV& X,int color1,int color2){return TV::Axis_Vector(1)*((X.x>0.5*m)?(T)(-1):(T)1)*(2*mu(1)+mu(0))/s;}
            };
            test=new ANALYTIC_TEST_0;
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
