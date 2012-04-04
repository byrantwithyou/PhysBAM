//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_INTEGRATION_CUTTING.h>
#include <PhysBAM_Geometry/Finite_Elements/BASIS_STENCIL_UNIFORM.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MAPPING.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_FLUID_SYSTEM.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

using namespace PhysBAM;

typedef float RW;
template<class TV>
struct ANALYTIC_TEST
{
    typedef typename TV::SCALAR T;
    T kg,m,s;
    VECTOR<T,2> mu;

    virtual TV u(const TV& X)=0;
    virtual T p(const TV& X)=0;
    virtual T phi(const TV& X)=0;
    virtual TV body(const TV& X,bool inside)=0;
    virtual TV interface(const TV& X)=0;
};

template<class TV>
void Analytic_Test(GRID<TV>& grid,GRID<TV>& coarse_grid,ANALYTIC_TEST<TV>& at)
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;

    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next())
        phi(it.index)=at.phi(it.Location());
    INTERFACE_FLUID_SYSTEM<TV> ifs(grid,coarse_grid,phi);
    ifs.Set_Matrix(at.mu);

    printf("\n");
    for(int i=0;i<TV::m;i++) printf("%c [%i %i) ", "uvw"[i], ifs.index_range_u[i].min_corner, ifs.index_range_u[i].max_corner);
    printf("p [%i %i) ", ifs.index_range_p.min_corner, ifs.index_range_p.max_corner);
    for(int i=0;i<TV::m;i++) printf("%cq [%i %i) ", "uvw"[i], ifs.index_range_q[i].min_corner, ifs.index_range_q[i].max_corner);
    printf("\n");

    ARRAY<TV,TV_INT> f_body[2];
    ARRAY<TV> f_interface;
    f_interface.Resize(ifs.object.mesh.elements.m);
    for(int s=0;s<2;s++) f_body[s].Resize(grid.Domain_Indices());

    KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > rhs,sol,kr_q,kr_s,kr_t,kr_r;
    ifs.Resize_Vector(sol);
    ifs.Resize_Vector(kr_q);
    ifs.Resize_Vector(kr_s);
    ifs.Resize_Vector(kr_t);
    ifs.Resize_Vector(kr_r);

    ARRAY<T,FACE_INDEX<TV::m> > exact_u,numer_u;
    ARRAY<T,TV_INT> exact_p,numer_p;

    for(int i=0; i<ifs.object.mesh.elements.m;i++)
        f_interface(i)=at.interface(ifs.object.Get_Element(i).Center());

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        for(int s=0;s<2;s++)
            f_body[s](it.index)=at.body(it.Location(),s);}

    ifs.Set_RHS(rhs,f_body,f_interface);

    CONJUGATE_RESIDUAL<T> cr;
    cr.Solve(ifs,sol,rhs,kr_q,kr_s,kr_t,kr_r,1e-10,0,100000);

    OCTAVE_OUTPUT<T>("M.txt").Write("M",ifs,kr_q,kr_s);
    OCTAVE_OUTPUT<T>("b.txt").Write("b",rhs);
    OCTAVE_OUTPUT<T>("x.txt").Write("x",sol);

    ifs.Get_U_Part(sol.v,numer_u);
    ifs.Get_P_Part(sol.v,numer_p);

    T u_linf=0,u_l2=0;
    int num_u=0;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        T a=numer_u(it.Full_Index()),b=at.u(it.Location())(it.Axis()),d=fabs(a-b);
        num_u++;
        u_linf=max(u_linf,d);
        u_l2+=sqr(d);}
    u_l2=sqrt(u_l2/num_u);

    T p_linf=0,p_l2=0;
    int num_p=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        T a=numer_p(it.index),b=at.p(it.Location()),d=fabs(a-b);
        num_p++;
        p_linf=max(p_linf,d);
        p_l2+=sqr(d);}
    p_l2=sqrt(p_l2/num_p);

    LOG::cout<<"U: linf "<<u_linf<<"  l2 "<<u_l2<<std::endl;
    LOG::cout<<"P: linf "<<p_linf<<"  l2 "<<p_l2<<std::endl;
}

template<class TV>
void Integration_Test(int argc,char* argv[])
{
    typedef typename TV::SCALAR T;
    static const int d=TV::m;
    typedef VECTOR<int,d> TV_INT;

    PARSE_ARGS parse_args;
    parse_args.Set_Extra_Arguments(-1,"<example number>");
    parse_args.Add_Double_Argument("-mu_i",1,"viscosity inside");
    parse_args.Add_Double_Argument("-mu_o",1,"viscosity outside");
    parse_args.Add_Double_Argument("-m",1,"meter scale");
    parse_args.Add_Double_Argument("-sec",1,"second scale");
    parse_args.Add_Double_Argument("-kg",1,"kilogram scale");
    parse_args.Add_Integer_Argument("-test",1,"test number");
    parse_args.Add_Integer_Argument("-resolution",4,"resolution");
    parse_args.Add_Integer_Argument("-cgf",2,"coarse grid factor");
    parse_args.Parse(argc,argv);

    int test_number;
    if(parse_args.Num_Extra_Args()<1){LOG::cerr<<"Test number is required."<<std::endl; exit(-1);}
    if(!STRING_UTILITIES::String_To_Value(parse_args.Extra_Arg(0),test_number)) throw VALUE_ERROR("The argument is not an integer.");

    ANALYTIC_TEST<TV>* test=0;

    switch(test_number){
        case 0:{
            struct ANALYTIC_TEST_0:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual TV u(const TV& X)
                {return TV::Axis_Vector(1)*(phi(X)<0?X.x-0.5*m:(X.x>0.5*m?m-X.x:-X.x))/s;}
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X){return -0.25*m+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV();}
                virtual TV interface(const TV& X)
                {return TV::Axis_Vector(1)*sign(X.x>0.5*m)*mu.Sum()/s;}
            };
            test=new ANALYTIC_TEST_0;
            break;}
        case 1:{
            struct ANALYTIC_TEST_1:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual TV u(const TV& X)
                {return TV::Axis_Vector(1)*(phi(X)<0)*(sqr(0.25*m)-sqr(X.x-0.5*m))/(m*s);}
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X){return -0.25*m+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV::Axis_Vector(1)*(inside*2*mu(1)/(m*s));}
                virtual TV interface(const TV& X)
                {return -(T)0.5*TV::Axis_Vector(1)*mu(1)/s;}
            };
            test=new ANALYTIC_TEST_1;
            break;}
        case 2:{
            struct ANALYTIC_TEST_2:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual TV u(const TV& X)
                {T x=abs(X.x-(T)0.5*m)-(T)0.25*m;return -TV::Axis_Vector(1)*x*(0.5*m-abs(x))/(m*s);}
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X){return -0.25*m+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV::Axis_Vector(1)*(inside?mu(1):-mu(0))*2/(m*s);}
                virtual TV interface(const TV& X)
                {return (T)0.5*TV::Axis_Vector(1)*(mu(0)-mu(1))/s;}
            };
            test=new ANALYTIC_TEST_2;
            break;}
        default:{
        LOG::cerr<<"Unknown test number."<<std::endl; exit(-1); break;}}

    test->m=parse_args.Get_Double_Value("-m");
    test->s=parse_args.Get_Double_Value("-sec");
    test->kg=parse_args.Get_Double_Value("-kg");
    test->mu(0)=parse_args.Get_Double_Value("-mu_o")*test->kg/(test->s*(d==3?test->m:1));
    test->mu(1)=parse_args.Get_Double_Value("-mu_i")*test->kg/(test->s*(d==3?test->m:1));
    int res=parse_args.Get_Integer_Value("-resolution");
    int cgf=parse_args.Get_Integer_Value("-cgf");

    if(res%cgf) PHYSBAM_FATAL_ERROR("Resolution must be divisible by coarse grid factor.");

    TV_INT counts=TV_INT()+res;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1)*test->m,true);
    GRID<TV> coarse_grid(grid.counts/cgf,grid.domain,true);
    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());

    std::string output_directory="output";
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    Analytic_Test(grid,coarse_grid,*test);
}

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    Integration_Test<VECTOR<double,2> >(argc,argv);

    return 0;
}
