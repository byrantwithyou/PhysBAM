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
std::string output_directory="output";

template<class TV>
GRID<TV>* Global_Grid(GRID<TV>* grid_in=0)
{
    static GRID<TV>* grid=0;
    GRID<TV>* old_grid=grid;
    if(grid_in) grid=grid_in;
    return old_grid;
}

//#################################################################################################################################################
// Debug Particles ################################################################################################################################
//#################################################################################################################################################

template<class TV> DEBUG_PARTICLES<TV>& Get_Debug_Particles()
{
    static DEBUG_PARTICLES<TV> debug_particles;
    return debug_particles;
}

template<class T,class TV>
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
void Flush_Frame()
{
    Dump_Frame<T,TV>(ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >(*Global_Grid<TV>()),"flush");
}

template<class T,class TV>
void Dump_Frame(const VECTOR_ND<T>& v,INTERFACE_FLUID_SYSTEM<TV>& ifs,const char* title)
{
    INTERPOLATED_COLOR_MAP<T> color_map;
    color_map.Initialize_Colors(1e-11,10,true,true,true);

    char buff[100];
    for(int inside=0;inside<2;inside++)
    {
        sprintf(buff, title, inside?'-':'+');
        ARRAY<T,FACE_INDEX<TV::m> > error(ifs.grid);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(ifs.grid);it.Valid();it.Next()){
            int index=ifs.index_map_u[it.Axis()]->Get_Index_Fixed(it.index,inside);
            if(index>=0)
                error(it.Full_Index())=v(index);}

        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(ifs.grid);it.Valid();it.Next()){
            int index=ifs.index_map_p->Get_Index_Fixed(it.index,inside);
            if(index>=0)
                Add_Debug_Particle(ifs.grid.Center(it.index),color_map(v(index)));}

        for(int i=0;i<ifs.object.mesh.elements.m;i++){
            Add_Debug_Particle(ifs.object.Get_Element(i).Center(),VECTOR<T,3>(0,1,1));
            TV N;
            for(int j=0;j<TV::m;j++) N(j)=v(ifs.index_range_q[j].min_corner+i);
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,N);}

        Dump_Frame<T,TV>(error,buff);}
}

//#################################################################################################################################################
// Analytic Test ##################################################################################################################################
//#################################################################################################################################################

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
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(coarse_grid);it.Valid();it.Next())
        phi(it.index)=at.phi(it.Location());
    INTERFACE_FLUID_SYSTEM<TV> ifs(grid,coarse_grid,phi);
    ifs.Set_Matrix(at.mu);

    for(int i=0;i<ifs.object.mesh.elements.m;i++){
        Add_Debug_Particle(ifs.object.Get_Element(i).X(0),VECTOR<T,3>(1,1,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,ifs.object.Get_Element(i).X(1)-ifs.object.Get_Element(i).X(0));
        Add_Debug_Particle(ifs.object.Get_Element(i).X(1),VECTOR<T,3>(1,1,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,ifs.object.Get_Element(i).X(0)-ifs.object.Get_Element(i).X(1));
        Add_Debug_Particle(ifs.object.Get_Element(i).Center(),VECTOR<T,3>(0,1,0));}
    Flush_Frame<T,TV>();

    printf("\n");
    for(int i=0;i<TV::m;i++) printf("%c [%i %i) ", "uvw"[i], ifs.index_range_u[i].min_corner, ifs.index_range_u[i].max_corner);
    printf("p [%i %i) ", ifs.index_range_p.min_corner, ifs.index_range_p.max_corner);
    for(int i=0;i<TV::m;i++) printf("%cq [%i %i) ", "uvw"[i], ifs.index_range_q[i].min_corner, ifs.index_range_q[i].max_corner);
    printf("\n");

    ARRAY<TV,TV_INT> f_body[2];
    ARRAY<TV> f_interface;
    f_interface.Resize(ifs.object.mesh.elements.m);
    for(int s=0;s<2;s++) f_body[s].Resize(grid.Domain_Indices());

    KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T> > rhs,sol,kr_p,kr_ap,kr_ar,kr_r,kr_z;
    ifs.Resize_Vector(sol);
    ifs.Resize_Vector(kr_p);
    ifs.Resize_Vector(kr_ap);
    ifs.Resize_Vector(kr_ar);
    ifs.Resize_Vector(kr_r);
    ifs.Resize_Vector(kr_z);

    ARRAY<T,FACE_INDEX<TV::m> > exact_u,numer_u,error_u;
    ARRAY<T,TV_INT> exact_p,numer_p,error_p;

    for(int i=0; i<ifs.object.mesh.elements.m;i++)
        f_interface(i)=at.interface(ifs.object.Get_Element(i).Center());

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        for(int s=0;s<2;s++)
            f_body[s](it.index)=at.body(it.Location(),s);}

    ifs.Set_RHS(rhs,f_body,f_interface);

    for(int i=0;i<TV::m;i++){
        char buff[100];
        sprintf(buff, "null %c %%c", "uvw"[i]);
        Dump_Frame(ifs.null_u[i],ifs,buff);
    }
    Dump_Frame(ifs.null_p,ifs,"null p %c");

    CONJUGATE_RESIDUAL<T> cr;
    cr.print_residuals=true;
    cr.Solve(ifs,sol,rhs,kr_p,kr_ap,kr_ar,kr_r,kr_z,0,0,100000);

    ifs.Multiply(sol,kr_r);
    kr_r.v-=rhs.v;
    LOG::cout<<"Residual: "<<ifs.Convergence_Norm(kr_r)<<std::endl;

    OCTAVE_OUTPUT<T>("M.txt").Write("M",ifs,kr_r,kr_z);
    OCTAVE_OUTPUT<T>("b.txt").Write("b",rhs);
    OCTAVE_OUTPUT<T>("x.txt").Write("x",sol);
    OCTAVE_OUTPUT<T>("r.txt").Write("r",kr_r);

    ifs.Get_U_Part(sol.v,numer_u);
    ifs.Get_P_Part(sol.v,numer_p);

    exact_u.Resize(grid);
    error_u.Resize(grid);
    exact_p.Resize(grid.Domain_Indices());
    error_p.Resize(grid.Domain_Indices());

    TV avg_u;
    TV_INT cnt_u;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face(it.Full_Index()); 
        if (abs(numer_u(face))<1e-10) numer_u(face)=0;
        exact_u(face)=at.u(it.Location())(face.axis);
        error_u(face)=numer_u(face)-exact_u(face);
        avg_u(face.axis)+=error_u(face);
        cnt_u(face.axis)++;}
    avg_u/=(TV)cnt_u;

    TV error_u_linf,error_u_l2;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        FACE_INDEX<TV::m> face(it.Full_Index()); 
        T d=error_u(face)-avg_u(face.axis);
        error_u_linf(face.axis)=max(error_u_linf(face.axis),abs(d));
        error_u_l2(face.axis)+=sqr(d);}
    error_u_l2=sqrt(error_u_l2/(TV)cnt_u);

    LOG::cout<<"U error:   linf="<<error_u_linf<<"   l2="<<error_u_l2<<std::endl;

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
        T d=error_p(it.index)-avg_p;
        error_p_linf=max(error_p_linf,abs(d));
        error_p_l2+=sqr(d);}
    error_p_l2=sqrt(error_p_l2/cnt_p);

    LOG::cout<<"P error:   linf "<<error_p_linf<<"   l2 "<<error_p_l2<<std::endl<<std::endl;
    
    LOG::cout<<"exact u"<<std::endl<<std::endl<<exact_u<<std::endl;
    LOG::cout<<"numer u"<<std::endl<<std::endl<<numer_u<<std::endl;
    LOG::cout<<"error u"<<std::endl<<std::endl<<error_u<<std::endl;
}

//#################################################################################################################################################
// Integration Test ###############################################################################################################################
//#################################################################################################################################################

template<class TV>
void Integration_Test(int argc,char* argv[])
{
    typedef typename TV::SCALAR T;
    static const int d=TV::m;
    typedef VECTOR<int,d> TV_INT;

    Get_Debug_Particles<TV>();

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
                {return TV::Axis_Vector(1)*((phi(X)<0)?(2*X.x-m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X){return -m/(T)6+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV();}
                virtual TV interface(const TV& X)
                {return TV::Axis_Vector(1)*((X.x>0.5*m)?1:-1)*(2*mu(1)+mu(0))/s;}
            };
            test=new ANALYTIC_TEST_0;
            break;}
        case 1:{
            struct ANALYTIC_TEST_1:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual TV u(const TV& X)
                {return TV::Axis_Vector(1)*((phi(X)<0)?(X.x-0.5*m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X){return -0.25*m+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV();}
                virtual TV interface(const TV& X)
                {return TV::Axis_Vector(1)*((X.x>0.5*m)?1:-1)*mu.Sum()/s;}
            };
            test=new ANALYTIC_TEST_1;
            break;}
        case 2:{
            struct ANALYTIC_TEST_2:public ANALYTIC_TEST<TV>
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
            test=new ANALYTIC_TEST_2;
            break;}
        case 3:{
            struct ANALYTIC_TEST_3:public ANALYTIC_TEST<TV>
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
            test=new ANALYTIC_TEST_3;
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

    Global_Grid(&grid);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    Analytic_Test(grid,coarse_grid,*test);
}

//#################################################################################################################################################
// Main ###########################################################################################################################################
//#################################################################################################################################################

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    Integration_Test<VECTOR<double,2> >(argc,argv);

    return 0;
}
