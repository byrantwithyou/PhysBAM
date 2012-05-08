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
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM.h>
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

//#################################################################################################################################################
// Debug Particles ################################################################################################################################
//#################################################################################################################################################

template<class TV> DEBUG_PARTICLES<TV>& Get_Debug_Particles()
{
    static DEBUG_PARTICLES<TV> debug_particles;
    return debug_particles;
}

template<class T,int d>
void Dump_Frame(const ARRAY<T,FACE_INDEX<d> >& u,const char* title)
{
    static int frame=0;
    char buff[100];
    sprintf(buff, "%s/%i", output_directory.c_str(), frame);
    FILE_UTILITIES::Create_Directory(buff);
    FILE_UTILITIES::Write_To_File<RW>((std::string)buff+"/mac_velocities.gz",u);
    if(title) FILE_UTILITIES::Write_To_Text_File((std::string)buff+"/frame_title",title);
    Get_Debug_Particles<VECTOR<T,d> >().Write_Debug_Particles(STREAM_TYPE((RW())),output_directory,frame);
    frame++;
}

template<class T,class TV>
void Flush_Frame(const char* title)
{
    Dump_Frame(ARRAY<typename TV::SCALAR,FACE_INDEX<TV::m> >(*Global_Grid<TV>()),title);
}

template<class T,class TV>
void Dump_Interface(const INTERFACE_STOKES_SYSTEM<TV>& iss,const bool arrows)
{
    for(int i=0;i<iss.object.mesh.elements.m;i++){
        Add_Debug_Particle(iss.object.Get_Element(i).X(0),VECTOR<T,3>(0,0.5,0));
        if(arrows) Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,(iss.object.Get_Element(i).X(1)-iss.object.Get_Element(i).X(0))/iss.grid.dX);
        Add_Debug_Particle(iss.object.Get_Element(i).X(1),VECTOR<T,3>(0,0.5,0));
        if(arrows) Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,(iss.object.Get_Element(i).X(0)-iss.object.Get_Element(i).X(1))/iss.grid.dX);
        Add_Debug_Particle(iss.object.Get_Element(i).Center(),VECTOR<T,3>(0,0.3,1));}
}

template<class T,class TV>
void Dump_System(const INTERFACE_STOKES_SYSTEM<TV>& iss,ANALYTIC_TEST<TV>& at)
{
    Dump_Interface<T,TV>(iss,true);
    Flush_Frame<T,TV>("interface");
    
    char buff[100];
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            Dump_Interface<T,TV>(iss,false);
            sprintf(buff,"dofs %c %c","uvw"[i],s?'-':'+');
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(iss.grid);it.Valid();it.Next())
                if(it.Axis()==i){
                    int index=iss.cm_u(i)->Get_Index(it.index,s);
                    if(index>=0){
                        if(iss.active_dofs.u(i)[s](index)==0)
                            Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
                        else if(iss.cm_u(i)->Get_Index(it.index,1-s)>=0)
                            Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));
                        else
                            Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0.5,0));
                        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV::Axis_Vector(i));}}                        
            Flush_Frame<T,TV>(buff);
        }

    for(int s=0;s<2;s++){
        Dump_Interface<T,TV>(iss,false);
        sprintf(buff,"dofs p %c",s?'-':'+');
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(iss.grid);it.Valid();it.Next()){
            int index=iss.cm_p->Get_Index(it.index,s);
            if(index>=0){
                if(iss.active_dofs.p[s](index)==0)
                    Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
                else if(iss.cm_p->Get_Index(it.index,1-s)>=0)
                    Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));
                else
                    Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0.5,0));}}
        Flush_Frame<T,TV>(buff);}

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(iss.grid);it.Valid();it.Next()){
        if(at.phi(it.Location())<0) Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,0.3,1));
        else Add_Debug_Particle(it.Location(),VECTOR<T,3>(0.9,0.2,0.2));}
    Flush_Frame<T,TV>("analytic level set");

    for(int i=0;i<iss.object.mesh.elements.m;i++){
        Add_Debug_Particle(iss.object.Get_Element(i).Center(),VECTOR<T,3>(0,0.1,0.5));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,at.interface(iss.object.Get_Element(i).Center()));}
    Flush_Frame<T,TV>("analytic interfacial forces");
    
    Dump_Interface<T,TV>(iss,false);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(iss.grid);it.Valid();it.Next()){
        int s=at.phi(it.Location())<0;
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(0.25,0.25,0.25));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,at.body(it.Location(),s));}
    Flush_Frame<T,TV>("analytic volumetric forces");
}

template<class T,class TV>
void Dump_Vector(const INTERFACE_STOKES_SYSTEM<TV>& iss,const INTERFACE_STOKES_SYSTEM_VECTOR<TV>& v,const char* title)
{
    char buff[100];
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            Dump_Interface<T,TV>(iss,false);
            sprintf(buff,"%s %c %c",title,"uvw"[i],s?'-':'+');
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(iss.grid);it.Valid();it.Next())
                if(it.Axis()==i){ 
                    int k=iss.cm_u(i)->Get_Index(it.index,s);
                    if(k>=0){
                        if(iss.cm_u(i)->Get_Index(it.index,1-s)>=0) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));
                        else Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0.5,0));
                        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV::Axis_Vector(i)*v.u(i)[s](k));}}
            Flush_Frame<T,TV>(buff);}

    for(int s=0;s<2;s++){
        Dump_Interface<T,TV>(iss,false);
        sprintf(buff,"%s p %c",title,s?'-':'+');
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(iss.grid);it.Valid();it.Next()){
            int k=iss.cm_p->Get_Index(it.index,s);
            if(k>=0){
                Add_Debug_Particle(it.Location(),v.p[s](k)>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,v.p[s](k));}}
        Flush_Frame<T,TV>(buff);}
}

template<class T,class TV>
void Dump_Vector2(const INTERFACE_STOKES_SYSTEM<TV>& iss,const INTERFACE_STOKES_SYSTEM_VECTOR<TV>& v,const char* title)
{
    char buff[100];
    for(int s=0;s<2;s++){
        ARRAY<T,FACE_INDEX<TV::m> > u(iss.grid);
        sprintf(buff,"%s %c",title,s?'-':'+');
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(iss.grid);it.Valid();it.Next()){
            int i=it.Axis();
            int k=iss.cm_u(i)->Get_Index(it.index,s);
            if(k>=0) u(it.Full_Index())=v.u(i)[s](k);}
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(iss.grid);it.Valid();it.Next()){
            int k=iss.cm_p->Get_Index(it.index,s);
            if(k>=0){
                Add_Debug_Particle(it.Location(),v.p[s](k)>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,v.p[s](k));}}
        for(int i=0;i<iss.object.mesh.elements.m;i++){
            Add_Debug_Particle(iss.object.Get_Element(i).Center(),VECTOR<T,3>(1,1,1));
            TV V;
            for(int j=0;j<TV::m;j++) V(j)=v.q(j)(i);
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,V);}
        Dump_Frame(u,buff);}
}

template<class T,class TV,class TV_INT>
void Dump_u_p(const INTERFACE_STOKES_SYSTEM<TV>& iss,const ARRAY<T,FACE_INDEX<TV::m> >& u,const ARRAY<T,TV_INT>& p,const char* title)
{
    char buff[100];
    for(int i=0;i<TV::m;i++){
        Dump_Interface<T,TV>(iss,false);
        sprintf(buff,"%s %c",title,"uvw"[i]);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(iss.grid);it.Valid();it.Next())
            if(it.Axis()==i){ 
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(0.25,0.25,0.25));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV::Axis_Vector(i)*u(it.Full_Index()));}
        Flush_Frame<T,TV>(buff);}

    Dump_Interface<T,TV>(iss,false);
    sprintf(buff,"%s p",title);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(iss.grid);it.Valid();it.Next()){
        Add_Debug_Particle(it.Location(),p(it.index)>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,p(it.index));}
    Flush_Frame<T,TV>(buff);
}

template<class T,class TV,class TV_INT>
void Dump_u_p(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& u,const ARRAY<T,TV_INT>& p,const char* title)
{
    INTERPOLATED_COLOR_MAP<T> color_map;
    color_map.Initialize_Colors(1e-12,1,true,true,true);

    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        Add_Debug_Particle(it.Location(),color_map(p(it.index)));

    Dump_Frame(u,title);
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

    virtual void Initialize()=0;
    virtual TV u(const TV& X,bool inside)=0;
    virtual T p(const TV& X)=0;
    virtual T phi(const TV& X)=0;
    virtual TV body(const TV& X,bool inside)=0;
    virtual TV interface(const TV& X)=0;

    TV u(const TV& X){return u(X,phi(X)<0);}
};

template<class TV>
void Analytic_Test(GRID<TV>& grid,GRID<TV>& coarse_grid,ANALYTIC_TEST<TV>& at,int max_iter,bool use_preconditioner,bool null,bool dump_matrix,bool debug_particles)
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;

    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(coarse_grid);it.Valid();it.Next()){
        T value=at.phi(it.Location());
        T tol=grid.dX.Min()*1e-6;
        if(abs(value)>tol) phi(it.index)=value;
        else phi(it.index)=((value>0)?(tol):(-tol));}
    INTERFACE_STOKES_SYSTEM<TV> iss(grid,coarse_grid,phi);
    iss.use_preconditioner=use_preconditioner;
    iss.Set_Matrix(at.mu);

    printf("\n");
    for(int i=0;i<TV::m;i++) for(int s=0;s<2;s++) printf("%c%c [%i] ","uvw"[i],"+-"[s],iss.cm_u(i)->dofs[s]);
    for(int s=0;s<2;s++) printf("p%c [%i] ","+-"[s],iss.cm_p->dofs[s]);
    printf("q [%i] ",iss.object.mesh.elements.m);
    printf("\n");

    INTERFACE_STOKES_SYSTEM_VECTOR<TV> rhs,sol;
    {
        VECTOR<ARRAY<TV,TV_INT>,2> f_body;
        ARRAY<TV> f_interface;
        VECTOR<ARRAY<T,FACE_INDEX<TV::m> >,2> u;
        
        f_interface.Resize(iss.object.mesh.elements.m);
        for(int s=0;s<2;s++) f_body[s].Resize(grid.Domain_Indices());
        for(int s=0;s<2;s++) u[s].Resize(grid);
        
        for(int i=0; i<iss.object.mesh.elements.m;i++)
            f_interface(i)=at.interface(iss.object.Get_Element(i).Center());
        
        for(int s=0;s<2;s++)
            for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
                f_body[s](it.index)=at.body(it.Location(),s);
        
        for(int s=0;s<2;s++)
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
                FACE_INDEX<TV::m> face(it.Full_Index()); 
                u[s](face)=at.u(it.Location(),s)(face.axis);}
        
        iss.Set_RHS(rhs,f_body,f_interface,u);
        iss.Resize_Vector(sol);
    }

    MINRES<T> mr;
    KRYLOV_SOLVER<T>* solver=&mr;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;

    solver->Solve(iss,sol,rhs,vectors,1e-10,0,max_iter);
    
    iss.Multiply(sol,*vectors(0));
    *vectors(0)-=rhs;
    LOG::cout<<"Residual: "<<iss.Convergence_Norm(*vectors(0))<<std::endl;

    ARRAY<T,FACE_INDEX<TV::m> > exact_u,numer_u,error_u;
    ARRAY<T,TV_INT> exact_p,numer_p,error_p;

    numer_u.Resize(iss.grid);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        int i=it.Axis();
        int s=at.phi(it.Location())<0;
        int k=iss.cm_u(it.Axis())->Get_Index(it.index,s);
        assert(k>=0);
        numer_u(it.Full_Index())=sol.u(i)[s](k);}

    numer_p.Resize(iss.grid.Domain_Indices());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int s=at.phi(it.Location())<0;
        int k=iss.cm_p->Get_Index(it.index,s);
        assert(k>=0);
        numer_p(it.index)=sol.p[s](k);}

    exact_u.Resize(grid);
    error_u.Resize(grid);
    exact_p.Resize(grid.Domain_Indices());
    error_p.Resize(grid.Domain_Indices());

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

    if(debug_particles){
        Dump_System<T,TV>(iss,at);
        Dump_Vector<T,TV>(iss,sol,"solution");
        Dump_u_p(iss,error_u,error_p,"error");
        Dump_u_p(iss.grid,error_u,error_p,"color mapped error");
        if(null&&iss.Nullspace_Check(rhs)){
            OCTAVE_OUTPUT<T>("n.txt").Write("n",rhs);
            iss.Multiply(rhs,*vectors(0));
            LOG::cout<<"nullspace found: "<<sqrt(iss.Inner_Product(*vectors(0),*vectors(0)))<<std::endl;
            rhs*=1/rhs.Max_Abs();
            Dump_Vector2<T,TV>(iss,rhs,"extra null mode");}}
        
    if(dump_matrix) OCTAVE_OUTPUT<T>("M.txt").Write("M",iss,*vectors(0),*vectors(1));
}

//#################################################################################################################################################
// Integration Test ###############################################################################################################################
//#################################################################################################################################################

template<class TV>
void Integration_Test(int argc,char* argv[],PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef VECTOR<int,d> TV_INT;

    Get_Debug_Particles<TV>().debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);

    int test_number;
    if(parse_args.Num_Extra_Args()<1){LOG::cerr<<"Test number is required."<<std::endl; exit(-1);}
    if(!STRING_UTILITIES::String_To_Value(parse_args.Extra_Arg(0),test_number)) throw VALUE_ERROR("The argument is not an integer.");

    ANALYTIC_TEST<TV>* test=0;

    switch(test_number){
        case 0:{ // Linear flow [u,v]=[0,v(x)] on [0,1/3],[1/3,2/3],[2/3,1], no pressure
            struct ANALYTIC_TEST_0:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){}
                virtual TV u(const TV& X,bool inside)
                {return TV::Axis_Vector(1)*(inside?(2*X.x-m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X){return -m/(T)6+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV();}
                virtual TV interface(const TV& X)
                {return TV::Axis_Vector(1)*((X.x>0.5*m)?(T)(-1):(T)1)*(2*mu(1)+mu(0))/s;}
            };
            test=new ANALYTIC_TEST_0;
            break;}
        case 1:{ // Linear flow [u,v]=[0,v(x)] on [0,0.25],[0.25,0.75],[0.75,1], no pressure
            struct ANALYTIC_TEST_1:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){}
                virtual TV u(const TV& X,bool inside)
                {return TV::Axis_Vector(1)*(inside?(X.x-0.5*m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X){return -0.25*m+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV();}
                virtual TV interface(const TV& X)
                {return TV::Axis_Vector(1)*((X.x>0.5*m)?(T)(-1):(T)1)*mu.Sum()/s;}
            };
            test=new ANALYTIC_TEST_1;
            break;}
        case 2:{ // Poiseuille flow (parabolic velocity profile) [u,v]=[u,v(x)] on [0.25,0.75], 0 velocity outside, no pressure
            struct ANALYTIC_TEST_2:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){}
                virtual TV u(const TV& X,bool inside)
                {return TV::Axis_Vector(1)*inside*(sqr(0.25*m)-sqr(X.x-0.5*m))/(m*s);}
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X){return -0.25*m+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV::Axis_Vector(1)*(inside*2*mu(1)/(m*s));}
                virtual TV interface(const TV& X)
                {return (T)0.5*TV::Axis_Vector(1)*mu(1)/s;}
            };
            test=new ANALYTIC_TEST_2;
            break;}
        case 3:{ // Opposite Poiseuille flows on [0.25,0.75] and [0,0.25],[0.75,1], no pressure
            struct ANALYTIC_TEST_3:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){}
                virtual TV u(const TV& X,bool inside)
                {return TV::Axis_Vector(1)*(inside?(sqr(0.25*m)-sqr(X.x-0.5*m)):((X.x>0.5*m)?(sqr(X.x-m)-sqr(0.25*m)):(sqr(X.x)-sqr(0.25*m))))/(m*s);}
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X){return -0.25*m+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV::Axis_Vector(1)*(inside?mu(1):-mu(0))*2/(m*s);}
                virtual TV interface(const TV& X)
                {return (T)0.5*TV::Axis_Vector(1)*(mu(1)-mu(0))/s;}
            };
            test=new ANALYTIC_TEST_3;
            break;}
        case 4:{ // Circular flow: linear growth for r<R1, 0 for r>R2, blend for R1<r<R2, no pressure
            struct ANALYTIC_TEST_4:public ANALYTIC_TEST<TV>
            {
                T ra,rb,ra2,rb2,rb2mra2,r_avg,r_avg2;
                using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize()
                {
                    ra=0.1666*m;
                    rb=0.3333*m;
                    ra2=sqr(ra);
                    rb2=sqr(rb);
                    rb2mra2=rb2-ra2;
                    r_avg=(ra+rb)/2;
                    r_avg2=sqr(r_avg);
                }
                virtual TV u(const TV& X,bool inside)
                {
                    T x=X.x-0.5*m;
                    T y=X.y-0.5*m;
                    T r2=VECTOR<T,2>(x,y).Magnitude_Squared();
                    TV velocity;
                    if(inside){
                        velocity.x=-y/s;
                        velocity.y=x/s;
                        velocity*=(rb2-r2)/(rb2mra2);}
                    else if(r2<r_avg2){
                        velocity.x=-y/s;
                        velocity.y=x/s;}
                    return velocity;
                }
                virtual T p(const TV& X){return T();}
                virtual T phi(const TV& X)
                {
                    T r=VECTOR<T,2>(X.x-0.5*m,X.y-0.5*m).Magnitude();
                    return (ra-rb)/2+abs(r-r_avg);
                }
                virtual TV body(const TV& X,bool inside)
                {
                    TV force;
                    if(!inside) return force;
                    force.x=-(X.y-0.5*m);
                    force.y=X.x-0.5*m;
                    force*=mu(1)*8/(rb2mra2)/s;
                    return force;
                }
                virtual TV interface(const TV& X)
                {
                    VECTOR<T,2> z(X.x-0.5*m,X.y-0.5*m);
                    T r=z.Normalize();
                    z=z.Rotate_Counterclockwise_90();
                    z*=mu(1)*2*sqr(r)/rb2mra2/s;
                    TV force;
                    force.x=z.x;
                    force.y=z.y;
                    if(r<r_avg) force*=-1;
                    return force;
                }
            };
            test=new ANALYTIC_TEST_4;
            break;}
        case 5:{ // Linear flow [u,v]=[0,v(x)] on [0,1/3],[1/3,2/3],[2/3,1], periodic pressure p(y)
            struct ANALYTIC_TEST_5:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){}
                virtual TV u(const TV& X,bool inside)
                {return TV::Axis_Vector(1)*(inside?(2*X.x-m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return (phi(X)<0)*sin(2*M_PI*X.y/m)*kg/(sqr(s)*(d==3?m:1));}
                virtual T phi(const TV& X){return -m/(T)6+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV::Axis_Vector(1)*2*M_PI*cos(2*M_PI*X.y/m)*kg/(sqr(s)*(d==3?sqr(m):m))*inside;}
                virtual TV interface(const TV& X)
                {return (TV::Axis_Vector(1)*(2*mu(1)+mu(0))/s-TV::Axis_Vector(0)*sin(2*M_PI*X.y/m)*kg/(sqr(s)*(d==3?m:1)))*((X.x>0.5*m)?(T)(-1):(T)1);}
            };
            test=new ANALYTIC_TEST_5;
            break;}
        case 6:{ // Linear flow [u,v]=[0,v(x)] on [0,1/3],[1/3,2/3],[2/3,1], linear pressure p(x) inside
            struct ANALYTIC_TEST_6:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){}
                virtual TV u(const TV& X,bool inside)
                {return TV::Axis_Vector(1)*(inside?(2*X.x-m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return (phi(X)<0)*(X.x-0.5*m)*kg/(sqr(s)*(d==3?sqr(m):m));}
                virtual T phi(const TV& X){return -m/(T)6+abs(X.x-0.5*m);}
                virtual TV body(const TV& X,bool inside){return TV::Axis_Vector(0)*kg/(sqr(s)*(d==3?sqr(m):m))*inside;}
                virtual TV interface(const TV& X)
                {return (TV::Axis_Vector(1)*(2*mu(1)+mu(0))/s-TV::Axis_Vector(0)*(X.x-0.5*m)*kg/(sqr(s)*(d==3?sqr(m):m)))*((X.x>0.5*m)?(T)(-1):(T)1);}
            };
            test=new ANALYTIC_TEST_6;
            break;}
        case 7:{ // No interface
            struct ANALYTIC_TEST_7:public ANALYTIC_TEST<TV>
            {
                virtual void Initialize(){}
                virtual TV u(const TV& X,bool inside){return TV();}
                virtual T p(const TV& X){return 0;}
                virtual T phi(const TV& X){return (-1);}
                virtual TV body(const TV& X,bool inside){return TV();}
                virtual TV interface(const TV& X){return TV();}
            };
            test=new ANALYTIC_TEST_7;
            break;}
        case 8:{ // Linear divergence field for r<R, 0 for r>R
            struct ANALYTIC_TEST_8:public ANALYTIC_TEST<TV>
            {
                T r;
                using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){r=m/3.0;}
                virtual TV u(const TV& X,bool inside){return (X-0.5*m)*inside;}
                virtual T p(const TV& X){return 0;}
                virtual T phi(const TV& X){return (X-0.5*m).Magnitude()-r;}
                virtual TV body(const TV& X,bool inside){return TV();}
                virtual TV interface(const TV& X){return -(X-0.5*m).Normalized()*2*mu(1);}
            };
            test=new ANALYTIC_TEST_8;
            break;}
        case 9:{ // Linear divergence field for r<R + quadratic pressure, 0 for r>R
            struct ANALYTIC_TEST_9:public ANALYTIC_TEST<TV>
            {
                T r;
                using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){r=m/M_PI;}
                virtual TV u(const TV& X,bool inside){return (X-0.5*m)*inside;}
                virtual T p(const TV& X){return (phi(X)<0)*(X-0.5*m).Magnitude_Squared();}
                virtual T phi(const TV& X){return (X-0.5*m).Magnitude()-r;}
                virtual TV body(const TV& X,bool inside){return (X-0.5*m)*2*inside;}
                virtual TV interface(const TV& X){return (X-0.5*m).Normalized()*((X-0.5*m).Magnitude_Squared()-2*mu(1));}
            };
            test=new ANALYTIC_TEST_9;
            break;}
        case 10:{
            struct ANALYTIC_TEST_10:public ANALYTIC_TEST<TV>
            {
                T r,m2,m4,u_term,p_term;
                using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){r=m/M_PI;m2=sqr(m);m4=sqr(m2);u_term=1;p_term=1;}
                virtual TV u(const TV& X,bool inside){TV x=X-0.5*m; return x*u_term*exp(-x.Magnitude_Squared()/m2)*inside;}
                virtual T p(const TV& X){return (phi(X)<0)*p_term*sin((X-0.5*m).Magnitude_Squared()/m2);}
                virtual T phi(const TV& X){return (X-0.5*m).Magnitude()-r;}
                virtual TV body(const TV& X,bool inside)
                {
                    TV x=X-0.5*m;
                    T x2=x.Magnitude_Squared();
                    T x2m2=x2/m2;
                    return x*(2*p_term*cos(x2m2)/m2+u_term*mu(1)*((TV::m+2)*m2-2*x2)*4*exp(-x2m2)/m4)*inside;
                }
                virtual TV interface(const TV& X)
                {
                    TV x=X-0.5*m;
                    T x2m2=x.Magnitude_Squared()/m2;
                    return (p_term*sin(x2m2)+u_term*2*mu(1)*exp(-x2m2)*(2*x2m2-1))*x.Normalized();  
                }
            };
            test=new ANALYTIC_TEST_10;
            break;}
        case 11:{
            struct ANALYTIC_TEST_11:public ANALYTIC_TEST<TV>
            {
                T r;
                using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){r=1.0/3.0;}
                virtual TV u(const TV& X,bool inside){TV x=X-0.5; return x*exp(x.Magnitude_Squared())*inside;}
                virtual T p(const TV& X){return 0;}
                virtual T phi(const TV& X){return (X-0.5).Magnitude()-r;}
                virtual TV body(const TV& X,bool inside)
                {
                    TV x=X-0.5;
                    T x2=x.Magnitude_Squared();
                    return x*(-1)*4*(2+x2)*exp(x2)*inside;
                }
                virtual TV interface(const TV& X)
                {
                    TV x=X-0.5;
                    T x2=x.Magnitude_Squared();
                    return x.Normalized()*(-1)*2*exp(x2)*(2*x2+1);  
                }
            };
            test=new ANALYTIC_TEST_11;
            break;}
        default:{
        LOG::cerr<<"Unknown test number."<<std::endl; exit(-1); break;}}

    output_directory=parse_args.Get_String_Value("-o");
    test->m=parse_args.Get_Double_Value("-m");
    test->s=parse_args.Get_Double_Value("-sec");
    test->kg=parse_args.Get_Double_Value("-kg");
    test->mu(0)=parse_args.Get_Double_Value("-mu_o")*test->kg/(test->s*(d==3?test->m:1));
    test->mu(1)=parse_args.Get_Double_Value("-mu_i")*test->kg/(test->s*(d==3?test->m:1));
    int res=parse_args.Get_Integer_Value("-resolution");
    int cgf=parse_args.Get_Integer_Value("-cgf");
    int max_iter=parse_args.Get_Integer_Value("-max_iter");
    bool use_preconditioner=parse_args.Get_Option_Value("-use_preconditioner");
    bool null=parse_args.Get_Option_Value("-null");
    bool dump_matrix=parse_args.Get_Option_Value("-dump_matrix");
    bool debug_particles=parse_args.Get_Option_Value("-debug_particles");

    if(res%cgf) PHYSBAM_FATAL_ERROR("Resolution must be divisible by coarse grid factor.");

    test->Initialize();

    TV_INT counts=TV_INT()+res;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1)*test->m,true);
    GRID<TV> coarse_grid(grid.counts/cgf,grid.domain,true);
    ARRAY<T,TV_INT> phi(coarse_grid.Node_Indices());

    Global_Grid(&grid);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    Analytic_Test(grid,coarse_grid,*test,max_iter,use_preconditioner,null,dump_matrix,debug_particles);
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
    parse_args.Add_Double_Argument("-mu_i",1,"viscosity inside");
    parse_args.Add_Double_Argument("-mu_o",1,"viscosity outside");
    parse_args.Add_Double_Argument("-m",1,"meter scale");
    parse_args.Add_Double_Argument("-sec",1,"second scale");
    parse_args.Add_Double_Argument("-kg",1,"kilogram scale");
    parse_args.Add_Integer_Argument("-test",1,"test number");
    parse_args.Add_Integer_Argument("-resolution",4,"resolution");
    parse_args.Add_Integer_Argument("-cgf",2,"coarse grid factor");
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
