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
#include <PhysBAM_Geometry/Finite_Elements/BOUNDARY_CONDITIONS_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/CELL_MANAGER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/VOLUME_FORCE_COLOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;

typedef float RW;
std::string output_directory="output";

typedef VECTOR<double,3> TV3;
TV3 color_map[4]={TV3(0,0.7,0),TV3(0.8,0.8,0),TV3(0,0.4,1),TV3(0.8,0.2,0)};

template<class TV>
GRID<TV>* Global_Grid(GRID<TV>* grid_in=0)
{
    static GRID<TV>* grid=0;
    GRID<TV>* old_grid=grid;
    if(grid_in) grid=grid_in;
    return old_grid;
}

template<class TV> struct ANALYTIC_TEST;


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
void Dump_Interface(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss)
{
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::CELL_ELEMENTS CELL_ELEMENTS;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::T_FACE T_FACE;
    typedef VECTOR<int,TV::m> TV_INT;

    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(iss.cdi->index_to_cell_elements);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<INTERFACE_ELEMENT>& interface_elements=cell_elements.interface;
        for(int i=0;i<interface_elements.m;i++){
            const INTERFACE_ELEMENT& V=interface_elements(i);
            if(V.color_pair.y>=0){
                if(V.color_pair.y>=0) Add_Debug_Object(V.face.X-V.face.Normal()*(T).003*iss.grid.dX.Min(),color_map[V.color_pair.y]);
                if("Alexey was here") Add_Debug_Object(V.face.X+V.face.Normal()*(T).003*iss.grid.dX.Min(),color_map[V.color_pair.x]);}
            else if(V.color_pair.x>=0) Add_Debug_Object(V.face.X-V.face.Normal()*(T).003*iss.grid.dX.Min(),color_map[V.color_pair.x]);}}
}

template<class T,class TV>
void Dump_System(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss,ANALYTIC_TEST<TV>& at)
{
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::CELL_ELEMENTS CELL_ELEMENTS;
    typedef typename CELL_DOMAIN_INTERFACE_COLOR<TV>::T_FACE T_FACE;
    typedef VECTOR<int,TV::m> TV_INT;

    Dump_Interface<T,TV>(iss);
    Flush_Frame<T,TV>("interface");
    const char* names[3]={"(neumann)","(dirichlet)","(slip)"};
    char buff[100];
    for(int c=0;c<iss.cdi->colors;c++){
        Dump_Interface<T,TV>(iss);
        sprintf(buff,"dofs %i %s",c,c>=0?"":names[-c-1]);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(iss.grid);it.Valid();it.Next()){
            int index=iss.cm_u(it.Axis())->Get_Index(it.index,c);
            if(index>=0){
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV::Axis_Vector(it.Axis()));}}
        Flush_Frame<T,TV>(buff);}

    for(int c=0;c<iss.cdi->colors;c++){
        Dump_Interface<T,TV>(iss);
        sprintf(buff,"dofs p %i %s",c,c>=0?"":names[-c-1]);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(iss.grid);it.Valid();it.Next()){
            int index=iss.cm_p->Get_Index(it.index,c);
            if(index>=0) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));}
        Flush_Frame<T,TV>(buff);}

    Dump_Interface<T,TV>(iss);
    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(iss.cdi->index_to_cell_elements);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<INTERFACE_ELEMENT>& interface_elements=cell_elements.interface;
        for(int i=0;i<interface_elements.m;i++){
            const INTERFACE_ELEMENT& V=interface_elements(i);
            if(V.color_pair.x<0) continue;
            Add_Debug_Particle(V.face.Center(),VECTOR<T,3>(0,0.1,0.5));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,at.j_surface(V.face.Center(),V.color_pair.x,V.color_pair.y));}}
    Flush_Frame<T,TV>("analytic interfacial forces");

    Dump_Interface<T,TV>(iss);
    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(iss.cdi->index_to_cell_elements);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<INTERFACE_ELEMENT>& interface_elements=cell_elements.interface;
        for(int i=0;i<interface_elements.m;i++){
            const INTERFACE_ELEMENT& V=interface_elements(i);
            if(V.color_pair.x!=-2 && V.color_pair.x!=-3) continue;
            Add_Debug_Particle(V.face.Center(),VECTOR<T,3>(0,0.1,0.5));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,at.d_surface(V.face.Center(),V.color_pair.x,V.color_pair.y));}}
    Flush_Frame<T,TV>("analytic dirichlet conditions");

    Dump_Interface<T,TV>(iss);
    for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(iss.cdi->index_to_cell_elements);it.Valid();it.Next()){
        const CELL_ELEMENTS& cell_elements=it.Data();
        const ARRAY<INTERFACE_ELEMENT>& interface_elements=cell_elements.interface;
        for(int i=0;i<interface_elements.m;i++){
            const INTERFACE_ELEMENT& V=interface_elements(i);
            if(V.color_pair.x!=-1 && V.color_pair.x!=-3) continue;
            Add_Debug_Particle(V.face.Center(),VECTOR<T,3>(0,0.1,0.5));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,at.n_surface(V.face.Center(),V.color_pair.x,V.color_pair.y));}}
    Flush_Frame<T,TV>("analytic neumann forces");
}

template<class T,class TV>
void Dump_Vector(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss,const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& v,const char* title)
{
    char buff[100];
    for(int c=0;c<iss.cdi->colors;c++){
        Dump_Interface<T,TV>(iss);
        sprintf(buff,"%s %i",title,c);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(iss.grid);it.Valid();it.Next()){
            int index=iss.cm_u(it.Axis())->Get_Index(it.index,c);
            if(index>=0){
                Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV::Axis_Vector(it.Axis())*v.u(it.Axis())(c)(index));}}
        Flush_Frame<T,TV>(buff);}

    for(int c=0;c<iss.cdi->colors;c++){
        Dump_Interface<T,TV>(iss);
        sprintf(buff,"%s p %i",title,c);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(iss.grid);it.Valid();it.Next()){
            int k=iss.cm_p->Get_Index(it.index,c);
            if(k>=0){
                Add_Debug_Particle(it.Location(),v.p(c)(k)>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,v.p(c)(k));}}
        Flush_Frame<T,TV>(buff);}
}

template<class T,class TV>
void Dump_Vector2(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss,const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& v,const char* title)
{
    char buff[100];
    for(int c=0;c<iss.cdi->colors;c++){
        sprintf(buff,"%s %i",title,c);
        Dump_Interface<T,TV>(iss);
        ARRAY<T,FACE_INDEX<TV::m> > u(iss.grid);
        for(UNIFORM_GRID_ITERATOR_FACE<TV> it(iss.grid);it.Valid();it.Next()){
            int i=it.Axis();
            int k=iss.cm_u(i)->Get_Index(it.index,c);
            if(k>=0) u(it.Full_Index())=v.u(i)(c)(k);}

        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(iss.grid);it.Valid();it.Next()){
            int k=iss.cm_p->Get_Index(it.index,c);
            if(k>=0){
                Add_Debug_Particle(it.Location(),v.p(c)(k)>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(1,0,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,v.p(c)(k));}}

        Dump_Frame(u,buff);}
}

template<class T,class TV,class TV_INT>
void Dump_u_p(const INTERFACE_STOKES_SYSTEM_COLOR<TV>& iss,const ARRAY<T,FACE_INDEX<TV::m> >& u,const ARRAY<T,TV_INT>& p,const char* title)
{
    char buff[100];
    Dump_Interface<T,TV>(iss);
    sprintf(buff,"%s",title);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(iss.grid);it.Valid();it.Next()){
        Add_Debug_Particle(it.Location(),VECTOR<T,3>(0.25,0.25,0.25));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV::Axis_Vector(it.Axis())*u(it.Full_Index()));}
    Flush_Frame<T,TV>(buff);

    Dump_Interface<T,TV>(iss);
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
struct ANALYTIC_TEST:public BOUNDARY_CONDITIONS_COLOR<TV>,public VOLUME_FORCE_COLOR<TV>
{
    typedef typename TV::SCALAR T;
    T kg,m,s;

    bool wrap;
    ARRAY<T> mu;

    ANALYTIC_TEST(): kg(1),m(1),s(1),wrap(false) {}
    virtual ~ANALYTIC_TEST(){}

    virtual void Initialize()=0;
    virtual T phi_value(const TV& X)=0;
    virtual int phi_color(const TV& X)=0;
    virtual TV u(const TV& X,int color)=0;
    virtual T p(const TV& X)=0;
    virtual TV F(const TV& X,int color)=0;

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
    
    INTERFACE_STOKES_SYSTEM_COLOR<TV> iss(grid,phi_value,phi_color,false);
    iss.use_preconditioner=use_preconditioner;
//    iss.use_p_null_mode=true;
//    iss.use_u_null_mode=true;
    iss.Set_Matrix(at.mu,at.wrap,&at,0,0);

    printf("\n");
    for(int i=0;i<TV::m;i++){for(int c=0;c<iss.cdi->colors;c++) printf("%c%d [%i]\t","uvw"[i],c,iss.cm_u(i)->dofs(c));printf("\n");}
    for(int c=0;c<iss.cdi->colors;c++) printf("p%d [%i]\t",c,iss.cm_p->dofs(c));printf("\n");
    printf("qn [%i]\t",iss.cdi->constraint_base_n);
    printf("qt [%i] ",iss.cdi->constraint_base_t);
    printf("\n");

    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> rhs,sol;
    {
        ARRAY<ARRAY<T,FACE_INDEX<TV::m> > > u;

        u.Resize(iss.cdi->colors);
        
        for(int c=0;c<iss.cdi->colors;c++) u(c).Resize(grid);

        for(int c=0;c<iss.cdi->colors;c++)
            for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
                FACE_INDEX<TV::m> face(it.Full_Index()); 
                u(c)(face)=at.u(it.Location(),c)(face.axis);}
        
        iss.Set_RHS(rhs,&at,&u,false);
        iss.Resize_Vector(sol);
    }

    MINRES<T> mr;
    CONJUGATE_RESIDUAL<T> cr;
    KRYLOV_SOLVER<T>* solver=&mr;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;
    KRYLOV_VECTOR_BASE<T>& tr=*rhs.Clone_Default();
    tr=rhs;

    if(dump_matrix){
        KRYLOV_SOLVER<T>::Ensure_Size(vectors,rhs,2);
        OCTAVE_OUTPUT<T>("M.txt").Write("M",iss,*vectors(0),*vectors(1));
        OCTAVE_OUTPUT<T>("b.txt").Write("b",rhs);}

//    solver->print_residuals=true;
    solver->Solve(iss,sol,rhs,vectors,1e-10,0,max_iter);
    if(dump_matrix) OCTAVE_OUTPUT<T>("x.txt").Write("x",sol);

    iss.Multiply(sol,*vectors(0));
    *vectors(0)-=rhs;
    LOG::cout<<"Residual: "<<iss.Convergence_Norm(*vectors(0))<<std::endl;

    for(int i=0;i<iss.null_modes.m;i++){
        iss.Multiply(*iss.null_modes(i),*vectors(0));
        LOG::cout<<"null mode["<<i<<"] "<<iss.Convergence_Norm(*vectors(0))<<std::endl;}

    ARRAY<T,FACE_INDEX<TV::m> > exact_u,numer_u,error_u;
    ARRAY<T,TV_INT> exact_p,numer_p,error_p;

    numer_u.Resize(iss.grid);
    exact_u.Resize(iss.grid);
    error_u.Resize(iss.grid);
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        int i=it.Axis();
        int c=at.phi_color(it.Location());
        if(c<0) continue;
        int k=iss.cm_u(it.Axis())->Get_Index(it.index,c);
        assert(k>=0);
        numer_u(it.Full_Index())=sol.u(i)(c)(k);}

    numer_p.Resize(iss.grid.Domain_Indices());
    exact_p.Resize(iss.grid.Domain_Indices());
    error_p.Resize(iss.grid.Domain_Indices());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        if(c<0) continue;
        int k=iss.cm_p->Get_Index(it.index,c);
        assert(k>=0);
        numer_p(it.index)=sol.p(c)(k);}

    TV avg_u;
    TV_INT cnt_u;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        if(c<0) continue;
        FACE_INDEX<TV::m> face(it.Full_Index()); 
        exact_u(face)=at.u(it.Location())(face.axis);
        error_u(face)=numer_u(face)-exact_u(face);
//        LOG::cout<<numer_u(face)<<"   "<<exact_u(face)<<"   "<<error_u(face)<<"   "<<at.phi_value(it.Location())<<std::endl;
        avg_u(face.axis)+=error_u(face);
        cnt_u(face.axis)++;}
    avg_u/=(TV)cnt_u;

    TV error_u_linf,error_u_l2;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        if(c<0) continue;
        FACE_INDEX<TV::m> face(it.Full_Index()); 
        error_u(face)-=avg_u(face.axis);
        error_u_linf(face.axis)=max(error_u_linf(face.axis),abs(error_u(face)));
        error_u_l2(face.axis)+=sqr(error_u(face));}
    error_u_l2=sqrt(error_u_l2/(TV)cnt_u);

    LOG::cout<<iss.grid.counts<<" U error:   linf="<<error_u_linf<<"   l2="<<error_u_l2<<std::endl;

    T avg_p=0;
    int cnt_p=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        int c=at.phi_color(it.Location());
        if(c<0) continue;
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
    enum WORKAROUND{SLIP=-3,DIRICHLET=-2,NEUMANN=-1};

    Get_Debug_Particles<TV>().debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);
    Get_Debug_Particles<TV>().debug_particles.template Add_Array<TV>(ATTRIBUTE_ID_V);

    T opt_s=1,opt_m=1,opt_kg=1;
    int res=4,max_iter=1000000;
    bool use_preconditioner=false,null=false,dump_matrix=false,debug_particles=false,bc_n=false,bc_d=false,bc_s=false;
    int test_number;
    parse_args.Add("-o",&output_directory,"dir","output directory");
    parse_args.Add("-m",&opt_m,"scale","meter scale");
    parse_args.Add("-s",&opt_s,"scale","second scale");
    parse_args.Add("-kg",&opt_kg,"scale","kilogram scale");
    parse_args.Add("-resolution",&res,"resolution","resolution");
    parse_args.Add("-max_iter",&max_iter,"iterations","max number of interations");
    parse_args.Add("-use_preconditioner",&use_preconditioner,"use Jacobi preconditioner");
    parse_args.Add("-null",&null,"find extra null modes of the matrix");
    parse_args.Add("-dump_matrix",&dump_matrix,"dump system matrix");
    parse_args.Add("-debug_particles",&debug_particles,"dump debug particles");
    parse_args.Add("-bc_n",&bc_n,"use Neumann boundary conditions");
    parse_args.Add("-bc_d",&bc_d,"use Dirichlet boundary conditions");
    parse_args.Add("-bc_s",&bc_s,"use slip boundary conditions");
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Parse();
    PHYSBAM_ASSERT(bc_n+bc_d+bc_s<2);
    int bc_type=bc_n?NEUMANN:(bc_s?SLIP:DIRICHLET);

    ANALYTIC_TEST<TV>* test=0;

    switch(test_number){
        case 0:{ // Two colors, periodic. Linear flow [u,v]=[0,v(x)] on [0,1/3],[1/3,2/3],[2/3,1], no pressure
            struct ANALYTIC_TEST_0:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);}
                virtual T phi_value(const TV& X){return abs(-m/(T)6+abs(X.x-0.5*m));}
                virtual int phi_color(const TV& X){return (-m/(T)6+abs(X.x-0.5*m))<0;}
                virtual TV u(const TV& X,int color){return TV::Axis_Vector(1)*(color?(2*X.x-m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return T();}
                virtual TV F(const TV& X,int color){return TV();}
                virtual TV j_surface(const TV& X,int color0,int color1){return TV::Axis_Vector(1)*((X.x>0.5*m)?(T)(-1):(T)1)*(2*mu(1)+mu(0))/s;}
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_0;
            break;}
        case 1:{ // Two colors, periodic. Linear flow [u,v]=[0,v(x)] on [0,0.25],[0.25,0.75],[0.75,1], no pressure
            struct ANALYTIC_TEST_1:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);}
                virtual T phi_value(const TV& X){return abs(-0.25*m+abs(X.x-0.5*m));}
                virtual int phi_color(const TV& X){return (-0.25*m+abs(X.x-0.5*m))<0;}
                virtual TV u(const TV& X,int color){return TV::Axis_Vector(1)*(color?(X.x-0.5*m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return T();}
                virtual TV F(const TV& X,int color){return TV();}
                virtual TV j_surface(const TV& X,int color0,int color1){return TV::Axis_Vector(1)*((X.x>0.5*m)?(T)(-1):(T)1)*mu.Sum()/s;}
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_1;
            break;}
        case 2:{ // Two colors, periodic. Poiseuille flow (parabolic velocity profile) [u,v]=[u,v(x)] on [0.25,0.75], 0 velocity outside, no pressure
            struct ANALYTIC_TEST_2:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);}
                virtual T phi_value(const TV& X){return abs(-0.25*m+abs(X.x-0.5*m));}
                virtual int phi_color(const TV& X){return (-0.25*m+abs(X.x-0.5*m))<0;}
                virtual TV u(const TV& X,int color)
                {return TV::Axis_Vector(1)*color*(sqr(0.25*m)-sqr(X.x-0.5*m))/(m*s);}
                virtual T p(const TV& X){return T();}
                virtual TV F(const TV& X,int color){return TV::Axis_Vector(1)*(color*2*mu(1)/(m*s));}
                virtual TV j_surface(const TV& X,int color0,int color1){return (T)0.5*TV::Axis_Vector(1)*mu(1)/s;}
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_2;
            break;}
        case 3:{ // Two colors, periodic. Opposite Poiseuille flows on [0.25,0.75] and [0,0.25],[0.75,1], no pressure
            struct ANALYTIC_TEST_3:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);}
                virtual T phi_value(const TV& X){return abs(-0.25*m+abs(X.x-0.5*m));}
                virtual int phi_color(const TV& X){return (-0.25*m+abs(X.x-0.5*m))<0;}
                virtual TV u(const TV& X,int color)
                {return TV::Axis_Vector(1)*(color?(sqr(0.25*m)-sqr(X.x-0.5*m)):((X.x>0.5*m)?(sqr(X.x-m)-sqr(0.25*m)):(sqr(X.x)-sqr(0.25*m))))/(m*s);}
                virtual T p(const TV& X){return T();}
                virtual TV F(const TV& X,int color){return TV::Axis_Vector(1)*(color?mu(1):-mu(0))*2/(m*s);}
                virtual TV j_surface(const TV& X,int color0,int color1){return (T)0.5*TV::Axis_Vector(1)*(mu(1)-mu(0))/s;}
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_3;
            break;}
        case 4:{ // Two colors, periodic. Circular flow: linear growth for r<ra, 0 for r>rb, blend for ra<r<rb, no pressure
            struct ANALYTIC_TEST_4:public ANALYTIC_TEST<TV>
            {
                T ra,rb,ra2,rb2,rb2mra2,r_avg,r_avg2;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize()
                {
                    wrap=true;mu.Append(1);mu.Append(2);
                    ra=0.1666*m;rb=0.3333*m;ra2=sqr(ra);rb2=sqr(rb);
                    rb2mra2=rb2-ra2;r_avg=(ra+rb)/2;r_avg2=sqr(r_avg);
                }
                virtual T phi_value(const TV& X){return abs((ra-rb)/2+abs(VECTOR<T,2>(X.x-0.5*m,X.y-0.5*m).Magnitude()-r_avg));}
                virtual int phi_color(const TV& X){return ((ra-rb)/2+abs(VECTOR<T,2>(X.x-0.5*m,X.y-0.5*m).Magnitude()-r_avg))<0;}
                virtual TV u(const TV& X,int color)
                {
                    T x=X.x-0.5*m,y=X.y-0.5*m;
                    T r2=VECTOR<T,2>(x,y).Magnitude_Squared();
                    TV velocity;
                    if(color){
                        velocity.x=-y/s;
                        velocity.y=x/s;
                        velocity*=(rb2-r2)/(rb2mra2);}
                    else if(r2<r_avg2){
                        velocity.x=-y/s;
                        velocity.y=x/s;}
                    return velocity;
                }
                virtual T p(const TV& X){return T();}
                virtual TV F(const TV& X,int color)
                {
                    TV force;
                    if(!color) return force;
                    force.x=-(X.y-0.5*m);
                    force.y=X.x-0.5*m;
                    force*=mu(1)*8/(rb2mra2)/s;
                    return force;
                }
                virtual TV j_surface(const TV& X,int color0,int color1)
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
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_4;
            break;}
        case 5:{ // Two colors, periodic. Linear flow [u,v]=[0,v(x)] on [0,1/3],[1/3,2/3],[2/3,1], periodic pressure p(y)
            struct ANALYTIC_TEST_5:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);}
                virtual T phi_value(const TV& X){return abs(-m/(T)6+abs(X.x-0.5*m));}
                virtual int phi_color(const TV& X){return (-m/(T)6+abs(X.x-0.5*m))<0;}
                virtual TV u(const TV& X,int color){return TV::Axis_Vector(1)*(color?(2*X.x-m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return phi_color(X)*sin(2*M_PI*X.y/m)*kg/(sqr(s)*(TV::m==3?m:1));}
                virtual TV F(const TV& X,int color){return TV::Axis_Vector(1)*2*M_PI*cos(2*M_PI*X.y/m)*kg/(sqr(s)*(TV::m==3?sqr(m):m))*color;}
                virtual TV j_surface(const TV& X,int color0,int color1)
                {return (TV::Axis_Vector(1)*(2*mu(1)+mu(0))/s-TV::Axis_Vector(0)*sin(2*M_PI*X.y/m)*kg/(sqr(s)*(TV::m==3?m:1)))*((X.x>0.5*m)?(T)(-1):(T)1);}
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_5;
            break;}
        case 6:{ // Two colors, periodic. Linear flow [u,v]=[0,v(x)] on [0,1/3],[1/3,2/3],[2/3,1], linear pressure p(x) color
            struct ANALYTIC_TEST_6:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);}
                virtual T phi_value(const TV& X){return abs(-m/(T)6+abs(X.x-0.5*m));}
                virtual int phi_color(const TV& X){return (-m/(T)6+abs(X.x-0.5*m))<0;}
                virtual TV u(const TV& X,int color){return TV::Axis_Vector(1)*(color?(2*X.x-m):((X.x>0.5*m)?(m-X.x):(-X.x)))/s;}
                virtual T p(const TV& X){return phi_color(X)*(X.x-0.5*m)*kg/(sqr(s)*(TV::m==3?sqr(m):m));}
                virtual TV F(const TV& X,int color){return TV::Axis_Vector(0)*kg/(sqr(s)*(TV::m==3?sqr(m):m))*color;}
                virtual TV j_surface(const TV& X,int color0,int color1)
                {return (TV::Axis_Vector(1)*(2*mu(1)+mu(0))/s-TV::Axis_Vector(0)*(X.x-0.5*m)*kg/(sqr(s)*(TV::m==3?sqr(m):m)))*((X.x>0.5*m)?(T)(-1):(T)1);}
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_6;
            break;}
        case 7:{ // One color, periodic. No interface, no forces
            struct ANALYTIC_TEST_7:public ANALYTIC_TEST<TV>
            {
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);}
                virtual T phi_value(const TV& X){return 1;}
                virtual int phi_color(const TV& X){return 0;}
                virtual TV u(const TV& X,int color){return TV();}
                virtual T p(const TV& X){return 0;}
                virtual TV F(const TV& X,int color){return TV();}
                virtual TV j_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_7;
            break;}
        case 8:{ // Two colors, periodic. Linear divergence field for r<R, 0 for r>R, no pressure
            struct ANALYTIC_TEST_8:public ANALYTIC_TEST<TV>
            {
                T r;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);r=m/M_PI;}
                virtual T phi_value(const TV& X){return abs((X-0.5*m).Magnitude()-r);}
                virtual int phi_color(const TV& X){return ((X-0.5*m).Magnitude()-r)<0;}
                virtual TV u(const TV& X,int color){return (X-0.5*m)*color;}
                virtual T p(const TV& X){return 0;}
                virtual TV F(const TV& X,int color){return TV();}
                virtual TV j_surface(const TV& X,int color0,int color1){return -(X-0.5*m).Normalized()*2*mu(1);}
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_8;
            break;}
        case 9:{ // Two colors, periodic. Linear divergence field for r<R + quadratic pressure, 0 for r>R
            struct ANALYTIC_TEST_9:public ANALYTIC_TEST<TV>
            {
                T r;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);r=m/M_PI;}
                virtual T phi_value(const TV& X){return abs((X-0.5*m).Magnitude()-r);}
                virtual int phi_color(const TV& X){return ((X-0.5*m).Magnitude()-r)<0;}
                virtual TV u(const TV& X,int color){return (X-0.5*m)*color;}
                virtual T p(const TV& X){return phi_color(X)*(X-0.5*m).Magnitude_Squared();}
                virtual TV F(const TV& X,int color){return (X-0.5*m)*2*color;}
                virtual TV j_surface(const TV& X,int color0,int color1){return (X-0.5*m).Normalized()*((X-0.5*m).Magnitude_Squared()-2*mu(1));}
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_9;
            break;}
        case 10:{ // Two colors, periodic. Exponential divergence velocity field and radial sine pressure field for r<R, 0 for r>R
            struct ANALYTIC_TEST_10:public ANALYTIC_TEST<TV>
            {
                T r,m2,m4,u_term,p_term;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                virtual void Initialize(){wrap=true;mu.Append(1);mu.Append(2);r=m/M_PI;m2=sqr(m);m4=sqr(m2);u_term=1;p_term=1;}
                virtual T phi_value(const TV& X){return abs((X-0.5*m).Magnitude()-r);}
                virtual int phi_color(const TV& X){return ((X-0.5*m).Magnitude()-r)<0;}
                virtual TV u(const TV& X,int color){TV x=X-0.5*m; return x*u_term*exp(-x.Magnitude_Squared()/m2)*color;}
                virtual T p(const TV& X){return phi_color(X)*p_term*sin((X-0.5*m).Magnitude_Squared()/m2);}
                virtual TV F(const TV& X,int color)
                {
                    TV x=X-0.5*m;
                    T x2=x.Magnitude_Squared();
                    T x2m2=x2/m2;
                    return x*(2*p_term*cos(x2m2)/m2+u_term*mu(1)*((TV::m+2)*m2-2*x2)*4*exp(-x2m2)/m4)*color;
                }
                virtual TV j_surface(const TV& X,int color0,int color1)
                {
                    TV x=X-0.5*m;
                    T x2m2=x.Magnitude_Squared()/m2;
                    return (p_term*sin(x2m2)+u_term*2*mu(1)*exp(-x2m2)*(2*x2m2-1))*x.Normalized();  
                }
                virtual TV d_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_10;
            break;}
        case 11:{ // if r<R then u=x,v=-y.  Otherwise, boundary condition from commandline (default Dirichlet)
            struct ANALYTIC_TEST_11:public ANALYTIC_TEST<TV>
            {
                T r,m2,m4,u_term,p_term;
                int bc_type;
                using ANALYTIC_TEST<TV>::kg;using ANALYTIC_TEST<TV>::m;using ANALYTIC_TEST<TV>::s;using ANALYTIC_TEST<TV>::wrap;using ANALYTIC_TEST<TV>::mu;
                ANALYTIC_TEST_11(int bc):bc_type(bc){}
                virtual void Initialize(){wrap=true;mu.Append(1);r=m/pi;m2=sqr(m);m4=sqr(m2);u_term=1;p_term=1;}
                virtual T phi_value(const TV& X){return abs((X-0.5*m).Magnitude()-r);}
                virtual int phi_color(const TV& X){return ((X-0.5*m).Magnitude()-r)<0?0:bc_type;}
                virtual TV u(const TV& X,int color){TV x=X-0.5*m;TV u;u(0)=x.x;u(1)=-x.y;return u;}
                virtual T p(const TV& X){return 0;}
                virtual TV F(const TV& X,int color){return TV();}
                virtual TV j_surface(const TV& X,int color0,int color1){return TV();}
                virtual TV n_surface(const TV& X,int color0,int color1){TV x=X-0.5*m;x.Normalize();x(1)*=-1;return x*2;}
                virtual TV d_surface(const TV& X,int color0,int color1){return u(X,max(color0,color1));}
                virtual TV u_jump(const TV& X,int color0,int color1){return TV();}
            };
            test=new ANALYTIC_TEST_11(bc_type);
            break;}
        default:{
        LOG::cerr<<"Unknown test number."<<std::endl; exit(-1); break;}}

    test->m=opt_m;
    test->s=opt_s;
    test->kg=opt_kg;
    test->Initialize();

    TV_INT counts=TV_INT()+res;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1)*test->m,true);

    Global_Grid(&grid);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    Analytic_Test(grid,*test,max_iter,use_preconditioner,null,dump_matrix,debug_particles);
    LOG::Finish_Logging();
    delete test;
}

//#################################################################################################################################################
// Main ###########################################################################################################################################
//#################################################################################################################################################

int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool opt_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&opt_3d,"use 3D");
    parse_args.Parse(true);

    if(opt_3d)
        Integration_Test<VECTOR<double,3> >(argc,argv,parse_args);
    else
        Integration_Test<VECTOR<double,2> >(argc,argv,parse_args);

    return 0;
}
