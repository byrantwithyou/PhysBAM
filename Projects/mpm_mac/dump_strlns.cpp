//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <limits>
#include <string>
using namespace PhysBAM;

template<typename TV>
struct FRAME_DATA
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    GRID<TV> grid;
    ARRAY<T,FACE_INDEX<TV::m> > velocity;
    ARRAY<T,TV_INT> s;
};

template<typename TV>
void Read_Output_Files(FRAME_DATA<TV>& fd,const std::string& input,int frame)
{
    Read_From_File(LOG::sprintf("%s/common/grid",input.c_str()),fd.grid);
    Read_From_File(LOG::sprintf("%s/%d/mac_velocities",input.c_str(),frame),fd.velocity);
}

template<class TV>
class STREAMLINE_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T> KV;
    mutable ARRAY<T> temp_vector;
    KRYLOV_VECTOR_WRAPPER<T,KV> null_phi;
    const FRAME_DATA<TV>& fd;
    
    STREAMLINE_SYSTEM(const FRAME_DATA<TV>& g);
    virtual ~STREAMLINE_SYSTEM();

    void Compute_Ones_Nullspace();
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result,bool transpose=false) const override;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const override;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const override;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const override;
    //void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const override;
    void Test() const;
};

template<class TV> void STREAMLINE_SYSTEM<TV>::
Test() const
{
    int m=fd.grid.Numbers_Of_Nodes().Product();
    ARRAY<T> x(m),y(m);
    RANDOM_NUMBERS<T> random;
    random.Fill_Uniform(x,-1,1);
    random.Fill_Uniform(y,-1,1);
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > v,z,w,u(y);
    z.v.Resize(m);
    w.v.Resize(m);
    v.v=x;
    Multiply(v,z);
    Multiply(u,w);
    LOG::printf("sym %P %P\n",Inner_Product(z,u),Inner_Product(v,w));
    LOG::printf("pd %P\n",Inner_Product(v,z));
}
template<class TV> STREAMLINE_SYSTEM<TV>::
STREAMLINE_SYSTEM(const FRAME_DATA<TV>& f)
    :KRYLOV_SYSTEM_BASE<T>(false,false),fd(f)
{
}
template<class TV> STREAMLINE_SYSTEM<TV>::
~STREAMLINE_SYSTEM()
{
}
template<class TV> void STREAMLINE_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
template<class TV> void STREAMLINE_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result,bool transpose) const
{
    const KRYLOV_VECTOR_WRAPPER<T,KV>& vx=dynamic_cast<const KRYLOV_VECTOR_WRAPPER<T,KV>&>(x);
    KRYLOV_VECTOR_WRAPPER<T,KV>& vresult=dynamic_cast<KRYLOV_VECTOR_WRAPPER<T,KV>&>(result);
    RANGE<TV_INT> node_indices=fd.grid.Node_Indices();
    T dx=fd.grid.dX.Max();
    for(RANGE_ITERATOR<TV::m> it(node_indices);it.Valid();it.Next()){
        int ci=fd.s.Standard_Index(it.index);
        vresult.v(ci)=0;
        for(int j=0;j<2*TV::m;j++){
            TV_INT n=GRID<TV>::Node_Neighbor(it.index,j);
            if(node_indices.Lazy_Inside_Half_Open(n)){
                int cn=fd.s.Standard_Index(n);
                vresult.v(ci)+=-vx.v(cn);
                vresult.v(ci)+=vx.v(ci);}}
        vresult.v(ci)/=dx*dx;}
}
template<class TV> double STREAMLINE_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const KRYLOV_VECTOR_WRAPPER<T,KV>& vx=dynamic_cast<const KRYLOV_VECTOR_WRAPPER<T,KV>&>(x);
    const KRYLOV_VECTOR_WRAPPER<T,KV>& vy=dynamic_cast<const KRYLOV_VECTOR_WRAPPER<T,KV>&>(y);
    double r=0;
#pragma omp parallel for reduction(+:r)
    for(int i=0;i<vx.v.m;i++)
        r+=vx.v(i)*vy.v(i);
    return r;
}
template<class TV> typename TV::SCALAR STREAMLINE_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    const KRYLOV_VECTOR_WRAPPER<T,KV>& vx=dynamic_cast<const KRYLOV_VECTOR_WRAPPER<T,KV>&>(x);
    T r=0;
#pragma omp parallel for reduction(max:r)
    for(int i=0;i<vx.v.m;i++)
        r=std::max(r,std::abs(vx.v(i)));
    return r;
}
template<class TV> void STREAMLINE_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project_Nullspace(x);
}
template<class TV> void STREAMLINE_SYSTEM<TV>::
Compute_Ones_Nullspace()
{
    null_phi.v.Resize(fd.grid.Numbers_Of_Nodes().Product(),no_init);
    null_phi.v.Fill(1/sqrt(null_phi.v.m));
}
template<class TV> void STREAMLINE_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    KRYLOV_VECTOR_WRAPPER<T,KV>& v=dynamic_cast<KRYLOV_VECTOR_WRAPPER<T,KV>&>(x);
    v.Copy(-Inner_Product(v,null_phi),null_phi,v);
}

template<typename TV>
VECTOR<typename TV::SCALAR,2> Compute_Streamline(FRAME_DATA<TV>& fd)
{
    PHYSBAM_ASSERT(TV::m==2);
    typedef typename TV::SCALAR T;
    typedef ARRAY<T> KV;
    fd.s.Resize(fd.grid.Node_Indices());
    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    KRYLOV_VECTOR_WRAPPER<T,KV> sol,rhs;
    sol.v.Resize(fd.grid.Numbers_Of_Nodes().Product());
    rhs.v.Resize(fd.grid.Numbers_Of_Nodes().Product());
    for(RANGE_ITERATOR<TV::m> it(fd.grid.Node_Indices());it.Valid();it.Next()){
        int ci=fd.s.Standard_Index(it.index);
        rhs.v(ci)=0;
        for(int a=0;a<TV::m;a++) for(int f=0;f<(1<<(TV::m-1));f++){
                FACE_INDEX<TV::m> face(a,GRID<TV>::Node_Face_Index(a,it.index,f));
                if(fd.grid.domain.Lazy_Inside(fd.grid.Face(face)))
                    rhs.v(ci)+=(-2*f+1)*(-2*a+1)*fd.velocity(face);}}
    STREAMLINE_SYSTEM<TV> sys(fd);
    sys.Compute_Ones_Nullspace();

    CONJUGATE_GRADIENT<T> cg;
    cg.finish_before_indefiniteness=true;
    cg.relative_tolerance=true;
    T tol=1e-6;
    int iterations=500;
    bool converged=cg.Solve(sys,sol,rhs,av,tol,0,iterations);
    if(!converged) LOG::printf("SOLVER DID NOT CONVERGE.\n");

    sys.Multiply(sol,*av(0));
    *av(0)-=rhs;
    T r=sys.Convergence_Norm(*av(0));
    LOG::printf("residual: %P\n",r);
    //sys.Test();

    T smin=std::numeric_limits<T>::max(),smax=std::numeric_limits<T>::min();
    for(RANGE_ITERATOR<TV::m> it(fd.grid.Node_Indices());it.Valid();it.Next()){
        fd.s(it.index)=sol.v(fd.s.Standard_Index(it.index));
        smin=min(smin,fd.s(it.index));
        smax=max(smax,fd.s(it.index));}
    LOG::printf("smin: %P\nsmax: %P\n",smin,smax);
    av.Delete_Pointers_And_Clean_Memory();
    return VECTOR<T,2>(smin,smax);
}

template<typename T>
struct STYLE
{
    T dir_scale=1,width_scale=0.3; // times dx
    VECTOR<T,3> contour_color=VECTOR<T,3>(0,1,1);
    VECTOR<T,3> arrow_color=VECTOR<T,3>(1,0,0);
};

template<typename T,typename TV>
void Dump_Streamline(FRAME_DATA<TV>& fd,T c,const STYLE<T>& style)
{
    typedef VECTOR<int,TV::m> TV_INT;
    static const int num_corners=1<<TV::m;
    T dx=fd.grid.dX.Max();
    auto l2g=[&fd](const auto& X,const TV_INT& cell)
    {
        return fd.grid.Center(cell)-fd.grid.dX/2+X*fd.grid.dX;
    };
    auto arrow=[&style,dx](const TV& x,const TV& v)
    {
        TV vn=v.Normalized();
        TV u;
        u(0)=-vn(1);
        u(1)=vn(0);
        vn*=style.dir_scale*dx;
        u*=style.width_scale*dx;
        return VECTOR<TV,3>(x+vn,x+u,x-u);
    };

    LINEAR_INTERPOLATION_MAC<TV,T> li(fd.grid);
    VECTOR<T,num_corners> phis;
    ARRAY<VECTOR<TV,TV::m> > surface;
    VECTOR<ARRAY<VECTOR<TV,TV::m> >,2*TV::m> boundary;
    VECTOR<VECTOR<ARRAY<VECTOR<TV,TV::m> >*,2>,2*TV::m> pboundary;
    for(int s=0;s<2*TV::m;s++) pboundary(s)(1)=&boundary(s);
    int cnt=0;
    for(CELL_ITERATOR<TV> it(fd.grid);it.Valid();it.Next()){
        MARCHING_CUBES<TV>::Compute_Phis_For_Cell(phis,fd.s,it.index);
        phis+=c;
        surface.Remove_All();
        for(int s=0;s<2*TV::m;s++) boundary(s).Remove_All();
        MARCHING_CUBES<TV>::Get_Elements_For_Cell(surface,pboundary,phis);
        for(const auto& t:surface){
            TV X=l2g(.5*(t(0)+t(1)),it.index);
            TV v=li.Clamped_To_Array(fd.velocity,X);
            if(cnt==0)
                Add_Debug_Object(arrow(X,v),style.arrow_color);
            Add_Debug_Object(l2g(t,it.index),style.contour_color);
            cnt=(cnt+1)%5;}}
}

template<typename T,typename TV>
void Dump_Streamlines(FRAME_DATA<TV>& fd,T smin,T smax,int n,const STYLE<T>& style)
{
    T ds=(smax-smin)/(n-1);
    for(int i=0;i<n;i++)
        Dump_Streamline(fd,smin+ds*i,style);
    Flush_Frame<TV>("surface");
}

int main(int argc, char* argv[])
{
    typedef double T;
    typedef VECTOR<T,2> TV;
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);

    PARSE_ARGS parse_args(argc,argv);
    int frame=0,ncontours=2;
    bool bw=false;
    STYLE<T> style;
    std::string input="input",output="output";
    parse_args.Extra(&input,"input","Input directory");
    parse_args.Add("-o",&output,"output","Output directory");
    parse_args.Add("-frame",&frame,"frame","Frame");
    parse_args.Add("-nc",&ncontours,"contour","Contours");
    parse_args.Add("-bw",&bw,"Blackwhite");
    parse_args.Add("-arrow_d",&style.dir_scale,"scale","Arrow scale");
    parse_args.Add("-arrow_w",&style.width_scale,"scale","Arrow width scale");
    parse_args.Parse();
    if(bw) style.contour_color=style.arrow_color=VECTOR<T,3>(0,0,0);

    FRAME_DATA<TV> fd;
    Read_Output_Files(fd,input,frame);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE(T(0)),fd.grid,output);
    Create_Directory(output+"/common");
    Flush_Frame(fd.velocity,"init");
    VECTOR<T,2> r=Compute_Streamline(fd);
    Dump_Streamlines(fd,r(0),r(1),ncontours,style);
    return 0;
}
