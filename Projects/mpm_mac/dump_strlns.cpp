//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/SORT.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/UNION_FIND.h>
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
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <limits>
#include <map>
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
void Read_Output_Files(FRAME_DATA<TV>& fd,const VIEWER_DIR& viewer_dir)
{
    Read_From_File(viewer_dir.output_directory+"/common/grid",fd.grid);
    Read_From_File(viewer_dir.current_directory+"/mac_velocities",fd.velocity);
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

template<typename TV>
struct STYLE
{
    typedef typename TV::SCALAR T;
    T dir_scale=1,width_scale=0.3,arrow_sep=8,radius=5; // times dx
    VECTOR<T,3> contour_color=VECTOR<T,3>(0,1,1);
    VECTOR<T,3> arrow_color=VECTOR<T,3>(1,0,0);
    RANGE<TV> domain=RANGE<TV>::Unit_Box();
};

template<typename T,typename TV,typename TV_INT>
void Cover(FRAME_DATA<TV>& fd,const TV_INT& cell,
    std::map<TV_INT,T,LEXICOGRAPHIC_COMPARE>& q,HASHTABLE<TV_INT>& covered,const STYLE<TV>& style)
{
    RANGE<TV_INT> box=RANGE<TV_INT>(cell).Thickened(style.radius).Intersect(fd.grid.Domain_Indices());
    for(RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next()){
        auto qitem=q.find(it.index);
        if(qitem!=q.end())
            q.erase(qitem);
        if(!covered.Contains(it.index))
            covered.Set(it.index);}
    for(RANGE_ITERATOR<TV::m> it(box,1,0,RI::ghost);it.Valid();it.Next()){
        if(fd.grid.Domain_Indices().Lazy_Inside_Half_Open(it.index)){
            if(covered.Contains(it.index) || q.find(it.index)!=q.end()) continue;
            TV_INT nodes[1<<TV::m];
            fd.grid.Nodes_In_Cell_From_Minimum_Corner_Node(it.index,nodes);
            T s=0;
            for(int j=0;j<(1<<TV::m);j++)
                s+=fd.s(nodes[j]);
            s/=(1<<TV::m);
            q.insert(std::make_pair(it.index,s));}}
}

template<typename T,typename TV,typename TV_INT>
void Dump_Streamline(FRAME_DATA<TV>& fd,T c,
    std::map<TV_INT,T,LEXICOGRAPHIC_COMPARE>& q,HASHTABLE<TV_INT>& covered,const STYLE<TV>& style)
{
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    T dx=fd.grid.dX.Max();
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
    auto finish_contour=[&fd,&style](TV x,T mag)
    {
        x=fd.grid.Clamp(x);
        LOG::printf("COORD %.16P %.16P %P\n",x(0),x(1),mag);
        LOG::printf("COORD \n",x(0),x(1));
    };
    auto emit_segment=[&fd,&style](const TV& x0,const TV& x1,T mag)
    {
        if(!style.domain.Lazy_Inside(x0) && !style.domain.Lazy_Inside(x1)){
            LOG::printf("COORD \n");
            return;}
        TV y0=fd.grid.Clamp(x0);
        TV y1=fd.grid.Clamp(x1);
        Add_Debug_Object(VECTOR<TV,2>(y0,y1),style.contour_color);
        LOG::printf("COORD %.16P %.16P %P\n",y0(0),y0(1),mag);
        if(y0!=x0 || y1!=x1)
            LOG::printf("COORD \n");
    };
    auto emit_arrow=[&fd,&style,arrow](const TV& x,const TV& d)
    {
        if(!style.domain.Lazy_Inside(x)) return;
        TV dn=d.Normalized();
        LOG::printf("ARROW %.16P %.16P %.16P %.16P\n",x(0),x(1),dn(0),dn(1));
        Add_Debug_Object(arrow(x,d),style.arrow_color);
    };

    T_SURFACE surface;
    MARCHING_CUBES<TV>::Create_Surface(surface,fd.grid,fd.s,c);
    ARRAY_VIEW<TV> X=surface.particles.X;
    ARRAY<VECTOR<int,2> >& E=surface.mesh.elements;
    HASHTABLE<int,int> next;
    UNION_FIND<int> contours(X.m+1);
    ARRAY<bool> indegree0(X.m,use_init,true);
    LINEAR_INTERPOLATION_MAC<TV,T> li(fd.grid);
    for(int i=0;i<E.m;i++){
        TV u=X(E(i)(1))-X(E(i)(0));
        TV y=0.5*(X(E(i)(0))+X(E(i)(1)));
        TV v=li.Clamped_To_Array(fd.velocity,y);
        if(u.Dot(v)<0) std::swap(E(i)(0),E(i)(1));
        PHYSBAM_ASSERT(!next.Contains(E(i)(0)));
        next.Set(E(i)(0),E(i)(1));
        contours.Union(E(i)(0),E(i)(1));
        indegree0(E(i)(1))=false;}
    HASHTABLE<int,int> head;
    for(int i=0;i<X.m;i++){
        if(indegree0(i))
            head.Set(contours.Find(i),i);}

    for(int i=0;i<E.m;i++){
        int start=head.Get_Default(contours.Find(E(i)(0)),E(i)(0));
        if(contours.Find(start)==contours.Find(X.m)) continue;
        contours.Union(start,X.m);
        int it=start,n=next.Get(start);
        T len=0;
        do{
            TV y=0.5*(X(it)+X(n));
            TV d=X(n)-X(it);
            T l=d.Normalize();
            if(len+l>=style.arrow_sep*dx){
                TV z=X(it)+d*(style.arrow_sep*dx-len);
                TV v=li.Clamped_To_Array(fd.velocity,z);
                len-=style.arrow_sep*dx;
                emit_arrow(z,v);}
            len+=l;
            Cover(fd,fd.grid.Clamp_To_Cell(y),q,covered,style);
            emit_segment(X(it),X(n),li.Clamped_To_Array(fd.velocity,X(it)).Magnitude());
            it=n;
            if(next.Contains(it))
                n=next.Get(it);
            else break;
        }while(it!=start);
        finish_contour(X(it),li.Clamped_To_Array(fd.velocity,X(it)).Magnitude());}
}

template<typename T,typename TV>
void Dump_Streamlines(FRAME_DATA<TV>& fd,T sc,const STYLE<TV>& style)
{
    typedef VECTOR<int,TV::m> TV_INT;
    std::map<TV_INT,T,LEXICOGRAPHIC_COMPARE> q; // cell index -> contour value
    HASHTABLE<TV_INT> covered;
    Dump_Streamline(fd,sc,q,covered,style);
    while(!q.empty()){
        T c=q.begin()->second;
        Dump_Streamline(fd,c,q,covered,style);}
    Flush_Frame("surface");
}

int main(int argc, char* argv[])
{
    typedef double T;
    typedef VECTOR<T,2> TV;
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);

    PARSE_ARGS parse_args(argc,argv);
    int frame=0;
    T seed_contour=0;
    STYLE<TV> style;
    VIEWER_DIR viewer_in("input");
    VIEWER_DIR viewer_out("output");
    parse_args.Extra(&viewer_in.output_directory,"input","Input directory");
    parse_args.Add("-o",&viewer_out.output_directory,"output","Output directory");
    parse_args.Add("-frame",&frame,"frame","Frame");
    parse_args.Add("-c",&seed_contour,"contour","Starting contour value");
    parse_args.Add("-r",&style.radius,"scale","Contour span radius");
    parse_args.Add("-arrow_d",&style.dir_scale,"scale","Arrow scale");
    parse_args.Add("-arrow_w",&style.width_scale,"scale","Arrow width scale");
    parse_args.Add("-arrow_sep",&style.arrow_sep,"scale","Arrow separation");
    parse_args.Parse();

    FRAME_DATA<TV> fd;
    viewer_in.Set(frame);
    Read_Output_Files(fd,viewer_in);
    VIEWER_OUTPUT vo(STREAM_TYPE(T(0)),viewer_out);
    Use_Debug_Particles<TV>();
    vo.Add_Common("grid",fd.grid);
    vo.Add("mac_velocities",fd.velocity);
    Flush_Frame("init");
    VECTOR<T,2> r=Compute_Streamline(fd);
    Dump_Streamlines(fd,seed_contour,style);
    return 0;
}
