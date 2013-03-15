#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <iomanip>

using namespace PhysBAM;

typedef float RW;
typedef double T;

T f(const VECTOR<T,1>& X) {return cos(X.x);}
T f(const VECTOR<T,2>& X) {return cos(X.x)*sin(X.y);};
T f(const VECTOR<T,3>& X) {return cos(X.x)*sin(X.y)*sin(X.z);}
template<class TV> T phi_f(const TV& X) {return X.Magnitude()-sqrt((T)2.0);};
//template<class TV> T phi_f(const TV& X) {return X.Max_Abs()-sqrt((T)2.0);};

template<class TV>
void Test(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    LOG::cout<<std::setprecision(2*sizeof(T));

    int resolution=32;
    int iterations=10;
    int order=3;
    int ghost=3;
    parse_args.Add("-resolution",&resolution,"res","resolution");
    parse_args.Add("-iterations",&iterations,"iter","iterations");
    parse_args.Add("-order",&order,"iter","order of accuracy");
    parse_args.Add("-ghost",&ghost,"ghost","number of layers to fill");
    parse_args.Parse();

    GRID<TV> grid(TV_INT()+resolution,RANGE<TV>::Centered_Box()*(T)pi,true);
    ARRAY<T,TV_INT> x(grid.Domain_Indices(ghost));
    ARRAY<bool,TV_INT> inside(grid.Domain_Indices());
    ARRAY<T,TV_INT> p(grid.Domain_Indices());
    LEVELSET<TV> phi(grid,p);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),grid,"output");
    vo.debug_particles.edge_separation=(T).02;

    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        TV X=it.Location();
        T r=phi_f(X);
        p(it.index)=r;
        if(r<=0){
            x(it.index)=f(X);
            inside(it.index)=true;}}

    EXTRAPOLATION_HIGHER_ORDER<TV,T>(grid,phi,iterations,order,ghost).Extrapolate_Cell(inside,x);

    T m=0;
    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        TV X(it.Location()),Y=X;
        T r=phi_f(X);
        T e=x(it.index)-f(X);
        if(r<grid.dX.Max()*ghost) m=max(abs(e),m);}
    LOG::cout<<"L-inf "<<m<<std::endl;
    vo.Flush_Frame("Finish");
}

int main(int argc,char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    int dim=2;
    parse_args.Add("-dim",&dim,"dimension","run in 1D / 2D / 3D");
    parse_args.Parse(true);
    if(dim==1) Test<VECTOR<T,1> >(parse_args);
    else if(dim==2) Test<VECTOR<T,2> >(parse_args);
    else if(dim==3) Test<VECTOR<T,3> >(parse_args);
    else PHYSBAM_FATAL_ERROR();

    return 0;
}
