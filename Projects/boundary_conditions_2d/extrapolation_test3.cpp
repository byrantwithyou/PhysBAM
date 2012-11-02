#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <iomanip>

using namespace PhysBAM;

template<class T> T f(const VECTOR<T,3>& r){return cos(r.x)*sin(r.y)*sin(r.z);}

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> TV_INT;

int main(int argc,char* argv[])
{
    LOG::cout<<std::setprecision(2*sizeof(T));

    int resolution=atoi(argv[1]);
    int iterations=atoi(argv[2]);
    int order=atoi(argv[3]);
    int ghost=atoi(argv[4]);

    GRID<TV> grid(TV_INT()+resolution,RANGE<TV>::Centered_Box()*(T)pi,true);
    ARRAY<T,FACE_INDEX<TV::m> > x(grid.Domain_Indices(ghost));
    ARRAY<bool,FACE_INDEX<TV::m> > inside(grid.Domain_Indices());
    ARRAY<T,TV_INT> p(grid.Domain_Indices());
    LEVELSET<TV> phi(grid,p);

    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        TV X(it.Location());
        T r=X.Magnitude()-sqrt((T)2.0);
        p(it.index)=r;
        if(r<=0){
            x(it.Full_Index())=f(it.Location());
            inside(it.Full_Index())=true;}}

    EXTRAPOLATION_HIGHER_ORDER<TV,T>::Extrapolate_Face(grid,phi,inside,ghost,x,iterations,order,ghost);

    T m=0;
    for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid);it.Valid();it.Next()){
        TV X(it.Location()),Y=X;
        T r=Y.Normalize(),c=sqrt((T)2.0);
        T e=x(it.Full_Index())-f(X);
        if(r-c<grid.dX.Max()*ghost) m=max(abs(e),m);}
    LOG::cout<<"L-inf "<<m<<std::endl;

    return 0;
}
